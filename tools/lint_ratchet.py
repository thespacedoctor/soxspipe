#!/usr/bin/env python
"""*fail a change when a ruff finding lands on a line that change touched*

Ruff has no baseline or ratchet feature, so the gate is built here. It mirrors the
`diff-cover --fail-under=80` precedent already used for coverage: the roughly 1,380
pre-existing findings in `soxspipe/` are ignored entirely, and only findings on added
or modified lines fail.

Two modes:

```bash
python tools/lint_ratchet.py --compare-branch origin/develop
python tools/lint_ratchet.py --staged
```

The first compares against the merge base with the given branch, which is what the CI
job runs on a pull request. The second compares the index against `HEAD`, which is what
the pre-commit hook runs locally. Note that `--staged` picks the *lines* from the index
but hands ruff the files on disk, so an unstaged edit is linted along with the staged
one. That is what pre-commit wants, since it stashes unstaged changes before running a
hook, but a direct `--staged` run in a dirty tree judges the worktree, not the commit.

Three rules are exempt from the ratchet and gated absolutely, anywhere in the tree:
`E722` (bare except), `S110` (`try`/`except`/`pass`) and `S112` (`try`/`except`/`continue`).
The phase-1 exception sweep removed every one of them from the package, so there is no
debt left for a ratchet to tolerate, and a new one is a silent swallow rather than a
style lapse. `SIM105` stays ratcheted: preferring `contextlib.suppress` is a style call,
and a handler whose `try` body is the logger itself has a reason to be written the long way.

A finding is matched on the line ruff anchors it to, not on the whole span it covers.
A multi-line finding such as `PLR0915` is anchored at the `def` line, so adding a
statement to a function that is already too long does not fail the gate, while writing
a new function that is too long does. The alternative, intersecting the whole span,
would let any edit inside a long function resurrect a pre-existing finding, which is
the thing the ratchet exists to prevent. Whole-file counts are printed as context so
the debt stays visible without being gated.

**Usage:**

```bash
python tools/lint_ratchet.py --compare-branch origin/develop
```
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import NamedTuple

# THE +++ LINE OF A UNIFIED DIFF: +++ b/soxspipe/recipes/soxs_mbias.py
TARGET_FILE_PATTERN = re.compile(r"^\+\+\+ (?:b/)?(.*)$")

# THE HUNK HEADER OF A UNIFIED DIFF: @@ -12,3 +14,5 @@
HUNK_HEADER_PATTERN = re.compile(r"^@@ -\d+(?:,\d+)? \+(\d+)(?:,(\d+))? @@")

# THE +++ PATH GIT USES WHEN A FILE WAS DELETED, WHICH HAS NO ADDED LINES TO LINT
DEV_NULL = "/dev/null"

# THE EXTENSIONS RUFF LINTS. NAMING A PATH OUTSIDE THIS SET IS NOT SAFE: RUFF PARSES AN
# EXPLICITLY NAMED .yaml FILE AS PYTHON AND REPORTS THIRTEEN SYNTAX ERRORS FOR IT, SO
# THE DIFF IS NARROWED HERE RATHER THAN LEFT FOR --force-exclude TO SORT OUT.
LINTABLE_PATHSPECS = ("*.py", "*.pyi", "*.pyw", "*.ipynb")

# THE RULES NO CHANGE MAY REINTRODUCE ANYWHERE, NOT EVEN ON A LINE IT DID NOT TOUCH
HARD_RULE_CODES = ("E722", "S110", "S112")

# THE HARD RULES ARE CHECKED FROM THE REPOSITORY ROOT, SO A PYTHON FILE OUTSIDE soxspipe,
# tests AND tools — docs/source/conf.py TODAY, ANYTHING ADDED TOMORROW — IS GATED TOO.
# RUFF'S OWN DEFAULT EXCLUDES KEEP build/, dist/ AND THE LIKE OUT OF IT.
HARD_RULE_PATHS = (".",)

EXIT_CLEAN = 0
EXIT_FINDINGS = 1
EXIT_TOOL_ERROR = 2


class Finding(NamedTuple):
    """*one ruff finding, anchored to the line ruff reported it on*"""

    path: str
    line: int
    code: str
    message: str

    def __str__(self) -> str:
        return f"{self.path}:{self.line}: {self.code} {self.message}"


def parse_changed_lines(diffText: str) -> dict[str, set[int]]:
    """*the line numbers each file gained or had modified, read from a unified diff*

    The diff must be generated with `-U0`, so every line inside a hunk is an added line
    and the hunk header alone gives the range. A deleted file contributes nothing. A
    renamed file contributes under its new path, since that is the path ruff will report.

    **Key Arguments:**

    - ``diffText`` -- the output of `git diff -U0`

    **Return:**

    - ``changedLines`` -- new-side line numbers, keyed by repository-relative path
    """
    changedLines: dict[str, set[int]] = {}
    currentPath: str | None = None

    for line in diffText.splitlines():
        targetFile = TARGET_FILE_PATTERN.match(line)
        if targetFile:
            path = targetFile.group(1)
            currentPath = None if path == DEV_NULL else path
            continue

        hunkHeader = HUNK_HEADER_PATTERN.match(line)
        if not hunkHeader or currentPath is None:
            continue

        start = int(hunkHeader.group(1))
        # AN ABSENT LENGTH MEANS ONE LINE; A LENGTH OF 0 MEANS A PURE DELETION
        length = int(hunkHeader.group(2)) if hunkHeader.group(2) is not None else 1
        changedLines.setdefault(currentPath, set()).update(range(start, start + length))

    return changedLines


def parse_ruff_findings(ruffJson: str, repoRoot: Path) -> list[Finding]:
    """*the findings in a ruff JSON report, with paths made repository-relative*

    **Key Arguments:**

    - ``ruffJson`` -- the output of `ruff check --output-format json`
    - ``repoRoot`` -- the repository root the reported absolute paths are relative to

    **Return:**

    - ``findings`` -- the parsed findings, in the order ruff reported them
    """
    findings = []
    for entry in json.loads(ruffJson or "[]"):
        path = Path(entry["filename"])
        try:
            relativePath = path.relative_to(repoRoot).as_posix()
        except ValueError:
            # A FINDING OUTSIDE THE REPOSITORY CANNOT MATCH A DIFF HUNK, BUT KEEP IT LEGIBLE
            relativePath = path.as_posix()

        findings.append(Finding(relativePath, entry["location"]["row"], entry["code"] or "", entry["message"]))

    return findings


def findings_on_changed_lines(findings: list[Finding], changedLines: dict[str, set[int]]) -> list[Finding]:
    """*the findings anchored to a line the change touched*

    **Key Arguments:**

    - ``findings`` -- every ruff finding in the tree
    - ``changedLines`` -- new-side line numbers, keyed by repository-relative path

    **Return:**

    - ``gated`` -- the findings that fail the ratchet, sorted by path and line
    """
    gated = [finding for finding in findings if finding.line in changedLines.get(finding.path, ())]

    return sorted(gated, key=lambda finding: (finding.path, finding.line, finding.code))


def hard_rule_findings(findings: list[Finding]) -> list[Finding]:
    """*the findings that break a hard rule, wherever in the tree they sit*

    **Key Arguments:**

    - ``findings`` -- the findings from a hard-rule ruff run

    **Return:**

    - ``hard`` -- the hard-rule findings, sorted by path and line
    """
    hard = [finding for finding in findings if finding.code in HARD_RULE_CODES]

    return sorted(hard, key=lambda finding: (finding.path, finding.line, finding.code))


def run_git_diff(compareBranch: str | None, repoRoot: Path) -> str:
    """*the unified diff of the change under test, at zero context*

    **Key Arguments:**

    - ``compareBranch`` -- the branch to compare against, or None to compare the index against HEAD
    - ``repoRoot`` -- the repository to run git in

    **Return:**

    - ``diffText`` -- the raw diff output
    """
    # THE TWO -c SETTINGS PIN THE DIFF FORMAT THIS MODULE PARSES, SO A CONTRIBUTOR'S
    # OWN git CONFIG CANNOT DROP THE b/ PREFIX OR QUOTE A PATH AND SILENTLY UNGATE A FILE
    command = ["git", "-c", "diff.noprefix=false", "-c", "core.quotePath=false"]
    command += ["diff", "-U0", "--no-color", "--find-renames", "--no-ext-diff"]
    command += ["--staged"] if compareBranch is None else [f"{compareBranch}...HEAD"]
    command += ["--", *LINTABLE_PATHSPECS]

    return _run(command, repoRoot)


def ruff_command() -> list[str]:
    """*the command that runs ruff, without depending on `PATH`*

    The hermetic test suite runs pytest through `env` with an explicit variable list that
    leaves `PATH` at its default, so the environment's `bin` directory is not on it. Running
    ruff as a module of the running interpreter sidesteps that, and pins the run to the
    version installed beside the interpreter rather than whichever one `PATH` happens to
    reach first. Only when ruff is not installed there does this fall back to `PATH`.

    **Return:**

    - ``command`` -- the command prefix to run ruff with
    """
    if importlib.util.find_spec("ruff") is not None:
        return [sys.executable, "-m", "ruff"]

    return ["ruff"]


def run_ruff(paths: list[str], repoRoot: Path) -> str:
    """*the ruff report for the given paths, as JSON*

    **Key Arguments:**

    - ``paths`` -- the repository-relative paths to check. An empty list checks nothing
    - ``repoRoot`` -- the repository to run ruff in

    **Return:**

    - ``ruffJson`` -- the raw JSON report, or an empty report when there was nothing to check
    """
    if not paths:
        return "[]"

    # THE -- STOPS A PATH THAT BEGINS WITH A DASH BEING READ AS AN OPTION
    command = [*ruff_command(), "check", "--output-format", "json", "--force-exclude", "--", *paths]

    # RUFF EXITS 1 WHEN IT FINDS SOMETHING, WHICH IS NOT AN ERROR HERE
    return _run(command, repoRoot, allowedStatuses=(0, 1))


def run_ruff_selected(codes: tuple[str, ...], paths: list[str], repoRoot: Path) -> str:
    """*the ruff report for the given rules over the given paths, as JSON*

    Unlike `run_ruff`, this narrows the rule set and widens the paths: the hard rules are
    judged over the whole repository, so a finding on an untouched line still fails. It also
    runs ruff `--isolated` and `--ignore-noqa`, so neither a `# noqa` comment nor a
    `per-file-ignores` entry can silence a hard rule.

    **Key Arguments:**

    - ``codes`` -- the ruff rule codes to select
    - ``paths`` -- the repository-relative paths to check
    - ``repoRoot`` -- the repository to run ruff in

    **Return:**

    - ``ruffJson`` -- the raw JSON report
    """
    # --isolated DROPS THE PROJECT CONFIGURATION, SO A per-file-ignores ENTRY ADDED LATER
    # CANNOT QUIETLY REOPEN THE HOLE, AND --ignore-noqa STOPS A # noqa COMMENT DOING THE SAME.
    # A HARD RULE THAT CAN BE SWITCHED OFF IN THE FILE THAT BREAKS IT IS NOT A HARD RULE.
    command = [*ruff_command(), "check", "--isolated", "--ignore-noqa", "--select", ",".join(codes)]
    command += ["--output-format", "json", "--force-exclude"]
    # THE -- STOPS A PATH THAT BEGINS WITH A DASH BEING READ AS AN OPTION
    command += ["--", *paths]

    # RUFF EXITS 1 WHEN IT FINDS SOMETHING, WHICH IS NOT AN ERROR HERE
    return _run(command, repoRoot, allowedStatuses=(0, 1))


def _run(command: list[str], repoRoot: Path, allowedStatuses: tuple[int, ...] = (0,)) -> str:
    """*run a command in the repository and return its standard output*

    **Key Arguments:**

    - ``command`` -- the command and its arguments
    - ``repoRoot`` -- the working directory to run in
    - ``allowedStatuses`` -- the exit statuses that are not failures. Default *(0,)*

    **Return:**

    - ``output`` -- the command's standard output
    """
    try:
        # NO SHELL, A FIXED EXECUTABLE, AND EVERY PATH ARGUMENT PLACED AFTER A -- SEPARATOR
        completed = subprocess.run(command, cwd=repoRoot, capture_output=True, text=True, check=False)  # noqa: S603
    except OSError as error:
        raise RatchetError(f"could not run {_command_name(command)}: {error}") from error

    if completed.returncode not in allowedStatuses:
        failure = f"{' '.join(command)} failed with status {completed.returncode}"
        raise RatchetError(f"{failure}:\n{completed.stderr.strip()}")

    return completed.stdout


def _command_name(command: list[str]) -> str:
    """*the name to report a command by in an error message*

    **Key Arguments:**

    - ``command`` -- the command and its arguments

    **Return:**

    - ``name`` -- the tool's own name, not the interpreter that runs it as a module
    """
    # A COMMAND RUN AS python -m <tool> IS THE TOOL'S FAILURE TO REPORT, NOT THE INTERPRETER'S
    if len(command) > 2 and command[1] == "-m":
        return command[2]

    return command[0]


class RatchetError(RuntimeError):
    """*the ratchet could not do its job, which is not the same as finding lint*"""


def repository_root(start: Path) -> Path:
    """*the root of the git repository containing a path*

    **Key Arguments:**

    - ``start`` -- a path inside the repository

    **Return:**

    - ``root`` -- the absolute repository root
    """
    return Path(_run(["git", "rev-parse", "--show-toplevel"], start).strip())


def main(argv: list[str] | None = None) -> int:
    """*run the ratchet from the command line*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments, excluding the program name. Default *None*, i.e. `sys.argv[1:]`

    **Return:**

    - ``status`` -- 0 when no finding lands on a changed line, 1 when one does, 2 when the ratchet itself failed
    """
    parser = argparse.ArgumentParser(description="Fail a change when a ruff finding lands on a line it touched.")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--compare-branch", help="compare against the merge base with this branch, e.g. origin/develop")
    source.add_argument("--staged", action="store_true", help="compare the staged index against HEAD")
    arguments = parser.parse_args(argv)

    try:
        repoRoot = repository_root(Path.cwd())
        diffText = run_git_diff(None if arguments.staged else arguments.compare_branch, repoRoot)
        changedLines = parse_changed_lines(diffText)
        # A FILE WHOSE ONLY HUNK IS A DELETION HAS NO ADDED LINE TO GATE, SO RUFF SKIPS IT
        changedFiles = sorted(path for path, lines in changedLines.items() if lines)
        findings = parse_ruff_findings(run_ruff(changedFiles, repoRoot), repoRoot)
        hardReport = run_ruff_selected(HARD_RULE_CODES, list(HARD_RULE_PATHS), repoRoot)
        hard = hard_rule_findings(parse_ruff_findings(hardReport, repoRoot))
    except RatchetError as error:
        print(f"lint ratchet: {error}", file=sys.stderr)
        return EXIT_TOOL_ERROR
    except (json.JSONDecodeError, KeyError) as error:
        print(f"lint ratchet: could not read the ruff report: {error}", file=sys.stderr)
        return EXIT_TOOL_ERROR

    # A HARD-RULE FINDING ON A CHANGED LINE IS REPORTED ONCE, IN THE HARD-RULE SECTION,
    # RATHER THAN IN BOTH: THE PROJECT'S OWN RULE SET SELECTS THESE CODES TOO
    ratchetFindings = [finding for finding in findings if finding.code not in HARD_RULE_CODES]
    gated = findings_on_changed_lines(ratchetFindings, changedLines)
    _print_report(gated, ratchetFindings, changedLines)
    _print_hard_rule_report(hard)

    return EXIT_FINDINGS if gated or hard else EXIT_CLEAN


def _print_report(gated: list[Finding], findings: list[Finding], changedLines: dict[str, set[int]]) -> None:
    """*print the failing findings, then the whole-file context they sit in*

    **Key Arguments:**

    - ``gated`` -- the findings that fail the ratchet
    - ``findings`` -- every finding in the changed files
    - ``changedLines`` -- new-side line numbers, keyed by repository-relative path
    """
    changedLineCount = sum(len(lines) for lines in changedLines.values())
    print(f"lint ratchet: {len(changedLines)} changed Python files, {changedLineCount} changed lines")

    for finding in gated:
        print(finding)

    if gated:
        plural = "finding" if len(gated) == 1 else "findings"
        print(f"\n{len(gated)} {plural} on changed lines. Fix these before merging.")
    else:
        print("\nNo findings on changed lines.")

    # THE PRE-EXISTING DEBT IS REPORTED BUT NEVER GATED, SO THE RATCHET ONLY EVER TIGHTENS
    preExisting = len(findings) - len(gated)
    print(f"{preExisting} pre-existing findings in the same files, not gated.")


def _print_hard_rule_report(hard: list[Finding]) -> None:
    """*print the hard-rule failures, which are gated wherever they sit*

    **Key Arguments:**

    - ``hard`` -- the findings that break a hard rule
    """
    codes = ", ".join(HARD_RULE_CODES)
    print(f"\nhard rules ({codes}), checked over {', '.join(HARD_RULE_PATHS)}:")

    for finding in hard:
        print(finding)

    if hard:
        plural = "finding" if len(hard) == 1 else "findings"
        print(f"\n{len(hard)} hard rule {plural}. These fail wherever they sit, not only on changed lines.")
    else:
        print("clean")


if __name__ == "__main__":
    sys.exit(main())
