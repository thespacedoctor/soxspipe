#!/usr/bin/env python
"""*survey one soxspipe module and report the numbers a refactor ticket opens with*

The per-module refactor template requires each module ticket to start from generated
numbers rather than hand-counted ones. This script produces them, as Markdown bullets
that paste into a Linear ticket body or a pull-request description without editing.

It reports the module's line count, its coverage at the branch point, its functions
over 50 lines and over 200 lines, its duplication hits, its bare excepts, its
remaining whole-file ruff findings, and its docstring drift count.

**Duplication hits** are counted under one definition, which is the definition the
shared helpers in `soxspipe/commonutils/toolkit.py` were built to absorb:

- `naive-utcnow` -- a call to `utcnow()`, which `toolkit.utcnow_string` replaces.
- `qc-row-builder` -- `self.qc = pd.concat([..., pd.DataFrame([{...}]), ...])`, a
  single-row append to the QC table, which `base_recipe.add_qc` replaces.
- `products-row-builder` -- the same shape assigned to `self.products`, which
  `base_recipe.add_product` replaces.
- `savefig` -- a call to `savefig()`, which `toolkit.save_qc_plot` replaces.

A concatenation that appends a whole table rather than one dict literal is **not** a
row builder and is not counted, because no helper replaces it.

The `utcnow` and `savefig` rules match on the called name and ignore the receiver, on
purpose: `datetime.utcnow()` and `datetime.datetime.utcnow()` are both naive, and
`plt.savefig()`, `fig.savefig()` and `ax.figure.savefig()` are all calls the plot
helper replaces. The cost is that an unrelated method of the same name would be
counted; no such method exists in this package.

"At the branch point" means the working tree the surveying agent is about to branch
from, which is the tree this script measures. It does not check out another revision.

**Usage:**

```bash
python tools/survey_module.py soxspipe/recipes/soxs_mbias.py
python tools/survey_module.py soxspipe/recipes/soxs_mbias.py --no-coverage
```

Exit status is 0 when the survey was produced and 2 when it could not be, matching
`tools/check_docstrings.py` and `tools/lint_ratchet.py`. The survey is a measurement,
not a gate, so it never exits 1.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import subprocess
import sys
import tempfile
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import NamedTuple

# THE HOUSE RULE THE SURVEY MEASURES AGAINST, AND THE THRESHOLD ABOVE WHICH THE
# TEMPLATE REQUIRES A WRITTEN REASON FOR LEAVING A FUNCTION UNSPLIT
OVERSIZED_FUNCTION_LINES = 50
UNSPLIT_REASON_LINES = 200

# THE DUPLICATION KINDS, IN THE ORDER THE REPORT LISTS THEM
DUPLICATION_KINDS = ("naive-utcnow", "qc-row-builder", "products-row-builder", "savefig")

DUPLICATION_LABELS = {
    "naive-utcnow": ("naive `utcnow` timestamp site", "naive `utcnow` timestamp sites"),
    "qc-row-builder": ("QC-table `pd.concat` row builder", "QC-table `pd.concat` row builders"),
    "products-row-builder": ("products-table `pd.concat` row builder", "products-table `pd.concat` row builders"),
    "savefig": ("`plt.savefig` call", "`plt.savefig` calls"),
}

# THE PACKAGE COVERAGE MEASURES. A SINGLE MODULE PATH PASSED TO --cov MEASURES NOTHING.
COVERAGE_PACKAGE = "soxspipe"

# THE TABLE ATTRIBUTES A SINGLE-ROW pd.concat APPENDS TO, AND THE KIND EACH REPORTS AS
ROW_BUILDER_ATTRIBUTES = {"qc": "qc-row-builder", "products": "products-row-builder"}

EXIT_CLEAN = 0
EXIT_TOOL_ERROR = 2


class SurveyError(RuntimeError):
    """*the survey could not be produced, which is not the same as finding nothing*"""


class FunctionLength(NamedTuple):
    """*one function that exceeds the 50-line house rule*"""

    name: str
    line: int
    length: int


class Hit(NamedTuple):
    """*one duplication hit, of one of the four counted kinds*"""

    kind: str
    line: int


class StaticSurvey(NamedTuple):
    """*everything the survey reads straight out of the module's source*"""

    path: Path
    lineCount: int
    functions: list[FunctionLength]
    hits: list[Hit]
    bareExcepts: list[int]


class Survey(NamedTuple):
    """*the complete survey of one module, ready to render*"""

    path: Path
    lineCount: int
    functions: list[FunctionLength]
    hits: list[Hit]
    bareExcepts: list[int]
    coveragePercent: float | None
    ruffCounts: dict[str, int]
    docstringFindings: int
    driftedFunctions: int


def static_survey(path: Path) -> StaticSurvey:
    """*read a module's source and measure everything that needs no other tool*

    **Key Arguments:**

    - ``path`` -- the path of the module to survey

    **Return:**

    - ``survey`` -- the line count, oversized functions, duplication hits and bare excepts
    """
    path = Path(path)
    try:
        source = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        raise SurveyError(f"could not read {path}: {error}") from error

    try:
        tree = ast.parse(source, filename=str(path))
    except (SyntaxError, ValueError) as error:
        raise SurveyError(f"could not parse {path}: {error}") from error

    return StaticSurvey(
        path=path,
        lineCount=len(source.splitlines()),
        functions=oversized_functions(tree),
        hits=duplication_hits(tree),
        bareExcepts=bare_except_lines(tree),
    )


def oversized_functions(tree: ast.Module) -> list[FunctionLength]:
    """*every function longer than the 50-line house rule, in source order*

    Length is measured from the `def` line to the last line of the body, so a
    decorator is excluded and the reported line is the one a reader jumps to.

    **Key Arguments:**

    - ``tree`` -- the parsed module

    **Return:**

    - ``functions`` -- the oversized functions, in source order
    """
    functions = []
    for node, qualifiedName in _walk_functions(tree):
        length = (node.end_lineno or node.lineno) - node.lineno + 1
        if length > OVERSIZED_FUNCTION_LINES:
            functions.append(FunctionLength(qualifiedName, node.lineno, length))

    return sorted(functions, key=lambda function: function.line)


def duplication_hits(tree: ast.Module) -> list[Hit]:
    """*every duplication hit in a module, under the definition in this file's docstring*

    **Key Arguments:**

    - ``tree`` -- the parsed module

    **Return:**

    - ``hits`` -- the hits found, in source order
    """
    hits = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            kind = _call_kind(node)
            if kind:
                hits.append(Hit(kind, node.lineno))
        elif isinstance(node, (ast.Assign, ast.AugAssign, ast.AnnAssign)):
            kind = _row_builder_kind(node)
            if kind:
                hits.append(Hit(kind, node.lineno))

    return sorted(hits, key=lambda hit: (hit.line, hit.kind))


def bare_except_lines(tree: ast.Module) -> list[int]:
    """*the line of every bare `except:` in a module*

    **Key Arguments:**

    - ``tree`` -- the parsed module

    **Return:**

    - ``lines`` -- the line numbers, in source order
    """
    return sorted(node.lineno for node in ast.walk(tree) if isinstance(node, ast.ExceptHandler) and node.type is None)


def parse_ruff_counts(ruffJson: str) -> dict[str, int]:
    """*count a ruff JSON report by rule code*

    **Key Arguments:**

    - ``ruffJson`` -- the raw `ruff check --output-format json` report

    **Return:**

    - ``counts`` -- the number of findings per rule code
    """
    try:
        findings = json.loads(ruffJson)
    except json.JSONDecodeError as error:
        raise SurveyError(f"could not read the ruff report: {error}") from error

    # WELL-FORMED JSON OF THE WRONG SHAPE IS A TOOL ERROR, NOT A TRACEBACK: THE SURVEY
    # EXITS 0 OR 2 AND NEVER 1, SO NOTHING HERE MAY ESCAPE AS AN UNCAUGHT EXCEPTION
    if not isinstance(findings, list) or any(not isinstance(finding, dict) for finding in findings):
        raise SurveyError("the ruff report is not a list of findings")

    counts: dict[str, int] = {}
    for finding in findings:
        # A SYNTAX ERROR HAS NO RULE CODE, SO IT IS COUNTED UNDER A PLACEHOLDER
        code = finding.get("code") or "syntax-error"
        if not isinstance(code, str):
            raise SurveyError(f"a ruff finding carries a non-string rule code: {code!r}")
        counts[code] = counts.get(code, 0) + 1

    return counts


def read_coverage_percent(reportPath: Path, modulePath: str) -> float:
    """*one module's covered percentage, out of a coverage JSON report*

    **Key Arguments:**

    - ``reportPath`` -- the path of the JSON report written by `pytest --cov-report=json`
    - ``modulePath`` -- the repository-relative path of the module, as coverage keys it

    **Return:**

    - ``percent`` -- the covered percentage of the surveyed module
    """
    try:
        report = json.loads(Path(reportPath).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise SurveyError(f"could not read the coverage report: {error}") from error

    try:
        return float(report["files"][modulePath]["summary"]["percent_covered"])
    except (KeyError, TypeError, ValueError) as error:
        # AN ABSENT ENTRY MEANS COVERAGE NEVER SAW THE MODULE, WHICH IS NOT 0% COVERAGE
        raise SurveyError(f"the coverage report has no entry for {modulePath}: {error}") from error


def measure_coverage(modulePath: str, repoRoot: Path) -> float:
    """*run the offline suite against one module and report its coverage*

    This is the same measurement CI would produce, since it runs in the repository
    with the repository's own pytest and coverage configuration.

    Coverage measures the whole `soxspipe` package and the module's own figure is read
    back out of the report. Passing the module's path to `--cov` instead measures
    nothing at all, silently: the run passes, the terminal report is empty, and no JSON
    is written.

    **Key Arguments:**

    - ``modulePath`` -- the repository-relative path of the module to measure
    - ``repoRoot`` -- the repository to run pytest in

    **Return:**

    - ``percent`` -- the covered percentage of the module
    """
    # A TEMPORARY DIRECTORY THAT CANNOT BE MADE OR REMOVED IS A TOOL ERROR LIKE ANY
    # OTHER, SO IT LEAVES HERE AS A SurveyError RATHER THAN AS A TRACEBACK
    try:
        with tempfile.TemporaryDirectory() as workspace:
            reportPath = Path(workspace) / "coverage.json"
            command = [sys.executable, "-m", "pytest", "tests/unit", "tests/integration"]
            command += [f"--cov={COVERAGE_PACKAGE}", f"--cov-report=json:{reportPath}", "--cov-fail-under=0", "-q"]
            # A FAILING SUITE STILL MEASURES COVERAGE, AND A BRANCH POINT WITH A RED SUITE
            # IS A FINDING FOR THE TICKET RATHER THAN A REASON TO REFUSE THE SURVEY
            _run(command, repoRoot, allowedStatuses=(0, 1))

            return read_coverage_percent(reportPath, modulePath)
    except OSError as error:
        raise SurveyError(f"could not use a temporary directory for the coverage report: {error}") from error


def count_ruff_findings(modulePath: str, repoRoot: Path) -> dict[str, int]:
    """*the whole-file ruff findings for one module, under the house configuration*

    **Key Arguments:**

    - ``modulePath`` -- the repository-relative path of the module to check
    - ``repoRoot`` -- the repository to run ruff in

    **Return:**

    - ``counts`` -- the number of findings per rule code
    """
    # THE -- STOPS A PATH THAT BEGINS WITH A DASH BEING READ AS AN OPTION
    command = ["ruff", "check", "--output-format", "json", "--force-exclude", "--", modulePath]

    # RUFF EXITS 1 WHEN IT FINDS SOMETHING, WHICH IS THE NORMAL CASE HERE
    return parse_ruff_counts(_run(command, repoRoot, allowedStatuses=(0, 1)))


def measure_docstring_drift(path: Path, repoRoot: Path) -> tuple[int, int]:
    """*the docstring findings and drifted-function count for one module*

    **Key Arguments:**

    - ``path`` -- the path of the module to check
    - ``repoRoot`` -- the repository holding `tools/check_docstrings.py`

    **Return:**

    - ``counts`` -- (findings, driftedFunctions) for the module
    """
    checker = _load_checker(repoRoot / "tools" / "check_docstrings.py")
    findings = checker.check_file(Path(path))

    toolErrors = [finding for finding in findings if finding.kind in checker.TOOL_ERROR_KINDS]
    if toolErrors:
        raise SurveyError(f"the docstring checker could not read {path}: {toolErrors[0].detail}")

    driftKinds = ("undocumented-argument", "phantom-argument")
    drifted = {finding.function for finding in findings if finding.kind in driftKinds}

    return len(findings), len(drifted)


def render_report(survey: Survey) -> str:
    """*render a survey as the Markdown bullets a ticket body carries*

    **Key Arguments:**

    - ``survey`` -- the survey to render

    **Return:**

    - ``report`` -- the Markdown bullet list, without a trailing newline
    """
    huge = [function for function in survey.functions if function.length > UNSPLIT_REASON_LINES]
    savefigLines = [hit.line for hit in survey.hits if hit.kind == "savefig"]

    coverage = "not measured" if survey.coveragePercent is None else f"{survey.coveragePercent:.1f}%"

    lines = [
        f"- {survey.lineCount} lines.",
        f"- Coverage {coverage}.",
        f"- Functions over 50 lines: {_functions_phrase(survey.functions)}",
        f"- Functions over 200 lines: {_functions_phrase(huge)}",
        f"- Duplication hits: {_hits_phrase(survey.hits)}",
        f"- Bare excepts: {_lines_phrase(survey.bareExcepts)}. `plt.savefig` calls: {_lines_phrase(savefigLines)}.",
        f"- Ruff findings: {_ruff_phrase(survey.ruffCounts)}",
        f"- Docstring drift: {_drift_phrase(survey.docstringFindings, survey.driftedFunctions)}",
    ]

    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    """*run the survey from the command line*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments, excluding the program name. Default *None*, i.e. `sys.argv[1:]`

    **Return:**

    - ``status`` -- 0 when the survey was produced, 2 when it could not be
    """
    parser = argparse.ArgumentParser(description="Survey one soxspipe module for its refactor ticket.")
    parser.add_argument("module", type=Path, help="the module to survey, e.g. soxspipe/recipes/soxs_mbias.py")
    parser.add_argument(
        "--no-coverage",
        action="store_true",
        help="skip the coverage measurement, which runs the offline suite and takes about a minute",
    )
    arguments = parser.parse_args(argv)

    try:
        repoRoot = repository_root(Path.cwd())
        modulePath = _repository_relative(arguments.module, repoRoot)
        static = static_survey(repoRoot / modulePath)
        findings, drifted = measure_docstring_drift(repoRoot / modulePath, repoRoot)
        survey = Survey(
            path=Path(modulePath),
            lineCount=static.lineCount,
            functions=static.functions,
            hits=static.hits,
            bareExcepts=static.bareExcepts,
            coveragePercent=None if arguments.no_coverage else measure_coverage(modulePath, repoRoot),
            ruffCounts=count_ruff_findings(modulePath, repoRoot),
            docstringFindings=findings,
            driftedFunctions=drifted,
        )
    except SurveyError as error:
        print(f"module survey: {error}", file=sys.stderr)
        return EXIT_TOOL_ERROR

    print(render_report(survey))

    return EXIT_CLEAN


def repository_root(start: Path) -> Path:
    """*the root of the git repository containing a path*

    **Key Arguments:**

    - ``start`` -- a path inside the repository

    **Return:**

    - ``root`` -- the absolute repository root
    """
    return Path(_run(["git", "rev-parse", "--show-toplevel"], start).strip())


def _repository_relative(path: Path, repoRoot: Path) -> str:
    """*a module path expressed relative to the repository root*

    **Key Arguments:**

    - ``path`` -- the path given on the command line, absolute or relative
    - ``repoRoot`` -- the absolute repository root

    **Return:**

    - ``modulePath`` -- the path relative to the repository root, with forward slashes
    """
    resolved = Path(path).resolve()
    try:
        return resolved.relative_to(repoRoot.resolve()).as_posix()
    except ValueError as error:
        raise SurveyError(f"{path} is outside the repository at {repoRoot}") from error


def _call_kind(node: ast.Call) -> str | None:
    """*the duplication kind of a call, when it is one of the counted calls*

    **Key Arguments:**

    - ``node`` -- the call node

    **Return:**

    - ``kind`` -- the duplication kind, or None when the call is not counted
    """
    name = node.func.attr if isinstance(node.func, ast.Attribute) else getattr(node.func, "id", None)
    if name == "utcnow":
        return "naive-utcnow"
    if name == "savefig":
        return "savefig"

    return None


def _row_builder_kind(node: ast.Assign | ast.AugAssign | ast.AnnAssign) -> str | None:
    """*the duplication kind of an assignment, when it appends one row to a tracked table*

    **Key Arguments:**

    - ``node`` -- the assignment node

    **Return:**

    - ``kind`` -- the duplication kind, or None when the assignment is not a row builder
    """
    targets = node.targets if isinstance(node, ast.Assign) else [node.target]
    attributes = [
        target.attr
        for target in targets
        if isinstance(target, ast.Attribute) and isinstance(target.value, ast.Name) and target.value.id == "self"
    ]
    kinds = [ROW_BUILDER_ATTRIBUTES[attribute] for attribute in attributes if attribute in ROW_BUILDER_ATTRIBUTES]
    if not kinds or not _is_single_row_concat(node.value):
        return None

    return kinds[0]


def _is_single_row_concat(value: ast.expr | None) -> bool:
    """*is this a `pd.concat` that appends exactly one dict literal?*

    Two spellings count, because both appear in this codebase's history: the current
    `pd.DataFrame([{...}])`, and the `pd.Series({...}).to_frame().T` it replaced. A
    concatenation of whole tables has no dict literal to lift into a helper, so it is
    not a row builder and is deliberately not counted.

    **Key Arguments:**

    - ``value`` -- the assigned expression

    **Return:**

    - ``isRowBuilder`` -- True when the concatenation appends a one-row literal
    """
    if not isinstance(value, ast.Call) or _call_name(value) != "concat":
        return False

    return any(isinstance(node, ast.Call) and _is_single_row_literal(node) for node in ast.walk(value))


def _is_single_row_literal(node: ast.Call) -> bool:
    """*does this call build a one-row frame out of a dict literal?*

    **Key Arguments:**

    - ``node`` -- the call node

    **Return:**

    - ``isSingleRow`` -- True for `pd.DataFrame([{...}])` and for `pd.Series({...})`
    """
    name = _call_name(node)
    if not node.args:
        return False

    if name == "DataFrame":
        rows = node.args[0]
        return isinstance(rows, ast.List) and len(rows.elts) == 1 and isinstance(rows.elts[0], ast.Dict)

    return name == "Series" and isinstance(node.args[0], ast.Dict)


def _call_name(node: ast.Call) -> str | None:
    """*the bare name of a called function, whether or not it is attribute access*

    **Key Arguments:**

    - ``node`` -- the call node

    **Return:**

    - ``name`` -- the called name, or None when it cannot be read statically
    """
    if isinstance(node.func, ast.Attribute):
        return node.func.attr

    return getattr(node.func, "id", None)


def _functions_phrase(functions: list[FunctionLength]) -> str:
    """*the report phrase naming a set of oversized functions*

    **Key Arguments:**

    - ``functions`` -- the functions to name

    **Return:**

    - ``phrase`` -- the count and the named functions, or "none."
    """
    if not functions:
        return "none."

    named = ", ".join(f"`{function.name}` ({function.length}, L{function.line})" for function in functions)

    return f"{len(functions)} — {named}."


def _hits_phrase(hits: list[Hit]) -> str:
    """*the report phrase grouping duplication hits by kind*

    **Key Arguments:**

    - ``hits`` -- the hits to group

    **Return:**

    - ``phrase`` -- the total, then each kind with its lines, or "none."
    """
    if not hits:
        return "none."

    grouped = []
    for kind in DUPLICATION_KINDS:
        lines = [hit.line for hit in hits if hit.kind == kind]
        if lines:
            singular, plural = DUPLICATION_LABELS[kind]
            label = singular if len(lines) == 1 else plural
            grouped.append(f"{len(lines)} {label} ({', '.join(f'L{line}' for line in lines)})")

    return f"{len(hits)} — {', '.join(grouped)}."


def _lines_phrase(lines: list[int]) -> str:
    """*the report phrase counting something and listing where it is*

    **Key Arguments:**

    - ``lines`` -- the line numbers

    **Return:**

    - ``phrase`` -- the count and the lines, or "none"
    """
    if not lines:
        return "none"

    return f"{len(lines)} ({', '.join(f'L{line}' for line in lines)})"


def _ruff_phrase(counts: dict[str, int]) -> str:
    """*the report phrase listing the remaining ruff findings by code*

    **Key Arguments:**

    - ``counts`` -- the number of findings per rule code

    **Return:**

    - ``phrase`` -- the total and the per-code breakdown, or "none."
    """
    if not counts:
        return "none."

    ranked = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    named = ", ".join(f"{code} ×{count}" for code, count in ranked)

    return f"{sum(counts.values())} — {named}."


def _drift_phrase(findings: int, drifted: int) -> str:
    """*the report phrase stating the docstring drift*

    **Key Arguments:**

    - ``findings`` -- the total docstring findings
    - ``drifted`` -- the number of functions whose arguments disagree with their docstring

    **Return:**

    - ``phrase`` -- the drifted-function count and the finding total, or "none."
    """
    if not findings:
        return "none."

    return f"{drifted} drifted functions, {findings} findings."


def _walk_functions(tree: ast.Module) -> Iterator[tuple[ast.FunctionDef | ast.AsyncFunctionDef, str]]:
    """*walk every function in a module, carrying its qualified name*

    **Key Arguments:**

    - ``tree`` -- the parsed module

    **Return:**

    - ``functions`` -- a generator of (node, qualifiedName) pairs
    """
    yield from _walk_scope(tree, prefix="")


def _walk_scope(node: ast.AST, prefix: str) -> Iterator[tuple[ast.FunctionDef | ast.AsyncFunctionDef, str]]:
    """*walk one scope, recursing into classes and nested functions*

    **Key Arguments:**

    - ``node`` -- the scope node to walk
    - ``prefix`` -- the dotted name prefix for anything defined in this scope

    **Return:**

    - ``functions`` -- a generator of (node, qualifiedName) pairs
    """
    for child in ast.iter_child_nodes(node):
        if isinstance(child, ast.ClassDef):
            yield from _walk_scope(child, f"{prefix}{child.name}.")
        elif isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
            yield child, f"{prefix}{child.name}"
            yield from _walk_scope(child, f"{prefix}{child.name}.")
        else:
            # A DEFINITION INSIDE AN IF, TRY, FOR, WHILE, WITH OR MATCH BLOCK BELONGS TO THIS SCOPE
            yield from _walk_scope(child, prefix)


def _load_checker(toolPath: Path) -> ModuleType:
    """*import the docstring checker by path, since `tools/` is not an installed package*

    **Key Arguments:**

    - ``toolPath`` -- the path of `tools/check_docstrings.py`

    **Return:**

    - ``checker`` -- the imported module
    """
    spec = importlib.util.spec_from_file_location("check_docstrings", toolPath)
    if spec is None or spec.loader is None:
        raise SurveyError(f"could not import the docstring checker at {toolPath}")

    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
    except (OSError, ImportError, SyntaxError) as error:
        raise SurveyError(f"could not import the docstring checker at {toolPath}: {error}") from error

    return module


def _run(command: list[str], workingDirectory: Path, allowedStatuses: tuple[int, ...] = (0,)) -> str:
    """*run a command and return its standard output*

    **Key Arguments:**

    - ``command`` -- the command and its arguments
    - ``workingDirectory`` -- the directory to run in
    - ``allowedStatuses`` -- the exit statuses that are not failures. Default *(0,)*

    **Return:**

    - ``output`` -- the command's standard output
    """
    try:
        # NO SHELL, A FIXED EXECUTABLE, AND EVERY PATH ARGUMENT PLACED AFTER A -- SEPARATOR
        completed = subprocess.run(command, cwd=workingDirectory, capture_output=True, text=True, check=False)  # noqa: S603
    except OSError as error:
        raise SurveyError(f"could not run {command[0]}: {error}") from error

    if completed.returncode not in allowedStatuses:
        failure = f"{' '.join(command)} failed with status {completed.returncode}"
        raise SurveyError(f"{failure}:\n{completed.stderr.strip()}")

    return completed.stdout


if __name__ == "__main__":
    sys.exit(main())
