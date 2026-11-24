"""Runs the required-tests workflow's diff step in a scratch git repository.

The step decides whether the required suite runs: only a change whose every file is the root `CHANGES.md` skips it.
"""

import os
import subprocess
from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.integration

WORKFLOW_PATH = Path(__file__).parents[2] / ".github" / "workflows" / "tests.yml"


def load_diff_script() -> str:
    """Return the shell script of the step that decides whether the suite runs."""
    # BASELOADER BUILDS ONLY str/list/dict; safe_load WOULD TURN THE `on` KEY INTO True
    workflow = yaml.load(WORKFLOW_PATH.read_text(), Loader=yaml.BaseLoader)  # noqa: S506
    steps = workflow["jobs"]["changes"]["steps"]
    script = next((step["run"] for step in steps if step.get("id") == "diff"), None)
    assert script is not None, "the changes job has no step with id `diff`"
    return script


def isolated_env(homePath: Path) -> dict[str, str]:
    """Return an environment that hides the caller's git config and any inherited GIT_* variables."""
    env = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    return {**env, "HOME": str(homePath), "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_NOSYSTEM": "1"}


class ScratchRepo:
    """A git repository in a temporary folder, driven through an isolated environment."""

    def __init__(self, rootPath: Path) -> None:
        self.path = rootPath / "repo"
        self.path.mkdir()
        self.env = isolated_env(rootPath)
        self.git("init", "--quiet")
        self.git("config", "user.email", "ci@example.com")
        self.git("config", "user.name", "ci")
        self.commit_files(["CHANGES.md", "soxspipe/module.py"])

    def git(self, *args: str) -> str:
        """Run git in the repository and return its stripped standard output."""
        result = subprocess.run(  # noqa: S603
            ["git", *args],  # noqa: S607
            cwd=self.path,
            env=self.env,
            check=True,
            capture_output=True,
            text=True,
        )
        return result.stdout.strip()

    def commit_files(self, names: list[str]) -> str:
        """Append one line to each named file, commit them, and return the new commit SHA."""
        for name in names:
            filePath = self.path / name
            filePath.parent.mkdir(parents=True, exist_ok=True)
            with filePath.open("a") as handle:
                handle.write("line\n")
        self.git("add", "--all")
        self.git("commit", "--quiet", "--message", "change")
        return self.git("rev-parse", "HEAD")

    def run_diff_step(self, outputPath: Path, eventEnv: dict[str, str]) -> str:
        """Run the diff step's script in the repository and return the `code` value it writes."""
        outputPath.write_text("")
        env = {**self.env, "GITHUB_OUTPUT": str(outputPath), **eventEnv}
        subprocess.run(["bash", "-c", load_diff_script()], cwd=self.path, env=env, check=True, capture_output=True)  # noqa: S603, S607
        lines = outputPath.read_text().splitlines()
        assert len(lines) == 1, f"expected one output line, got {lines}"
        return lines[0].removeprefix("code=")


@pytest.fixture
def repo(tmp_path: Path) -> ScratchRepo:
    """Return a scratch repository with one commit that holds a changelog and a module."""
    return ScratchRepo(tmp_path)


@pytest.mark.parametrize(
    ("changedFiles", "expected"),
    [
        (["CHANGES.md"], "false"),
        (["CHANGES.md", "soxspipe/module.py"], "true"),
        (["soxspipe/module.py"], "true"),
        (["docs/CHANGES.md"], "true"),
    ],
)
def test_pull_request_skips_the_suite_only_for_a_changelog_only_change(
    repo: ScratchRepo, tmp_path: Path, changedFiles: list[str], expected: str
) -> None:
    """Only a pull request whose every changed file is the root CHANGES.md skips the suite."""
    baseSha = repo.git("rev-parse", "HEAD")
    headSha = repo.commit_files(changedFiles)

    code = repo.run_diff_step(
        tmp_path / "output", {"EVENT_NAME": "pull_request", "PR_BASE_SHA": baseSha, "PR_HEAD_SHA": headSha}
    )

    assert code == expected


def test_pull_request_ignores_changes_made_on_the_base_branch_after_it_forked(
    repo: ScratchRepo, tmp_path: Path
) -> None:
    """The three-dot range compares the head with the merge base, so later base-branch commits do not count."""
    forkSha = repo.git("rev-parse", "HEAD")
    baseBranch = repo.git("rev-parse", "--abbrev-ref", "HEAD")
    repo.git("switch", "--quiet", "--create", "pr-branch")
    headSha = repo.commit_files(["CHANGES.md"])
    repo.git("switch", "--quiet", baseBranch)
    baseSha = repo.commit_files(["soxspipe/module.py"])
    assert repo.git("merge-base", baseSha, headSha) == forkSha

    code = repo.run_diff_step(
        tmp_path / "output", {"EVENT_NAME": "pull_request", "PR_BASE_SHA": baseSha, "PR_HEAD_SHA": headSha}
    )

    assert code == "false"


def test_push_of_a_changelog_only_commit_skips_the_suite(repo: ScratchRepo, tmp_path: Path) -> None:
    """A push is judged by the range between the old and new branch tips."""
    beforeSha = repo.git("rev-parse", "HEAD")
    pushSha = repo.commit_files(["CHANGES.md"])

    code = repo.run_diff_step(
        tmp_path / "output", {"EVENT_NAME": "push", "PUSH_BEFORE_SHA": beforeSha, "PUSH_SHA": pushSha}
    )

    assert code == "false"


@pytest.mark.parametrize("beforeSha", ["0" * 40, "", "f" * 40])
def test_push_without_a_usable_base_runs_the_suite(repo: ScratchRepo, tmp_path: Path, beforeSha: str) -> None:
    """A new branch (all-zero base), a missing base or an unknown commit runs the suite rather than guess."""
    pushSha = repo.commit_files(["CHANGES.md"])

    code = repo.run_diff_step(
        tmp_path / "output", {"EVENT_NAME": "push", "PUSH_BEFORE_SHA": beforeSha, "PUSH_SHA": pushSha}
    )

    assert code == "true"
