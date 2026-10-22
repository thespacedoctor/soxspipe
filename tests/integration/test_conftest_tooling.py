"""*tests for the suite's own conftest hooks, run against a throwaway project*

Each test copies the real `tests/conftest.py` into a temporary project and runs pytest on it
in a subprocess, so the hooks run exactly as they do in the suite but can be pointed at a
module the real suite would never contain.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

REPO_ROOT = Path(__file__).resolve().parents[2]
CONFTEST_PATH = REPO_ROOT / "tests" / "conftest.py"

PYPROJECT = """[tool.pytest.ini_options]
addopts = ["--strict-markers"]
markers = ["unit", "integration", "real_data"]
"""

UNMARKED_MODULE = "def test_without_a_category():\n    assert True\n"

MARKED_MODULE = "import pytest\n\npytestmark = pytest.mark.unit\n\ndef test_with_a_category():\n    assert True\n"

RECORD_TMP_PATH_MODULE = """import os
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit


def test_records_its_tmp_path(tmp_path):
    with Path(os.environ["RECORD_FILE"]).open("a", encoding="utf-8") as handle:
        handle.write(f"{tmp_path}\\n")
"""

SCRUBBED_PREFIXES = ("PYTEST_ADDOPTS", "PYTEST_PLUGINS", "PYTEST_DISABLE_PLUGIN_AUTOLOAD", "PYTEST_XDIST_")

MISSING_CATEGORY_MESSAGE = "must declare exactly one execution category"


def _make_project(root: Path, **modules: str) -> Path:
    """*build a throwaway pytest project holding the real conftest and the given test modules*

    **Key Arguments:**

    - ``root`` -- the directory to create the project in
    - ``modules`` -- test module source, keyed by file name stem

    **Return:**

    - ``project`` -- the project directory
    """
    project = root / "project"
    project.mkdir()
    (project / "pyproject.toml").write_text(PYPROJECT, encoding="utf-8")
    shutil.copy(CONFTEST_PATH, project / "conftest.py")
    for stem, source in modules.items():
        (project / f"{stem}.py").write_text(source, encoding="utf-8")

    return project


def _run_pytest(
    project: Path, *arguments: str, tempRoot: Path, record: Path | None = None
) -> subprocess.CompletedProcess:
    """*run pytest on a throwaway project in a subprocess*

    **Key Arguments:**

    - ``project`` -- the project to run in
    - ``arguments`` -- extra pytest arguments
    - ``tempRoot`` -- the directory the subprocess treats as its system temporary directory
    - ``record`` -- a file the recording test module appends each `tmp_path` to. Default *None*

    **Return:**

    - ``completed`` -- the finished process, with output captured as text
    """
    # THE OUTER RUN MAY BE AN XDIST WORKER OR CARRY PYTEST OPTIONS, NONE OF WHICH THE INNER RUN MAY INHERIT
    inherited = {name: value for name, value in os.environ.items() if not name.startswith(SCRUBBED_PREFIXES)}
    # A WIDE TERMINAL STOPS PYTEST TRUNCATING THE FAILED LINE'S MESSAGE
    environment = {**inherited, "TMPDIR": str(tempRoot), "COLUMNS": "300"}
    if record is not None:
        environment["RECORD_FILE"] = str(record)

    command = [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", "-rfE", *arguments]

    return subprocess.run(command, cwd=project, env=environment, capture_output=True, text=True, check=False)  # noqa: S603


@pytest.fixture
def temp_root(tmp_path) -> Path:
    """*a directory for the subprocess to use as its system temporary directory*

    **Return:**

    - ``tempRoot`` -- the empty directory
    """
    tempRoot = tmp_path / "systmp"
    tempRoot.mkdir()

    return tempRoot


def test_an_unmarked_test_aborts_a_serial_run_with_a_usage_error(tmp_path, temp_root):
    # ARRANGE
    project = _make_project(tmp_path, test_plain=UNMARKED_MODULE)

    # ACT
    completed = _run_pytest(project, tempRoot=temp_root)

    # ASSERT
    assert completed.returncode == pytest.ExitCode.USAGE_ERROR
    assert MISSING_CATEGORY_MESSAGE in completed.stderr + completed.stdout


def test_an_unmarked_test_errors_alone_under_xdist_without_an_internal_error(tmp_path, temp_root):
    # ARRANGE
    project = _make_project(tmp_path, test_plain=UNMARKED_MODULE, test_marked=MARKED_MODULE)

    # ACT
    completed = _run_pytest(project, "-n", "2", tempRoot=temp_root)

    # ASSERT
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert "INTERNALERROR" not in output
    errorLines = [line for line in output.splitlines() if line.startswith("ERROR")]
    assert len(errorLines) == 1
    assert "test_plain.py::test_without_a_category" in errorLines[0]
    assert MISSING_CATEGORY_MESSAGE in errorLines[0]
    assert "1 passed, 1 error" in output


@pytest.mark.parametrize(
    "decorator",
    [
        "@pytest.mark.skip(reason='r')",
        "@pytest.mark.skipif(True, reason='r')",
        "@pytest.mark.xfail(reason='r')",
        "@pytest.mark.xfail(run=False, reason='r')",
        "@pytest.mark.xfail(strict=False, raises=ValueError, reason='r')",
    ],
)
def test_an_unmarked_test_cannot_hide_behind_a_skip_or_xfail_mark_under_xdist(tmp_path, temp_root, decorator):
    # ARRANGE
    source = f"import pytest\n\n{decorator}\ndef test_hidden():\n    assert True\n"
    project = _make_project(tmp_path, test_hidden=source)

    # ACT
    completed = _run_pytest(project, "-n", "2", tempRoot=temp_root)

    # ASSERT
    output = completed.stdout + completed.stderr
    assert completed.returncode != 0
    assert MISSING_CATEGORY_MESSAGE in output
    assert "INTERNALERROR" not in output


def test_an_unmarked_test_cannot_hide_behind_a_skipping_fixture_under_xdist(tmp_path, temp_root):
    # ARRANGE
    source = (
        "import pytest\n\n@pytest.fixture\ndef skipper():\n    pytest.skip('r')\n\n"
        "def test_hidden(skipper):\n    assert True\n"
    )
    project = _make_project(tmp_path, test_hidden=source)

    # ACT
    completed = _run_pytest(project, "-n", "2", tempRoot=temp_root)

    # ASSERT
    assert completed.returncode != 0
    assert MISSING_CATEGORY_MESSAGE in completed.stdout + completed.stderr


def test_an_explicit_basetemp_wins_over_the_checkout_default(tmp_path, temp_root):
    # ARRANGE
    project = _make_project(tmp_path, test_record=RECORD_TMP_PATH_MODULE)
    record = tmp_path / "record.txt"
    explicit = tmp_path / "explicit-base"

    # ACT
    completed = _run_pytest(project, f"--basetemp={explicit}", tempRoot=temp_root, record=record)

    # ASSERT
    assert completed.returncode == 0, completed.stdout
    recorded = Path(record.read_text(encoding="utf-8").strip())
    assert recorded.is_relative_to(explicit)
    assert list(temp_root.iterdir()) == []


def test_tmp_path_lands_in_a_directory_unique_to_the_checkout(tmp_path, temp_root):
    # ARRANGE
    project = _make_project(tmp_path, test_record=RECORD_TMP_PATH_MODULE)
    record = tmp_path / "record.txt"

    # ACT
    completed = _run_pytest(project, tempRoot=temp_root, record=record)

    # ASSERT
    assert completed.returncode == 0, completed.stdout
    recorded = Path(record.read_text(encoding="utf-8").strip())
    baseTemps = [path for path in temp_root.iterdir() if path.name.startswith("soxspipe-pytest-")]
    assert len(baseTemps) == 1
    # THE LAST NAME PART IS THE PROCESS ID OF THE RUN, SO TWO RUNS IN ONE CHECKOUT DO NOT SHARE A DIRECTORY
    assert baseTemps[0].name.rsplit("-", 1)[1].isdecimal()
    assert recorded.is_relative_to(baseTemps[0])
    assert not recorded.is_relative_to(project)


def test_xdist_workers_get_a_subfolder_of_the_checkout_basetemp(tmp_path, temp_root):
    # ARRANGE
    project = _make_project(tmp_path, test_record=RECORD_TMP_PATH_MODULE, test_marked=MARKED_MODULE)
    record = tmp_path / "record.txt"

    # ACT
    completed = _run_pytest(project, "-n", "2", tempRoot=temp_root, record=record)

    # ASSERT
    assert completed.returncode == 0, completed.stdout
    recorded = Path(record.read_text(encoding="utf-8").strip())
    baseTemp = next(path for path in temp_root.iterdir() if path.name.startswith("soxspipe-pytest-"))
    assert recorded.relative_to(baseTemp).parts[0].startswith("popen-gw")
