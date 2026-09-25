"""Unit tests for the docstring signature checker in `tools/check_docstrings.py`."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

pytestmark = pytest.mark.unit

TOOL_PATH = Path(__file__).resolve().parents[2] / "tools" / "check_docstrings.py"


def _load_checker():
    """*import the checker by path, since `tools/` is not an installed package*"""
    spec = importlib.util.spec_from_file_location("check_docstrings", TOOL_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


checker = _load_checker()


def _check_source(tmp_path: Path, source: str) -> list:
    modulePath = tmp_path / "sample.py"
    modulePath.write_text(source)
    return checker.check_file(modulePath)


def _kinds(findings: list) -> list[str]:
    return [finding.kind for finding in findings]


def _function(source: str):
    return checker.ast.parse(source).body[0]


def test_check_arguments_reports_only_argument_findings(tmp_path: Path) -> None:
    node = _function('''
def measure(frame, save):
    """*measure*"""
    return frame
''')

    findings = checker._check_arguments(node, tmp_path / "sample.py", "measure", "*measure*")

    assert _kinds(findings) == ["no-arguments-section"]
    assert findings[0].detail == "frame, save"


def test_check_return_reports_only_return_findings(tmp_path: Path) -> None:
    node = _function('''
def measure(frame):
    """*measure*"""
    return frame
''')

    findings = checker._check_return(node, tmp_path / "sample.py", "measure", "*measure*")

    assert _kinds(findings) == ["undocumented-return"]
    assert findings[0].detail == ""


def test_agreeing_docstring_produces_no_findings(tmp_path: Path) -> None:
    # ARRANGE
    source = '''
def prepare_frame(frame, save=False):
    """*prepare a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to prepare
    - ``save`` -- save the frame to file. Default *False*

    **Return:**

    - ``frame`` -- the prepared frame
    """
    return frame
'''

    # ACT
    findings = _check_source(tmp_path, source)

    # ASSERT
    assert findings == []


def test_argument_missing_from_docstring_is_a_finding(tmp_path: Path) -> None:
    source = '''
def prepare_frame(frame, save=False):
    """*prepare a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to prepare
    """
    return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["undocumented-argument"]
    assert findings[0].detail == "save"
    assert findings[0].function == "prepare_frame"


def test_documented_argument_not_in_signature_is_a_finding(tmp_path: Path) -> None:
    source = '''
def prepare_frame(frame):
    """*prepare a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to prepare
    - ``save`` -- removed from the signature long ago
    """
    return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["phantom-argument"]
    assert findings[0].detail == "save"


def test_missing_key_arguments_section_is_its_own_finding(tmp_path: Path) -> None:
    source = '''
def prepare_frame(frame):
    """*prepare a frame*"""
    return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["no-arguments-section"]


def test_self_and_cls_are_never_required_to_be_documented(tmp_path: Path) -> None:
    source = '''
class recipe(object):
    """*a recipe*"""

    def run(self):
        """*run the recipe*"""
        return None

    @classmethod
    def build(cls):
        """*build a recipe*"""
        return None
'''

    findings = _check_source(tmp_path, source)

    assert findings == []


def test_init_arguments_are_read_from_the_class_docstring(tmp_path: Path) -> None:
    source = '''
class extractor(object):
    """*extract a spectrum*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    """

    def __init__(self, log, settings):
        return None
'''

    findings = _check_source(tmp_path, source)

    assert findings == []


def test_init_missing_from_the_class_docstring_is_a_finding(tmp_path: Path) -> None:
    source = '''
class extractor(object):
    """*extract a spectrum*

    **Key Arguments:**

    - ``log`` -- logger
    """

    def __init__(self, log, settings):
        return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["undocumented-argument"]
    assert findings[0].detail == "settings"
    assert findings[0].function == "extractor.__init__"


def test_init_documenting_itself_is_checked_when_the_class_has_no_docstring(tmp_path: Path) -> None:
    source = '''
class extractor(object):

    def __init__(self, log, settings):
        """*build an extractor*

        **Key Arguments:**

        - ``log`` -- logger
        """
        return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["undocumented-argument"]
    assert findings[0].detail == "settings"


def test_varargs_and_kwargs_are_checked_by_name(tmp_path: Path) -> None:
    source = '''
def dispatch(*args, **kwargs):
    """*dispatch a call*

    **Key Arguments:**

    - ``args`` -- positional arguments
    """
    return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["undocumented-argument"]
    assert findings[0].detail == "kwargs"


def test_undocumented_return_value_is_a_finding(tmp_path: Path) -> None:
    source = '''
def measure(frame):
    """*measure a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to measure
    """
    return frame.mean()
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["undocumented-return"]


def test_an_empty_return_section_does_not_document_a_return(tmp_path: Path) -> None:
    source = '''
def measure(frame):
    """*measure a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to measure

    **Return:**
    """
    return frame.mean()
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["undocumented-return"]


def test_a_return_section_documented_without_the_canonical_bullet_is_accepted(tmp_path: Path) -> None:
    # SEVERAL MODULES DOCUMENT RETURNS WITHOUT THE ``--`` SEPARATOR, AND THAT STILL COUNTS AS DOCUMENTED
    source = '''
def measure(frame):
    """*measure a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to measure

    **Return:**

    - ``mean``
    """
    return frame.mean()
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == []


def test_a_yield_counts_as_a_returned_value(tmp_path: Path) -> None:
    source = '''
def rows(frame):
    """*walk the rows*

    **Key Arguments:**

    - ``frame`` -- the frame to walk
    """
    yield frame
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["undocumented-return"]


def test_documented_return_with_no_returned_value_is_a_finding(tmp_path: Path) -> None:
    source = '''
def measure(frame):
    """*measure a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to measure

    **Return:**

    - ``value`` -- the measurement
    """
    self.value = frame
    return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["phantom-return"]


def test_a_return_inside_a_nested_function_does_not_count(tmp_path: Path) -> None:
    source = '''
def measure(frame):
    """*measure a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to measure
    """

    def inner():
        """*inner helper*

        **Return:**

        - ``answer`` -- the answer
        """
        return 42

    inner()
    return None
'''

    findings = _check_source(tmp_path, source)

    assert findings == []


def test_section_header_variants_without_a_colon_are_recognised(tmp_path: Path) -> None:
    source = '''
def measure(frame):
    """*measure a frame*

    **Key Arguments**

    - ``frame`` -- the frame to measure

    **Returns:**

    - ``value`` -- the measurement
    """
    return frame.mean()
'''

    findings = _check_source(tmp_path, source)

    assert findings == []


def test_a_function_defined_inside_a_control_flow_block_is_checked(tmp_path: Path) -> None:
    source = '''
import sys

if sys.version_info >= (3, 12):

    def prepare_frame(frame, save):
        """*prepare a frame*

        **Key Arguments:**

        - ``frame`` -- the frame to prepare
        """
        return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["undocumented-argument"]
    assert findings[0].function == "prepare_frame"


def test_a_method_of_a_class_defined_inside_a_try_block_is_checked(tmp_path: Path) -> None:
    source = '''
try:

    class recipe(object):
        """*a recipe*"""

        def run(self, frame):
            """*run the recipe*"""
            return None

except ImportError:
    pass
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["no-arguments-section"]
    assert findings[0].function == "recipe.run"


def test_an_argument_bullet_missing_its_separator_is_a_finding(tmp_path: Path) -> None:
    source = '''
def prepare_frame(frame):
    """*prepare a frame*

    **Key Arguments:**

    - ``frame`` the frame to prepare
    """
    return None
'''

    findings = _check_source(tmp_path, source)

    assert _kinds(findings) == ["malformed-argument-bullet"]
    assert findings[0].detail == "frame"


def test_a_malformed_bullet_still_counts_as_documenting_its_argument(tmp_path: Path) -> None:
    source = '''
def prepare_frame(frame, save):
    """*prepare a frame*

    **Key Arguments:**

    - ``frame`` the frame to prepare
    - ``save`` -- save the frame. Default *False*
    """
    return None
'''

    findings = _check_source(tmp_path, source)

    # THE MALFORMED BULLET IS REPORTED AS PUNCTUATION, NOT AS A MISSING ARGUMENT
    assert _kinds(findings) == ["malformed-argument-bullet"]


def test_an_unreadable_file_is_reported_rather_than_raised(tmp_path: Path) -> None:
    findings = checker.check_file(tmp_path / "does_not_exist.py")

    assert _kinds(findings) == ["unreadable"]


def test_main_exits_two_when_a_file_cannot_be_read(tmp_path: Path) -> None:
    modulePath = tmp_path / "sample.py"
    modulePath.write_bytes(b'def measure():\n    """*\xff\xfe*"""\n    return None\n')

    assert checker.main([str(modulePath)]) == 2


def test_a_function_with_no_docstring_is_not_a_drift_finding(tmp_path: Path) -> None:
    source = """
def measure(frame):
    return frame.mean()
"""

    findings = _check_source(tmp_path, source)

    assert findings == []


def test_findings_carry_the_path_and_line_number_of_the_function(tmp_path: Path) -> None:
    source = '''
def measure(frame, save):
    """*measure a frame*

    **Key Arguments:**

    - ``frame`` -- the frame to measure
    """
    return None
'''

    findings = _check_source(tmp_path, source)

    assert findings[0].line == 2
    assert findings[0].path.name == "sample.py"


def test_a_syntax_error_is_reported_rather_than_raised(tmp_path: Path) -> None:
    findings = _check_source(tmp_path, "def broken(:\n    pass\n")

    assert _kinds(findings) == ["unparseable"]


def test_find_python_files_walks_a_directory_and_accepts_a_file(tmp_path: Path) -> None:
    (tmp_path / "package").mkdir()
    (tmp_path / "package" / "one.py").write_text("")
    (tmp_path / "package" / "two.py").write_text("")
    (tmp_path / "package" / "notes.txt").write_text("")

    fromDirectory = checker.find_python_files([tmp_path / "package"])
    fromFile = checker.find_python_files([tmp_path / "package" / "one.py"])

    assert [path.name for path in fromDirectory] == ["one.py", "two.py"]
    assert [path.name for path in fromFile] == ["one.py"]


def test_main_exits_non_zero_when_a_module_has_findings(tmp_path: Path) -> None:
    modulePath = tmp_path / "sample.py"
    modulePath.write_text('''
def measure(frame, save):
    """*measure*

    **Key Arguments:**

    - ``frame`` -- a frame
    """
    return None
''')

    assert checker.main([str(modulePath)]) == 1


def test_main_exits_zero_for_a_clean_module(tmp_path: Path) -> None:
    modulePath = tmp_path / "sample.py"
    modulePath.write_text('''
def measure(frame):
    """*measure*

    **Key Arguments:**

    - ``frame`` -- a frame
    """
    return None
''')

    assert checker.main([str(modulePath)]) == 0
