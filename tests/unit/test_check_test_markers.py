"""*tests for the static execution-category marker checker*

The checker is a pure AST pass, so these tests feed it source text directly and only touch
the filesystem for the CLI path.
"""

from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "tools" / "check_test_markers.py"
CONFTEST_PATH = REPO_ROOT / "tests" / "conftest.py"


def _load_check_test_markers() -> ModuleType:
    """*import the checker from tools/, which is not an installed package*

    **Return:**

    - ``module`` -- the imported `check_test_markers` module
    """
    specification = importlib.util.spec_from_file_location("check_test_markers", MODULE_PATH)
    if specification is None or specification.loader is None:
        raise AssertionError(f"the marker checker is not importable from {MODULE_PATH}")

    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)

    return module


check_test_markers = _load_check_test_markers()


def _names(source: str) -> list[str]:
    """*the names of the tests in a source text that lack exactly one category*

    **Key Arguments:**

    - ``source`` -- the Python source of a test module

    **Return:**

    - ``names`` -- the offending test names, in source order
    """
    return [finding.name for finding in check_test_markers.find_unmarked_tests(source, "test_x.py")]


SPELLING_MESSAGE = "use pytest.mark.<category> (house spelling)"
TUPLE_MESSAGE = "pytestmark must be one mark or a list of marks, not a tuple"


def _details(source: str) -> list[tuple[int, str]]:
    """*the line and message of every finding that is not a missing category*

    **Key Arguments:**

    - ``source`` -- the Python source of a test module

    **Return:**

    - ``details`` -- (line, message) pairs, in source order
    """
    findings = check_test_markers.find_unmarked_tests(source, "test_x.py")

    return [(finding.line, finding.detail) for finding in findings if finding.detail]


def test_accepts_a_module_level_single_mark():
    source = "import pytest\npytestmark = pytest.mark.unit\n\ndef test_a():\n    pass\n"

    assert _names(source) == []


def test_accepts_a_module_level_list_with_one_category_and_other_marks():
    source = "import pytest\npytestmark = [pytest.mark.integration, pytest.mark.slow]\n\ndef test_a():\n    pass\n"

    assert _names(source) == []


def test_rejects_a_module_level_tuple_because_pytest_cannot_collect_it():
    source = "import pytest\npytestmark = (pytest.mark.real_data,)\n\ndef test_a():\n    pass\n"

    assert _details(source) == [(2, TUPLE_MESSAGE)]


def test_rejects_a_class_level_tuple():
    source = (
        "import pytest\n\nclass TestThing:\n    pytestmark = (pytest.mark.unit,)\n\n"
        "    def test_a(self):\n        pass\n"
    )

    assert _details(source) == [(4, TUPLE_MESSAGE)]


def test_accepts_a_function_decorator():
    source = "import pytest\n\n@pytest.mark.unit\ndef test_a():\n    pass\n"

    assert _names(source) == []


def test_accepts_a_called_mark_decorator():
    source = "import pytest\n\n@pytest.mark.unit()\ndef test_a():\n    pass\n"

    assert _names(source) == []


def test_accepts_a_category_with_a_parametrize_decorator():
    source = "import pytest\n\n@pytest.mark.unit\n@pytest.mark.parametrize('x', [1, 2])\ndef test_a(x):\n    pass\n"

    assert _names(source) == []


def test_accepts_an_async_test_function():
    source = "import pytest\npytestmark = pytest.mark.unit\n\nasync def test_a():\n    pass\n"

    assert _names(source) == []


def test_ignores_a_decorator_that_is_not_an_attribute_access():
    source = "import functools\n\n@functools\ndef test_a():\n    pass\n"

    assert _names(source) == ["test_a"]


def test_reports_a_function_with_no_category():
    source = "def test_a():\n    pass\n"

    assert _names(source) == ["test_a"]


def test_reports_a_function_with_only_a_non_category_mark():
    source = "import pytest\npytestmark = pytest.mark.slow\n\ndef test_a():\n    pass\n"

    assert _names(source) == ["test_a"]


def test_reports_a_function_with_two_categories():
    source = "import pytest\npytestmark = pytest.mark.unit\n\n@pytest.mark.integration\ndef test_a():\n    pass\n"

    assert _names(source) == ["test_a"]


def test_accepts_the_same_category_declared_twice():
    source = "import pytest\npytestmark = pytest.mark.unit\n\n@pytest.mark.unit\ndef test_a():\n    pass\n"

    assert _names(source) == []


def test_ignores_functions_that_are_not_tests():
    source = (
        "def helper():\n    pass\n\ndef _test_private():\n    pass\n\n"
        "class Helper:\n    def test_b(self):\n        pass\n"
    )

    assert _names(source) == []


def test_accepts_a_class_decorator_for_its_methods():
    source = "import pytest\n\n@pytest.mark.unit\nclass TestThing:\n    def test_a(self):\n        pass\n"

    assert _names(source) == []


def test_accepts_a_class_level_pytestmark():
    source = (
        "import pytest\n\nclass TestThing:\n    pytestmark = [pytest.mark.integration]\n\n"
        "    def test_a(self):\n        pass\n"
    )

    assert _names(source) == []


def test_accepts_a_category_on_the_method():
    source = "import pytest\n\nclass TestThing:\n    @pytest.mark.unit\n    def test_a(self):\n        pass\n"

    assert _names(source) == []


def test_reports_a_method_with_no_category_by_its_qualified_name():
    source = "class TestThing:\n    def test_a(self):\n        pass\n"

    assert _names(source) == ["TestThing.test_a"]


def test_reports_a_method_that_adds_a_second_category_to_its_class():
    source = (
        "import pytest\n\n@pytest.mark.unit\nclass TestThing:\n"
        "    @pytest.mark.integration\n    def test_a(self):\n        pass\n"
    )

    assert _names(source) == ["TestThing.test_a"]


def test_module_category_covers_a_class_method():
    source = "import pytest\npytestmark = pytest.mark.unit\n\nclass TestThing:\n    def test_a(self):\n        pass\n"

    assert _names(source) == []


def test_a_module_mark_added_with_an_augmented_assignment_counts():
    source = (
        "import pytest\npytestmark = [pytest.mark.slow]\npytestmark += [pytest.mark.unit]\n\ndef test_a():\n    pass\n"
    )

    assert _names(source) == []


def test_reports_the_line_of_the_offending_definition():
    source = "import pytest\npytestmark = pytest.mark.unit\n\n@pytest.mark.integration\ndef test_a():\n    pass\n"

    findings = check_test_markers.find_unmarked_tests(source, "tests/unit/test_x.py")

    assert [finding.line for finding in findings] == [5]


def test_formats_a_finding_as_path_line_and_message():
    finding = check_test_markers.Finding("tests/unit/test_x.py", 7, "test_a")

    assert str(finding) == "tests/unit/test_x.py:7: test_a must declare exactly one execution category"


def test_a_later_assignment_replaces_an_earlier_one():
    source = (
        "import pytest\npytestmark = pytest.mark.unit\npytestmark = pytest.mark.integration\n\n"
        "def test_a():\n    pass\n"
    )

    assert _names(source) == []


def test_a_later_assignment_can_drop_the_category():
    source = (
        "import pytest\npytestmark = [pytest.mark.unit]\npytestmark = [pytest.mark.slow]\n\ndef test_a():\n    pass\n"
    )

    assert _names(source) == ["test_a"]


def test_an_annotated_assignment_replaces_the_value():
    source = (
        "import pytest\npytestmark = [pytest.mark.unit]\n"
        "pytestmark: list = [pytest.mark.integration]\n\ndef test_a():\n    pass\n"
    )

    assert _names(source) == []


def test_an_augmented_assignment_extends_the_value_and_can_add_a_second_category():
    source = (
        "import pytest\npytestmark = [pytest.mark.unit]\n"
        "pytestmark += [pytest.mark.integration]\n\ndef test_a():\n    pass\n"
    )

    assert _names(source) == ["test_a"]


def test_a_class_inherits_the_category_of_a_base_class_in_the_same_module():
    source = (
        "import pytest\n\n@pytest.mark.unit\nclass TestBase:\n    pass\n\n"
        "class TestChild(TestBase):\n    def test_a(self):\n        pass\n"
    )

    assert _names(source) == []


def test_a_class_inherits_the_class_level_pytestmark_of_a_base_through_a_chain():
    source = (
        "import pytest\n\nclass Root:\n    pytestmark = [pytest.mark.integration]\n\n"
        "class Middle(Root):\n    pass\n\nclass TestChild(Middle):\n    def test_a(self):\n        pass\n"
    )

    assert _names(source) == []


def test_a_base_category_plus_a_child_category_is_two_categories():
    source = (
        "import pytest\n\n@pytest.mark.unit\nclass TestBase:\n    pass\n\n"
        "@pytest.mark.integration\nclass TestChild(TestBase):\n    def test_a(self):\n        pass\n"
    )

    assert _names(source) == ["TestChild.test_a"]


def test_a_base_class_from_another_module_is_not_followed():
    source = "from elsewhere import TestBase\n\nclass TestChild(TestBase):\n    def test_a(self):\n        pass\n"

    assert _names(source) == ["TestChild.test_a"]


def test_a_cycle_of_base_classes_terminates_and_still_merges_the_category():
    cyclic = (
        "import pytest\n\n@pytest.mark.unit\nclass TestA(TestB):\n    pass\n\n"
        "class TestB(TestA):\n    def test_a(self):\n        pass\n"
    )
    withoutMerging = cyclic.replace("class TestB(TestA):", "class TestB:")

    assert _names(cyclic) == []
    assert _names(withoutMerging) == ["TestB.test_a"]


def test_a_nested_class_inherits_the_category_of_a_base_in_the_same_enclosing_class():
    source = (
        "import pytest\n\nclass TestOuter:\n    @pytest.mark.unit\n    class TestBase:\n        pass\n\n"
        "    class TestChild(TestBase):\n        def test_a(self):\n            pass\n"
    )

    assert _names(source) == []


def test_a_nested_class_inherits_the_category_of_a_module_level_base():
    source = (
        "import pytest\n\n@pytest.mark.unit\nclass TestBase:\n    pass\n\n"
        "class TestOuter:\n    class TestChild(TestBase):\n        def test_a(self):\n            pass\n"
    )

    assert _names(source) == []


def test_a_nested_class_does_not_see_a_base_defined_in_a_sibling_class():
    source = (
        "import pytest\n\nclass TestOne:\n    @pytest.mark.unit\n    class TestBase:\n        pass\n\n"
        "class TestTwo:\n    class TestChild(TestBase):\n        def test_a(self):\n            pass\n"
    )

    assert _names(source) == ["TestTwo.TestChild.test_a"]


def test_an_annotation_without_a_value_leaves_pytestmark_alone():
    source = "import pytest\npytestmark = pytest.mark.unit\npytestmark: list\n\ndef test_a():\n    pass\n"

    assert _names(source) == []


def test_an_augmented_assignment_to_an_attribute_is_not_the_module_pytestmark():
    source = (
        "import pytest\nimport other\npytestmark = pytest.mark.unit\n"
        "other.pytestmark += [pytest.mark.integration]\n\ndef test_a():\n    pass\n"
    )

    assert _names(source) == []


def test_the_house_spelling_gives_no_spelling_finding():
    source = (
        "import pytest\npytestmark = [pytest.mark.unit]\n\n"
        "@pytest.mark.parametrize('x', [1])\ndef test_a(x):\n    pass\n"
    )

    assert _details(source) == []


def test_rejects_a_mark_imported_from_pytest():
    source = "from pytest import mark\n\n@mark.unit\ndef test_a():\n    pass\n"

    assert _details(source) == [(1, SPELLING_MESSAGE), (3, SPELLING_MESSAGE)]


def test_rejects_a_mark_imported_from_pytest_under_another_name():
    source = "from pytest import mark as m\n\n@m.unit\ndef test_a():\n    pass\n"

    assert _details(source) == [(1, SPELLING_MESSAGE), (3, SPELLING_MESSAGE)]


def test_rejects_pytest_imported_under_an_alias():
    source = "import pytest as p\n\n@p.mark.unit\ndef test_a():\n    pass\n"

    assert _details(source) == [(1, SPELLING_MESSAGE), (3, SPELLING_MESSAGE)]


def test_rejects_a_category_decorator_on_an_unknown_owner():
    source = "import pytest\nimport other\n\n@other.mark.unit\ndef test_a():\n    pass\n"

    assert _details(source) == [(4, SPELLING_MESSAGE)]


def test_rejects_a_pytestmark_value_that_is_not_a_pytest_mark():
    source = "import pytest\nUNIT = pytest.mark.unit\npytestmark = UNIT\n\ndef test_a():\n    pass\n"

    assert _details(source) == [(3, SPELLING_MESSAGE)]


def test_rejects_a_pytestmark_list_element_that_is_not_a_pytest_mark():
    source = "import pytest\nUNIT = pytest.mark.unit\npytestmark = [UNIT]\n\ndef test_a():\n    pass\n"

    assert _details(source) == [(3, SPELLING_MESSAGE)]


def test_pytest_param_marks_are_ignored():
    source = (
        "import pytest\npytestmark = pytest.mark.unit\n\n"
        "@pytest.mark.parametrize('x', [pytest.param(1, marks=pytest.mark.integration)])\n"
        "def test_a(x):\n    pass\n"
    )

    assert _names(source) == [] and _details(source) == []


def test_a_class_marked_not_a_test_is_still_checked():
    source = "class TestThing:\n    __test__ = False\n\n    def test_a(self):\n        pass\n"

    assert _names(source) == ["TestThing.test_a"]


def test_formats_a_spelling_finding_as_path_line_and_its_problem():
    finding = check_test_markers.Finding("tests/unit/test_x.py", 3, "", SPELLING_MESSAGE)

    assert str(finding) == f"tests/unit/test_x.py:3: {SPELLING_MESSAGE}"


def test_the_category_set_matches_the_one_conftest_enforces_at_runtime():
    conftestTree = ast.parse(CONFTEST_PATH.read_text(encoding="utf-8"))
    assignments = [
        node.value
        for node in ast.walk(conftestTree)
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == "_CATEGORY_MARKERS" for t in node.targets)
    ]

    runtimeCategories = ast.literal_eval(assignments[0].args[0])

    assert set(runtimeCategories) == set(check_test_markers.CATEGORY_MARKERS)


def test_main_returns_zero_for_a_clean_file(tmp_path, capsys):
    path = tmp_path / "test_ok.py"
    path.write_text("import pytest\npytestmark = pytest.mark.unit\n\ndef test_a():\n    pass\n", encoding="utf-8")

    status = check_test_markers.main([str(path)])

    assert status == 0
    assert capsys.readouterr().out == ""


def test_main_returns_one_and_prints_each_finding(tmp_path, capsys):
    path = tmp_path / "test_bad.py"
    path.write_text("def test_a():\n    pass\n", encoding="utf-8")

    status = check_test_markers.main([str(path)])

    assert status == 1
    assert f"{path}:1: test_a must declare exactly one execution category" in capsys.readouterr().out


def test_main_checks_every_path_it_is_given(tmp_path, capsys):
    first = tmp_path / "test_one.py"
    second = tmp_path / "test_two.py"
    first.write_text("def test_a():\n    pass\n", encoding="utf-8")
    second.write_text("def test_b():\n    pass\n", encoding="utf-8")

    status = check_test_markers.main([str(first), str(second)])

    output = capsys.readouterr().out
    assert status == 1
    assert "test_a" in output and "test_b" in output


def test_main_returns_two_when_a_file_is_unreadable(tmp_path, capsys):
    status = check_test_markers.main([str(tmp_path / "missing.py")])

    assert status == 2
    assert "missing.py" in capsys.readouterr().err


def test_main_returns_two_when_a_file_does_not_parse(tmp_path, capsys):
    path = tmp_path / "test_broken.py"
    path.write_text("def test_a(:\n", encoding="utf-8")

    status = check_test_markers.main([str(path)])

    assert status == 2
    assert "test_broken.py" in capsys.readouterr().err


def test_main_with_no_paths_checks_nothing(capsys):
    status = check_test_markers.main([])

    assert status == 0
