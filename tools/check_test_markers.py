#!/usr/bin/env python
"""*fail when a test does not declare exactly one execution category*

`tests/conftest.py` refuses to run a suite in which a collected test lacks exactly one of
the `unit`, `integration` and `real_data` markers. That check happens at collection, so a
mistake only shows up once pytest runs. This tool finds the same mistake from the source
text alone, so a pre-commit hook can reject it first. The runtime check stays the
authority: this tool reads a deliberately small set of spellings and refuses to guess at
the rest.

It reads the source with `ast` and never imports the test module or `soxspipe`, because the
pre-commit hook runs with whatever `python` is on `PATH`, which is often too old to import
the package.

A test is a top-level `test_*` function, or a `test_*` method of a `Test*` class. Its
category comes from any of these, spelled `pytest.mark.<category>` after `import pytest`:

- a module-level `pytestmark`, as one mark or a list of marks
- a class decorator, or a class-level `pytestmark`
- a function or method decorator
- a `Test*` class gets the marks of its base classes defined in the same module, however
  far up the chain. A nested class sees the classes of its own enclosing class and of the
  module, as Python resolves a base name

The categories found are collected into a set, as the runtime check does, so one category
declared in two places is fine and two different categories are not. `pytestmark`
assignments are read in order: `=` and an annotated assignment replace the value, and `+=`
extends it.

Three things are reported rather than guessed at:

- any other spelling of a mark: `from pytest import mark`, `import pytest as p`, an alias
  such as `UNIT = pytest.mark.unit`, or `other.mark.unit`. Each is reported as
  `use pytest.mark.<category> (house spelling)`
- a tuple `pytestmark`. pytest wraps a value that is not a list in a list, so a tuple fails
  at collection
- a test with no category, or more than one

These are out of scope, and the runtime check judges them:

- marks on `pytest.param(..., marks=...)`
- a `__test__ = False` class, which pytest skips but this tool still checks
- a base class imported from another module, whose marks this tool cannot see
- a test method inherited from a base class, which is checked only where the base is
  itself a `Test*` class

**Usage:**

```bash
python tools/check_test_markers.py tests/unit/test_example.py
```

Each finding is printed as `path:line: name must declare exactly one execution category`,
or `path:line: <problem>` for a spelling or tuple finding. Exit status is 0 when every test
is marked, 1 when a test is not, and 2 when a file could not be read or parsed.
"""

from __future__ import annotations

import argparse
import ast
import sys
from pathlib import Path
from typing import NamedTuple

# THE EXECUTION CATEGORIES. tests/conftest.py ENFORCES THE SAME SET AT RUNTIME, AND A UNIT
# TEST FAILS IF THE TWO DRIFT APART, SINCE THIS TOOL CANNOT IMPORT THE CONFTEST
CATEGORY_MARKERS = frozenset({"unit", "integration", "real_data"})

TEST_FUNCTION_PREFIX = "test_"
TEST_CLASS_PREFIX = "Test"
PYTESTMARK_NAME = "pytestmark"
PYTEST_MODULE = "pytest"

SPELLING_MESSAGE = "use pytest.mark.<category> (house spelling)"
TUPLE_MESSAGE = "pytestmark must be one mark or a list of marks, not a tuple"

EXIT_CLEAN = 0
EXIT_FINDINGS = 1
EXIT_UNREADABLE = 2

_Definition = ast.FunctionDef | ast.AsyncFunctionDef | ast.ClassDef


class Finding(NamedTuple):
    """*one problem in a test module: a test without one category, or a mark spelling*"""

    path: str
    line: int
    name: str
    detail: str = ""

    def __str__(self) -> str:
        if self.detail:
            return f"{self.path}:{self.line}: {self.detail}"

        return f"{self.path}:{self.line}: {self.name} must declare exactly one execution category"


def _house_mark_name(node: ast.expr) -> str | None:
    """*the mark name in a `pytest.mark.<name>` expression, called or not*

    **Key Arguments:**

    - ``node`` -- a decorator or `pytestmark` element

    **Return:**

    - ``name`` -- the mark name, or None when the node is not spelled `pytest.mark.<name>`
    """
    if isinstance(node, ast.Call):
        node = node.func
    if not isinstance(node, ast.Attribute):
        return None

    owner = node.value
    isHouseSpelling = (
        isinstance(owner, ast.Attribute)
        and owner.attr == "mark"
        and isinstance(owner.value, ast.Name)
        and owner.value.id == PYTEST_MODULE
    )

    return node.attr if isHouseSpelling else None


def _is_foreign_category_mark(node: ast.expr) -> bool:
    """*whether a decorator names a category through anything but `pytest.mark.<category>`*

    **Key Arguments:**

    - ``node`` -- a decorator

    **Return:**

    - ``isForeign`` -- True for `mark.unit`, `p.mark.unit`, `other.mark.unit` and the like
    """
    if isinstance(node, ast.Call):
        node = node.func

    return isinstance(node, ast.Attribute) and node.attr in CATEGORY_MARKERS and _house_mark_name(node) is None


def _decorator_marks(node: _Definition, path: str) -> tuple[set[str], list[Finding]]:
    """*the categories a function or class declares through its decorators, and any spelling findings*

    **Key Arguments:**

    - ``node`` -- the decorated definition
    - ``path`` -- the file path to report findings under

    **Return:**

    - ``categories`` -- the category names found
    - ``findings`` -- one spelling finding per decorator that names a category the wrong way
    """
    categories: set[str] = set()
    findings: list[Finding] = []
    for decorator in node.decorator_list:
        name = _house_mark_name(decorator)
        if name in CATEGORY_MARKERS:
            categories.add(name)
        elif _is_foreign_category_mark(decorator):
            findings.append(Finding(path, decorator.lineno, "", SPELLING_MESSAGE))

    return categories, findings


def _evaluate_pytestmark(value: ast.expr, isExtension: bool) -> tuple[set[str], str | None]:
    """*the categories a `pytestmark` value holds, and the problem with it if it has one*

    **Key Arguments:**

    - ``value`` -- the right-hand side of the assignment
    - ``isExtension`` -- True for `+=`, where a tuple is fine because `list += tuple` is valid

    **Return:**

    - ``categories`` -- the category names found
    - ``problem`` -- the message to report, or None when the value is readable
    """
    if isinstance(value, ast.Tuple) and not isExtension:
        return set(), TUPLE_MESSAGE

    elements = value.elts if isinstance(value, ast.List | ast.Tuple) else [value]
    names = [_house_mark_name(element) for element in elements]
    if any(name is None for name in names):
        return set(), SPELLING_MESSAGE

    return {name for name in names if name in CATEGORY_MARKERS}, None


def _body_marks(body: list[ast.stmt], path: str) -> tuple[set[str], list[Finding]]:
    """*the categories a module or class body assigns to `pytestmark`, read in order*

    **Key Arguments:**

    - ``body`` -- the statements of a module or class
    - ``path`` -- the file path to report findings under

    **Return:**

    - ``categories`` -- the category names the final value holds
    - ``findings`` -- the spelling and tuple findings
    """
    categories: set[str] = set()
    findings: list[Finding] = []
    for statement in body:
        if isinstance(statement, ast.Assign):
            targets, value = statement.targets, statement.value
        elif isinstance(statement, ast.AnnAssign | ast.AugAssign) and statement.value is not None:
            targets, value = [statement.target], statement.value
        else:
            continue
        if not any(isinstance(target, ast.Name) and target.id == PYTESTMARK_NAME for target in targets):
            continue

        isExtension = isinstance(statement, ast.AugAssign)
        found, problem = _evaluate_pytestmark(value, isExtension)
        if problem is not None:
            findings.append(Finding(path, statement.lineno, "", problem))
        categories = categories | found if isExtension else found

    return categories, findings


def _class_marks(
    node: ast.ClassDef, classes: dict[str, ast.ClassDef], path: str, visiting: frozenset[str] = frozenset()
) -> tuple[set[str], list[Finding]]:
    """*the categories a class declares itself or inherits from base classes in the same module*

    **Key Arguments:**

    - ``node`` -- the class
    - ``classes`` -- the module's top-level classes, by name
    - ``path`` -- the file path to report findings under
    - ``visiting`` -- the classes already on the way up the chain, which stops a cycle. Default *frozenset()*

    **Return:**

    - ``categories`` -- the category names found
    - ``findings`` -- the spelling and tuple findings
    """
    categories, findings = _decorator_marks(node, path)
    bodyCategories, bodyFindings = _body_marks(node.body, path)
    categories |= bodyCategories
    findings += bodyFindings

    for base in node.bases:
        if not isinstance(base, ast.Name) or base.id not in classes or base.id in visiting:
            continue
        baseCategories, baseFindings = _class_marks(classes[base.id], classes, path, visiting | {node.name})
        categories |= baseCategories
        findings += baseFindings

    return categories, findings


def _check_body(
    body: list[ast.stmt], scope: set[str], classes: dict[str, ast.ClassDef], path: str, prefix: str
) -> list[Finding]:
    """*the findings in a module or class body*

    **Key Arguments:**

    - ``body`` -- the statements of a module or class
    - ``scope`` -- the categories declared by the enclosing module and classes
    - ``classes`` -- the module's top-level classes, by name. Classes defined in ``body`` are added to them
    - ``path`` -- the file path to report findings under
    - ``prefix`` -- the enclosing class qualifier, e.g. `TestThing.`, or an empty string

    **Return:**

    - ``findings`` -- the offending tests and mark spellings
    """
    findings: list[Finding] = []
    # A BASE IS VISIBLE IF IT IS DEFINED IN THIS BODY OR AT MODULE LEVEL, AS PYTHON RESOLVES IT
    visible = {**classes, **{node.name: node for node in body if isinstance(node, ast.ClassDef)}}

    for statement in body:
        if isinstance(statement, ast.FunctionDef | ast.AsyncFunctionDef):
            if not statement.name.startswith(TEST_FUNCTION_PREFIX):
                continue
            decorated, spellingFindings = _decorator_marks(statement, path)
            findings += spellingFindings
            if len(scope | decorated) != 1:
                findings.append(Finding(path, statement.lineno, f"{prefix}{statement.name}"))
        elif isinstance(statement, ast.ClassDef) and statement.name.startswith(TEST_CLASS_PREFIX):
            classCategories, classFindings = _class_marks(statement, visible, path)
            findings += classFindings
            findings += _check_body(
                statement.body, scope | classCategories, classes, path, f"{prefix}{statement.name}."
            )

    return findings


def _import_findings(tree: ast.Module, path: str) -> list[Finding]:
    """*the imports that give pytest's marks a name other than `pytest.mark`*

    **Key Arguments:**

    - ``tree`` -- the parsed module
    - ``path`` -- the file path to report findings under

    **Return:**

    - ``findings`` -- one spelling finding per `from pytest import mark` or `import pytest as ...`
    """
    findings = []
    for node in ast.walk(tree):
        isMarkImport = (
            isinstance(node, ast.ImportFrom)
            and node.module == PYTEST_MODULE
            and any(alias.name == "mark" for alias in node.names)
        )
        isAliasedImport = isinstance(node, ast.Import) and any(
            alias.name == PYTEST_MODULE and alias.asname not in (None, PYTEST_MODULE) for alias in node.names
        )
        if isMarkImport or isAliasedImport:
            findings.append(Finding(path, node.lineno, "", SPELLING_MESSAGE))

    return findings


def find_unmarked_tests(source: str, path: str) -> list[Finding]:
    """*find every test that does not declare exactly one execution category, and every unreadable mark*

    **Key Arguments:**

    - ``source`` -- the Python source of a test module
    - ``path`` -- the file path to report findings under

    **Return:**

    - ``findings`` -- the problems found, in source order and without repeats

    **Raises:**

    - ``SyntaxError`` -- when the source does not parse
    """
    tree = ast.parse(source, filename=path)
    classes = {node.name: node for node in tree.body if isinstance(node, ast.ClassDef)}
    moduleCategories, findings = _body_marks(tree.body, path)
    findings += _import_findings(tree, path)
    findings += _check_body(tree.body, moduleCategories, classes, path, "")

    return sorted(set(findings), key=lambda finding: (finding.line, finding.detail, finding.name))


def main(argv: list[str] | None = None) -> int:
    """*run the checker from the command line*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments, excluding the program name. Default *None*, i.e. `sys.argv[1:]`

    **Return:**

    - ``status`` -- 0 when every test is marked, 1 when one is not, 2 when a file could not be read or parsed
    """
    parser = argparse.ArgumentParser(description="Reject a test that lacks exactly one execution category.")
    parser.add_argument("paths", nargs="*", type=Path, help="the test files to check")
    arguments = parser.parse_args(argv)

    status = EXIT_CLEAN
    for path in arguments.paths:
        try:
            findings = find_unmarked_tests(path.read_text(encoding="utf-8-sig"), str(path))
        except (OSError, UnicodeError, SyntaxError) as error:
            print(f"{path}: could not check the file: {error}", file=sys.stderr)
            status = EXIT_UNREADABLE
            continue

        for finding in findings:
            print(finding)
        if findings and status == EXIT_CLEAN:
            status = EXIT_FINDINGS

    return status


if __name__ == "__main__":
    sys.exit(main())
