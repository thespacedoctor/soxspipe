#!/usr/bin/env python
"""*check that soxspipe docstrings still agree with the signatures they document*

The checker reads source with `ast` and never imports the module under test, so it
runs without astropy and without any import side effects.

Scope is signature agreement only: argument lists, argument names, and whether a
return value is documented. Whether the prose describes what a function actually
does is out of scope, and no checker can do it.

The docstring format is the bespoke `thespacedoctor`/`fundamentals` Markdown that
`sphinx-autodoc2` renders:

```
*a one-line summary*

**Key Arguments:**

- ``frame`` -- the frame to prepare
- ``save`` -- save the frame to file. Default *False*

**Return:**

- ``frame`` -- the prepared frame
```

**Usage:**

```bash
python tools/check_docstrings.py soxspipe/recipes/soxs_mbias.py
python tools/check_docstrings.py soxspipe
```

Exit status is 0 when nothing is found, 1 when drift is reported, and 2 when a file
could not be read or parsed.
"""

from __future__ import annotations

import argparse
import ast
import re
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import NamedTuple

# SECTION HEADERS, WITH THE COLON-LESS AND PLURAL VARIANTS THE CODEBASE ALSO USES
ARGUMENTS_HEADER_PATTERN = re.compile(r"^\s*\*\*Key Arguments:?\*\*\s*$", re.MULTILINE)
RETURN_HEADER_PATTERN = re.compile(r"^\s*\*\*Returns?:?\*\*\s*$", re.MULTILINE)
ANY_HEADER_PATTERN = re.compile(r"^[ \t]*\*\*[A-Z][^*\n]*\*\*[ \t]*$", re.MULTILINE)

# A DOCUMENTED ARGUMENT BULLET: - ``name`` -- description
ARGUMENT_BULLET_PATTERN = re.compile(r"^\s*-\s+``(\w+)``(?P<remainder>.*)$", re.MULTILINE)

# THE SEPARATOR THE FORMAT REQUIRES BETWEEN AN ARGUMENT NAME AND ITS DESCRIPTION
BULLET_SEPARATOR_PATTERN = re.compile(r"^\s+--\s+\S")

# ARGUMENT NAMES THAT ARE NEVER DOCUMENTED IN THIS CODEBASE
IMPLICIT_ARGUMENT_NAMES = ("self", "cls")

FINDING_KINDS = (
    "undocumented-argument",
    "phantom-argument",
    "no-arguments-section",
    "malformed-argument-bullet",
    "undocumented-return",
    "phantom-return",
    "unparseable",
    "unreadable",
)

# FINDING KINDS THAT MEAN THE CHECKER COULD NOT DO ITS JOB, RATHER THAN DRIFT
TOOL_ERROR_KINDS = ("unparseable", "unreadable")


class Finding(NamedTuple):
    """*one disagreement between a signature and the docstring that documents it*"""

    path: Path
    line: int
    function: str
    kind: str
    detail: str

    def __str__(self) -> str:
        detail = f" ({self.detail})" if self.detail else ""
        return f"{self.path}:{self.line}: {self.function}: {self.kind}{detail}"


def documented_arguments(docstring: str) -> list[str]:
    """*the argument names documented in a docstring's Key Arguments section*

    **Key Arguments:**

    - ``docstring`` -- the raw docstring text

    **Return:**

    - ``names`` -- the documented argument names, in the order they appear
    """
    return [name for name, wellFormed in documented_argument_bullets(docstring)]


def documented_argument_bullets(docstring: str) -> list[tuple[str, bool]]:
    """*the argument bullets in a docstring's Key Arguments section, and whether each is well formed*

    A bullet is well formed when the name is followed by the format's `` -- `` separator and a
    description. A malformed bullet still documents its argument, so it is reported as punctuation
    rather than as a missing argument.

    **Key Arguments:**

    - ``docstring`` -- the raw docstring text

    **Return:**

    - ``bullets`` -- (name, wellFormed) pairs, in the order they appear
    """
    section = _section_body(docstring, ARGUMENTS_HEADER_PATTERN)
    if section is None:
        return []

    return [
        (match.group(1), BULLET_SEPARATOR_PATTERN.match(match.group("remainder")) is not None)
        for match in ARGUMENT_BULLET_PATTERN.finditer(section)
    ]


def has_arguments_section(docstring: str) -> bool:
    """*does the docstring carry a Key Arguments section at all?*

    **Key Arguments:**

    - ``docstring`` -- the raw docstring text

    **Return:**

    - ``present`` -- True if a Key Arguments header is present
    """
    return ARGUMENTS_HEADER_PATTERN.search(docstring) is not None


def has_return_section(docstring: str) -> bool:
    """*does the docstring carry a Return section?*

    **Key Arguments:**

    - ``docstring`` -- the raw docstring text

    **Return:**

    - ``present`` -- True if a Return header is present
    """
    return RETURN_HEADER_PATTERN.search(docstring) is not None


def signature_arguments(node: ast.FunctionDef | ast.AsyncFunctionDef) -> list[str]:
    """*every argument name in a function signature, minus `self` and `cls`*

    **Key Arguments:**

    - ``node`` -- the function definition node

    **Return:**

    - ``names`` -- the argument names, in declaration order
    """
    arguments = node.args
    named = [*arguments.posonlyargs, *arguments.args, *arguments.kwonlyargs]
    names = [argument.arg for argument in named]

    # *args AND **kwargs ARE DOCUMENTED BY THEIR BARE NAMES IN THIS CODEBASE
    if arguments.vararg:
        names.append(arguments.vararg.arg)
    if arguments.kwarg:
        names.append(arguments.kwarg.arg)

    if names and names[0] in IMPLICIT_ARGUMENT_NAMES:
        names = names[1:]

    return names


def returns_a_value(node: ast.FunctionDef | ast.AsyncFunctionDef) -> bool:
    """*does the function return or yield a value of its own?*

    Returns and yields inside nested functions, lambdas, and nested classes belong to
    those definitions, not to this one, so they are not counted.

    **Key Arguments:**

    - ``node`` -- the function definition node

    **Return:**

    - ``returns`` -- True if the function's own body returns or yields a value
    """
    for child in _own_body_nodes(node):
        if isinstance(child, (ast.Yield, ast.YieldFrom)):
            return True
        if isinstance(child, ast.Return) and not _is_none_literal(child.value):
            return True

    return False


def check_function(
    node: ast.FunctionDef | ast.AsyncFunctionDef,
    path: Path,
    qualifiedName: str,
    classDocstring: str | None = None,
) -> list[Finding]:
    """*check one function against the docstring that documents it*

    **Key Arguments:**

    - ``node`` -- the function definition node
    - ``path`` -- the path of the file the function was read from
    - ``qualifiedName`` -- the function name as it should be reported, e.g. `recipe.run`
    - ``classDocstring`` -- the enclosing class docstring, for `__init__` only. Default *None*

    **Return:**

    - ``findings`` -- the disagreements found, empty when the docstring agrees
    """
    docstring = ast.get_docstring(node)

    # __init__ ARGUMENTS ARE DOCUMENTED ON THE CLASS DOCSTRING IN THIS CODEBASE, BUT A
    # CONSTRUCTOR THAT DOCUMENTS ITSELF IS STILL CHECKED RATHER THAN SILENTLY SKIPPED
    if node.name == "__init__" and classDocstring:
        docstring = classDocstring
    if not docstring:
        return []

    findings = []
    arguments = signature_arguments(node)

    if arguments and not has_arguments_section(docstring):
        findings.append(Finding(path, node.lineno, qualifiedName, "no-arguments-section", ", ".join(arguments)))
    else:
        bullets = documented_argument_bullets(docstring)
        documented = [name for name, _ in bullets]
        findings.extend(
            Finding(path, node.lineno, qualifiedName, "undocumented-argument", name) for name in arguments if name not in documented
        )
        findings.extend(
            Finding(path, node.lineno, qualifiedName, "phantom-argument", name) for name in documented if name not in arguments
        )
        findings.extend(
            Finding(path, node.lineno, qualifiedName, "malformed-argument-bullet", name) for name, wellFormed in bullets if not wellFormed
        )

    returns = returns_a_value(node)
    documentsReturn = has_return_section(docstring)
    if returns and not documentsReturn:
        findings.append(Finding(path, node.lineno, qualifiedName, "undocumented-return", ""))
    elif documentsReturn and not returns and node.name != "__init__":
        findings.append(Finding(path, node.lineno, qualifiedName, "phantom-return", ""))

    return findings


def check_source(source: str, path: Path) -> list[Finding]:
    """*check every function in a module's source*

    **Key Arguments:**

    - ``source`` -- the module source text
    - ``path`` -- the path the source was read from, used in the findings

    **Return:**

    - ``findings`` -- the disagreements found, in source order
    """
    try:
        tree = ast.parse(source, filename=str(path))
    except (SyntaxError, ValueError) as error:
        lineNumber = getattr(error, "lineno", 0) or 0
        return [Finding(path, lineNumber, "<module>", "unparseable", str(getattr(error, "msg", error)))]

    findings = []
    for node, qualifiedName, classDocstring in _walk_functions(tree):
        findings.extend(check_function(node, path, qualifiedName, classDocstring))

    return sorted(findings, key=lambda finding: (finding.line, finding.kind, finding.detail))


def check_file(path: Path) -> list[Finding]:
    """*check every function in one Python file*

    **Key Arguments:**

    - ``path`` -- the path of the file to check

    **Return:**

    - ``findings`` -- the disagreements found, in source order
    """
    path = Path(path)
    try:
        source = path.read_text(encoding="utf-8")
    except (OSError, UnicodeError) as error:
        # ONE UNREADABLE FILE MUST NOT ABORT A WHOLE-TREE SCAN
        return [Finding(path, 0, "<module>", "unreadable", str(error))]

    return check_source(source, path)


def find_python_files(paths: list[Path]) -> list[Path]:
    """*expand the given paths into Python files, walking any directory*

    **Key Arguments:**

    - ``paths`` -- the file and directory paths to expand

    **Return:**

    - ``files`` -- the Python files found, sorted and deduplicated
    """
    files = set()
    for path in paths:
        path = Path(path)
        if path.is_dir():
            files.update(path.rglob("*.py"))
        else:
            files.add(path)

    return sorted(files)


def main(argv: list[str] | None = None) -> int:
    """*run the checker from the command line*

    **Key Arguments:**

    - ``argv`` -- the command-line arguments, excluding the program name. Default *None*, i.e. `sys.argv[1:]`

    **Return:**

    - ``status`` -- 0 when nothing was found, 1 when drift was found, 2 when a file could not be read or parsed
    """
    parser = argparse.ArgumentParser(description="Check soxspipe docstrings against their signatures.")
    parser.add_argument("paths", nargs="+", type=Path, help="module files or directories to check")
    parser.add_argument("--quiet", action="store_true", help="report the summary only, not the individual findings")
    arguments = parser.parse_args(argv)

    files = find_python_files(arguments.paths)
    findings = [finding for path in files for finding in check_file(path)]

    toolErrors = [finding for finding in findings if finding.kind in TOOL_ERROR_KINDS]
    for finding in toolErrors:
        print(finding, file=sys.stderr)

    if not arguments.quiet:
        for finding in findings:
            if finding not in toolErrors:
                print(finding)

    _print_summary(findings, len(files))

    if toolErrors:
        return 2

    return 1 if findings else 0


def _print_summary(findings: list[Finding], fileCount: int) -> None:
    """*print the per-kind counts and the headline drift count*

    **Key Arguments:**

    - ``findings`` -- the findings to summarise
    - ``fileCount`` -- the number of files checked
    """
    print(f"\n{len(findings)} findings in {fileCount} files")
    for kind in FINDING_KINDS:
        count = sum(1 for finding in findings if finding.kind == kind)
        if count:
            print(f"  {kind}: {count}")

    driftedFunctions = {
        (finding.path, finding.function) for finding in findings if finding.kind in ("undocumented-argument", "phantom-argument")
    }
    print(f"  drifted functions: {len(driftedFunctions)}")


def _section_body(docstring: str, headerPattern: re.Pattern[str]) -> str | None:
    """*the text between a section header and the next header*

    **Key Arguments:**

    - ``docstring`` -- the raw docstring text
    - ``headerPattern`` -- the compiled pattern matching the section header line

    **Return:**

    - ``body`` -- the section body, or None when the header is absent
    """
    header = headerPattern.search(docstring)
    if not header:
        return None

    remainder = docstring[header.end():]
    nextHeader = ANY_HEADER_PATTERN.search(remainder)

    return remainder[: nextHeader.start()] if nextHeader else remainder


def _own_body_nodes(node: ast.AST) -> Iterator[ast.AST]:
    """*walk a function body, skipping the bodies of nested functions and classes*

    **Key Arguments:**

    - ``node`` -- the function definition node to walk

    **Return:**

    - ``nodes`` -- a generator of the nodes belonging to this function itself
    """
    for child in ast.iter_child_nodes(node):
        if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef, ast.Lambda)):
            continue
        yield child
        yield from _own_body_nodes(child)


def _walk_functions(tree: ast.Module) -> Iterator[tuple[ast.FunctionDef | ast.AsyncFunctionDef, str, str | None]]:
    """*walk every function in a module, carrying its qualified name and enclosing class docstring*

    **Key Arguments:**

    - ``tree`` -- the parsed module

    **Return:**

    - ``functions`` -- a generator of (node, qualifiedName, classDocstring) triples
    """
    yield from _walk_scope(tree, prefix="", classDocstring=None)


def _walk_scope(node: ast.AST, prefix: str, classDocstring: str | None) -> Iterator[tuple[ast.FunctionDef | ast.AsyncFunctionDef, str, str | None]]:
    """*walk one scope, recursing into classes and nested functions*

    **Key Arguments:**

    - ``node`` -- the scope node to walk
    - ``prefix`` -- the dotted name prefix for anything defined in this scope
    - ``classDocstring`` -- the docstring of the immediately enclosing class, if any

    **Return:**

    - ``functions`` -- a generator of (node, qualifiedName, classDocstring) triples
    """
    for child in ast.iter_child_nodes(node):
        if isinstance(child, ast.ClassDef):
            yield from _walk_scope(child, f"{prefix}{child.name}.", ast.get_docstring(child))
        elif isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
            yield child, f"{prefix}{child.name}", classDocstring
            yield from _walk_scope(child, f"{prefix}{child.name}.", None)
        else:
            # A DEFINITION INSIDE AN IF, TRY, FOR, WHILE, WITH OR MATCH BLOCK BELONGS TO THIS SCOPE
            yield from _walk_scope(child, prefix, classDocstring)


def _is_none_literal(node: ast.AST | None) -> bool:
    """*is this a bare `return` or an explicit `return None`?*

    **Key Arguments:**

    - ``node`` -- the returned expression, or None for a bare return

    **Return:**

    - ``isNone`` -- True when nothing of substance is returned
    """
    return node is None or (isinstance(node, ast.Constant) and node.value is None)


if __name__ == "__main__":
    sys.exit(main())
