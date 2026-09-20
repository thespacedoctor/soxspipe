#!/usr/bin/env python
"""
*Validate SQL identifiers at the trust boundary, for callers that cannot bind them as parameters*

Author
: David Young

Date Created
: September 20, 2026
"""

import re

# NO IMPORTS FROM OTHER SOXSPIPE MODULES: THE PACKAGE `__init__` IMPORT ORDER
# IS LOAD-BEARING, AND THIS MODULE MUST STAY SAFE TO IMPORT FROM ANYWHERE.
_IDENTIFIER_GRAMMAR = r"[A-Za-z_][A-Za-z0-9_]{0,63}"


class UnsafeSqlIdentifierError(ValueError):
    """*raised when a candidate SQL identifier (a table or column name) fails the safe-identifier grammar*"""


def validate_sql_identifier(identifier, label):
    """*validate a value before it is interpolated into SQL text as an identifier*

    SQLite cannot bind a table or column name as a query parameter, so any
    identifier built from something other than a literal in the source code
    must be checked against a safe grammar before it is interpolated.

    **Key Arguments:**

    - ``identifier`` -- the candidate identifier to validate
    - ``label`` -- a short, human-readable name for what ``identifier`` represents, used in the raised error message

    **Return:**

    - ``identifier`` -- the validated identifier, unchanged

    **Raises:**

    - ``UnsafeSqlIdentifierError`` -- if ``identifier`` is not a string matching ``[A-Za-z_][A-Za-z0-9_]{0,63}``

    **Usage:**

    ```python
    from soxspipe.commonutils.sql_identifiers import validate_sql_identifier
    tableName = validate_sql_identifier(table_name, "table name")
    ```
    """
    if not isinstance(identifier, str) or re.fullmatch(_IDENTIFIER_GRAMMAR, identifier) is None:
        raise UnsafeSqlIdentifierError(f"Unsafe {label}: {identifier!r}")

    return identifier
