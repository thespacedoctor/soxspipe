#!/usr/bin/env python
"""
*collect the warnings and errors an extraction logs so they can be printed on its QC plots*

Author
: David Young

Date Created
: October 14, 2026
"""

import os
from collections.abc import Mapping, Sequence
from typing import Any

os.environ["TERM"] = "vt100"


class extraction_note_logger:
    """
    *a logger proxy that remembers the warnings, errors and notes it is given, and forwards every call but ``note``*

    Only ``warning``, ``error`` and ``note`` are recorded (``note`` is not forwarded). A ``warning`` is forwarded to the
    wrapped logger's ``info``, because the QC plot already shows it as a warning note. Every other call
    (``debug``, ``info``, ``print``, ``warn``, ``critical``, ``exception``, ``log(level, ...)`` ...) and every
    attribute (``handlers`` ...) is forwarded to the wrapped logger without being recorded, so the proxy can stand
    in for it. ``notes`` keeps first-seen order and drops exact duplicates. A message given with arguments is
    recorded as ``msg % args``, as ``logging`` formats it.

    **Key Arguments:**

    - ``log`` -- the logger to wrap

    **Usage:**

    ```python
    from soxspipe.commonutils.extraction_notes import extraction_note_logger
    noteLog = extraction_note_logger(log)
    noteLog.warning("Order skipped: 12")
    noteLog.notes
    # [("WARNING", "Order skipped: 12")]
    ```
    """

    def __init__(self, log: Any) -> None:
        self._log = log
        self._notes = []

    @property
    def notes(self) -> list[tuple[str, str]]:
        """*the notes collected so far, in first-seen order*

        **Return:**

        - ``notes`` -- a new list of ``(level, message)`` pairs
        """
        return list(self._notes)

    def warning(self, msg: object, *args: object, **kwargs: Any) -> Any:
        """*record the message as a warning note, then forward it to the wrapped logger's ``info``*

        **Key Arguments:**

        - ``msg`` -- the message to log
        - ``args`` -- positional arguments for the wrapped logger
        - ``kwargs`` -- keyword arguments for the wrapped logger

        **Return:**

        - ``result`` -- whatever the wrapped logger's ``info`` returns
        """
        self._record("WARNING", msg, args)
        return self._log.info(msg, *args, **kwargs)

    def error(self, msg: object, *args: object, **kwargs: Any) -> Any:
        """*record the message as an error note, then forward it to the wrapped logger*

        **Key Arguments:**

        - ``msg`` -- the message to log
        - ``args`` -- positional arguments for the wrapped logger
        - ``kwargs`` -- keyword arguments for the wrapped logger

        **Return:**

        - ``result`` -- whatever the wrapped logger's ``error`` returns
        """
        self._record("ERROR", msg, args)
        return self._log.error(msg, *args, **kwargs)

    def note(self, msg: object) -> None:
        """*record an informational note for the QC plots; it is not forwarded to the wrapped logger*

        **Key Arguments:**

        - ``msg`` -- the note text
        """
        self._record("NOTE", msg, ())

    def _record(self, level: str, msg: object, args: tuple[object, ...] | Mapping[str, object]) -> None:
        text = str(msg)
        # A SINGLE NON-EMPTY MAPPING IS THE NAMED-FIELD ARGUMENT, AS LOGGING UNWRAPS IT
        if len(args) == 1 and isinstance(args[0], Mapping) and args[0]:
            args = args[0]
        if args:
            try:
                text = text % args
            except (TypeError, ValueError, KeyError):
                # KEEP THE RAW MESSAGE WHEN THE ARGUMENTS DO NOT FIT, AS LOGGING REPORTS THE FAULT ITSELF
                text = str(msg)
        note = (level, text)
        if note not in self._notes:
            self._notes.append(note)

    def __getattr__(self, name: str) -> Any:
        # ONLY CALLED FOR NAMES THE PROXY DOES NOT DEFINE ITSELF. GUARD _log SO COPYING OR UNPICKLING CANNOT RECURSE
        if name == "_log":
            raise AttributeError(name)
        return getattr(self._log, name)


def format_extraction_notes(notes: Sequence[tuple[str, str]], tag: str | None = None) -> list[str]:
    """*turn collected notes into the strings drawn on a QC plot*

    **Key Arguments:**

    - ``notes`` -- the ``(level, message)`` pairs to format
    - ``tag`` -- a label to put in front of every note, e.g. ``A1``. Default *None* (no prefix)

    **Return:**

    - ``lines`` -- one display string per note, e.g. ``[A1] WARNING: Order skipped: 12``

    **Usage:**

    ```python
    from soxspipe.commonutils.extraction_notes import format_extraction_notes
    lines = format_extraction_notes([("WARNING", "Order skipped: 12")], tag="A1")
    ```
    """
    prefix = f"[{tag}] " if tag else ""
    return [f"{prefix}{level}: {message}" for level, message in notes]
