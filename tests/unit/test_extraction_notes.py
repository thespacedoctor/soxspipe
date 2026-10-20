"""Contracts for the extraction note collector and its display formatter."""

from __future__ import annotations

import pytest

from soxspipe.commonutils.extraction_notes import (
    extraction_note_logger,
    format_extraction_notes,
)

pytestmark = pytest.mark.unit


def test_collector_forwards_warning_and_error_to_the_wrapped_logger(log: object) -> None:
    collector = extraction_note_logger(log)

    collector.warning("something odd")
    collector.error("something broken")

    assert log.messages == [("warning", "something odd"), ("error", "something broken")]


def test_collector_records_warning_and_error_in_the_order_they_arrive(log: object) -> None:
    collector = extraction_note_logger(log)

    collector.warning("first")
    collector.error("second")
    collector.warning("third")

    assert collector.notes == [("WARNING", "first"), ("ERROR", "second"), ("WARNING", "third")]


def test_collector_records_an_informational_note_without_forwarding_it(log: object) -> None:
    collector = extraction_note_logger(log)

    collector.warning("odd")
    collector.note("tuning skipped")
    collector.note("tuning skipped")

    assert collector.notes == [("WARNING", "odd"), ("NOTE", "tuning skipped")]
    assert log.messages == [("warning", "odd")]


def test_collector_ignores_debug_info_and_print_but_still_forwards_them(log: object) -> None:
    collector = extraction_note_logger(log)

    collector.debug("d")
    collector.info("i")
    collector.print("p")

    assert collector.notes == []
    assert log.messages == [("debug", "d"), ("info", "i"), ("print", "p")]


def test_collector_drops_exact_duplicates_but_forwards_every_call(log: object) -> None:
    collector = extraction_note_logger(log)

    collector.warning("same")
    collector.warning("same")
    collector.error("same")

    assert collector.notes == [("WARNING", "same"), ("ERROR", "same")]
    assert len(log.messages) == 3


def test_collector_records_the_text_of_a_non_string_message(log: object) -> None:
    collector = extraction_note_logger(log)

    collector.warning(ValueError("bad value"))

    assert collector.notes == [("WARNING", "bad value")]


def test_collector_exposes_the_wrapped_loggers_other_attributes(log: object) -> None:
    collector = extraction_note_logger(log)

    assert collector.handlers is log.handlers


def test_collector_notes_are_a_copy_so_callers_cannot_mutate_the_record(log: object) -> None:
    collector = extraction_note_logger(log)
    collector.warning("kept")

    collector.notes.append(("ERROR", "injected"))

    assert collector.notes == [("WARNING", "kept")]


def test_format_extraction_notes_prefixes_each_note_with_its_level() -> None:
    notes = [("WARNING", "Order skipped: 12"), ("ERROR", "spill")]

    assert format_extraction_notes(notes) == ["WARNING: Order skipped: 12", "ERROR: spill"]


def test_format_extraction_notes_adds_the_tag_prefix_when_given() -> None:
    notes = [("WARNING", "Order skipped: 12")]

    assert format_extraction_notes(notes, tag="A1") == ["[A1] WARNING: Order skipped: 12"]


def test_format_extraction_notes_returns_empty_list_for_no_notes() -> None:
    assert format_extraction_notes([], tag="A1") == []


class _FormattingLogger:
    """A module-level logger double that accepts logging-style arguments and can be pickled."""

    handlers: tuple = ()

    def warning(self, msg: object, *args: object, **kwargs: object) -> None:
        return None

    def error(self, msg: object, *args: object, **kwargs: object) -> None:
        return None


def test_collector_records_the_message_formatted_with_its_arguments_like_logging() -> None:
    collector = extraction_note_logger(_FormattingLogger())

    collector.warning("Order %s skipped (%d slices)", 12, 40)
    collector.error("spill in %s", "order 11")

    assert collector.notes == [("WARNING", "Order 12 skipped (40 slices)"), ("ERROR", "spill in order 11")]


def test_collector_keeps_a_percent_sign_literal_when_no_arguments_are_given() -> None:
    collector = extraction_note_logger(_FormattingLogger())

    collector.warning("50% of slices")

    assert collector.notes == [("WARNING", "50% of slices")]


def test_collector_forwards_the_original_message_and_arguments_unformatted(log: object) -> None:
    received: list[tuple[object, ...]] = []

    class ArgLogger:
        def warning(self, msg: object, *args: object, **kwargs: object) -> None:
            received.append((msg, args, kwargs))

    extraction_note_logger(ArgLogger()).warning("n=%d", 3, stacklevel=2)

    assert received == [("n=%d", (3,), {"stacklevel": 2})]


def test_collector_does_not_record_the_other_levels_but_forwards_them(log: object) -> None:
    collector = extraction_note_logger(log)

    collector.critical("fatal")

    assert collector.notes == []
    assert log.messages == [("critical", "fatal")]


def test_collector_survives_copy_and_pickle_round_trips_without_recursing() -> None:
    import copy
    import pickle

    collector = extraction_note_logger(_FormattingLogger())
    collector.warning("kept")

    copied = copy.copy(collector)
    deepCopied = copy.deepcopy(collector)
    unpickled = pickle.loads(pickle.dumps(collector))  # noqa: S301

    assert copied.notes == deepCopied.notes == unpickled.notes == [("WARNING", "kept")]


def test_collector_without_a_wrapped_logger_raises_attribute_error_instead_of_recursing() -> None:
    bare = object.__new__(extraction_note_logger)

    with pytest.raises(AttributeError):
        getattr(bare, "anything")  # noqa: B009


def test_collector_unwraps_a_single_mapping_argument_like_logging() -> None:
    collector = extraction_note_logger(_FormattingLogger())

    collector.warning("order %(order)s", {"order": 12})

    assert collector.notes == [("WARNING", "order 12")]


def test_collector_falls_back_to_the_raw_text_when_a_mapping_key_is_missing() -> None:
    collector = extraction_note_logger(_FormattingLogger())

    collector.error("order %(order)s", {"other": 1})

    assert collector.notes == [("ERROR", "order %(order)s")]
