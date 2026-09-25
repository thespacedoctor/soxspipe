"""Tests for polynomial-order override validation."""

import pytest

from soxspipe.recipes.poly_orders import validate_poly_orders


@pytest.mark.parametrize(
    ("value", "digit_count", "expected"),
    [(3454, 4, 3454), ("3454", 4, 3454), (34, 2, 34), (345435, 6, 345435)],
)
def test_valid_poly_orders(value, digit_count, expected):
    assert validate_poly_orders(value, digit_count) == expected


@pytest.mark.parametrize(
    ("value", "digit_count"),
    [
        (7, 4),
        (3454.9, 4),
        ("345", 4),
        ("34545", 4),
        ("34.5", 4),
        (True, 4),
        (0, 4),
        (None, 4),
        (-345, 4),
        (3454, 2),
        (3454, 6),
    ],
)
def test_invalid_poly_orders(value, digit_count):
    with pytest.raises(
        TypeError,
        match=rf"THE poly VALUE NEEDS TO BE A {digit_count} DIGIT INTEGER",
    ):
        validate_poly_orders(value, digit_count)
