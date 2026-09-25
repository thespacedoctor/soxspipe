"""Validation helpers for recipe polynomial-order overrides."""


def validate_poly_orders(poly_orders, digit_count):
    """Return a validated polynomial-order override as an integer.

    Command-line overrides arrive as strings, while callers of the Python API
    commonly provide integers. No other type is accepted: in particular,
    converting floats would silently discard their fractional part.
    """
    message = f"THE poly VALUE NEEDS TO BE A {digit_count} DIGIT INTEGER"
    if isinstance(poly_orders, bool) or not isinstance(poly_orders, (int, str)):
        raise TypeError(message)

    value = str(poly_orders)
    if len(value) != digit_count or not value.isascii() or not value.isdigit():
        raise TypeError(message)

    return int(value)
