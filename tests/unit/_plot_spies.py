"""Matplotlib spies shared by the DY-79 plot characterization tests."""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure


def spy_figures(monkeypatch: pytest.MonkeyPatch) -> list[Figure]:
    """Patch `matplotlib.pyplot.figure`, delegating to the real function while
    recording every `Figure` it returns.

    **Key Arguments:**

    - ``monkeypatch`` -- the pytest monkeypatch fixture

    **Return:**

    - a list that is appended to (in call order) with each `Figure` produced
    """
    original = plt.figure
    figures: list[Figure] = []

    def _wrapper(*args: object, **kwargs: object) -> Figure:
        figure = original(*args, **kwargs)
        figures.append(figure)
        return figure

    monkeypatch.setattr(plt, "figure", _wrapper)
    return figures


def quiet_show(monkeypatch: pytest.MonkeyPatch) -> list[tuple]:
    """Patch `matplotlib.pyplot.show` to a no-op recorder so `show=True` code
    paths can be exercised without a display, without clearing the figure."""
    calls: list[tuple] = []
    monkeypatch.setattr(plt, "show", lambda *args, **kwargs: calls.append((args, kwargs)))
    return calls
