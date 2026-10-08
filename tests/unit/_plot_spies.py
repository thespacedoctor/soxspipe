"""Matplotlib spies shared by the DY-79 plot characterization tests."""

from __future__ import annotations

from collections.abc import Callable

import matplotlib.pyplot as plt
import pytest
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.transforms import Bbox


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


def probe_on_savefig(monkeypatch: pytest.MonkeyPatch, probe: Callable[[Figure], object]) -> list[object]:
    """Patch `matplotlib.pyplot.savefig` to record `probe(figure)` for the current figure, then save as normal.

    QC plot functions clear or close their figure after saving, so the figure must be probed at save time.

    **Key Arguments:**

    - ``monkeypatch`` -- the pytest monkeypatch fixture
    - ``probe`` -- called with the current figure each time `savefig` runs

    **Return:**

    - a list that is appended to (in call order) with each probe result
    """
    probes: list[object] = []
    realSavefig = plt.savefig

    def recorder(*args: object, **kwargs: object) -> None:
        probes.append(probe(plt.gcf()))
        realSavefig(*args, **kwargs)

    monkeypatch.setattr(plt, "savefig", recorder)
    return probes


def image_interpolations(figure: Figure) -> list[str]:
    """Return the interpolation setting of every image on the figure, in axis order."""
    return [image.get_interpolation() for ax in figure.axes for image in ax.get_images()]


def table_extents(figure: Figure) -> list[Bbox]:
    """Draw the figure on an Agg canvas and return the window extent of every table, in axis order."""
    renderer = FigureCanvasAgg(figure).get_renderer()
    figure.draw(renderer)
    return [table.get_window_extent(renderer) for ax in figure.axes for table in ax.tables]
