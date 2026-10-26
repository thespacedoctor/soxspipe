"""Close the matplotlib figures that a plotting function opened, whether it returns or raises.

The long QC plot functions open a figure in the middle of a body that is hundreds of lines long. Wrapping
that body in ``try``/``finally`` would re-indent every line, so the guard is a decorator instead.

The pyplot figure registry belongs to the process and is shared by all of its threads. The guard is safe when
each process plots on its own, as the pipeline's worker processes do. It takes no lock, so two threads that
plot at the same time can close each other's figures.
"""

import functools
import logging
import sys

__all__ = ["closes_its_figures"]

logger = logging.getLogger(__name__)


def _open_figure_managers():
    """*return the figure managers that pyplot holds open*

    The managers are compared by identity and not by figure number, because pyplot hands out the number of a closed
    figure again.

    **Return:**

    - ``managers`` -- a list of figure managers. Empty if pyplot is not imported, because then no figure is open
    """
    if "matplotlib.pyplot" not in sys.modules:
        return []
    from matplotlib._pylab_helpers import Gcf

    return Gcf.get_all_fig_managers()


def _close_figures_opened_since(managersBefore):
    """*close each figure that is open now and was not open when ``managersBefore`` was taken*

    A figure that fails to close is logged and does not stop the others from closing.

    **Key Arguments:**

    - ``managersBefore`` -- the figure managers that were open before the call
    """
    leftOpen = [manager for manager in _open_figure_managers() if manager not in managersBefore]
    if not leftOpen:
        return
    import matplotlib.pyplot as plt

    for manager in leftOpen:
        try:
            plt.close(manager.canvas.figure)
        except Exception:
            logger.warning("could not close figure %s after a plotting function", manager.num, exc_info=True)


def closes_its_figures(func):
    """*close every figure that the decorated function opens and leaves open, on a normal return and on a raise*

    A figure that the function closes itself is not closed again, and a figure that was open before the call is left
    alone, even when the function closes it and pyplot reuses its number. An error while closing a figure is logged as a
    warning and never replaces the error that the function raised.

    **Key Arguments:**

    - ``func`` -- the plotting function to guard

    **Return:**

    - ``wrapper`` -- the guarded function, with the same signature

    **Usage:**

    ```python
    @closes_its_figures
    def plot_results(self, orderPixelTable):
        import matplotlib.pyplot as plt

        fig = plt.figure()
    ```
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        managersBefore = _open_figure_managers()
        try:
            return func(*args, **kwargs)
        finally:
            _close_figures_opened_since(managersBefore)

    return wrapper
