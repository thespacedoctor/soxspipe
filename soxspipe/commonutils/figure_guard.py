"""Close the matplotlib figures that a plotting function opened, whether it returns or raises.

The long QC plot functions open a figure in the middle of a body that is hundreds of lines long. Wrapping
that body in ``try``/``finally`` would re-indent every line, so the guard is a decorator instead.
"""

import functools
import sys

__all__ = ["closes_its_figures"]


def _open_figure_numbers():
    """*return the numbers of the figures that pyplot holds open*

    **Return:**

    - ``figureNumbers`` -- a set of figure numbers. Empty if pyplot is not imported, because then no figure is open
    """
    pyplot = sys.modules.get("matplotlib.pyplot")
    return set(pyplot.get_fignums()) if pyplot else set()


def closes_its_figures(func):
    """*close every figure that the decorated function opens and leaves open, on a normal return and on a raise*

    A figure that the function closes itself is not closed twice, and a figure that was open before the call is left
    alone.

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
        openBefore = _open_figure_numbers()
        try:
            return func(*args, **kwargs)
        finally:
            leftOpen = _open_figure_numbers() - openBefore
            if leftOpen:
                import matplotlib.pyplot as plt

                for figureNumber in leftOpen:
                    plt.close(figureNumber)

    return wrapper
