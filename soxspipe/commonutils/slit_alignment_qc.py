"""QC plot of the slit shifts that slit alignment applied to a group of frames.

The plot shows the shift of every frame against the first frame in time, one line per order (or one line for the
frame, when a frame has a single shift). A shift that was not measured and was taken from other orders or frames is
marked with a cross.
"""

from soxspipe.commonutils.slit_alignment import WHOLE_FRAME

__all__ = ["plot_slit_alignment_qc"]

PLOT_SIZE_INCHES = (8, 4.5)
FALLBACK_COLOUR = "tab:red"


def plot_slit_alignment_qc(shifts, title, filePath):
    """*plot the slit shift of every frame, with a cross on each shift that was not measured*

    **Key Arguments:**

    - ``shifts`` -- the shift table of ``slit_shifts``
    - ``title`` -- the figure title
    - ``filePath`` -- the PDF path to write

    **Return:**

    - ``filePath`` -- the path the plot was written to

    **Usage:**

    ```python
    from soxspipe.commonutils.slit_alignment_qc import plot_slit_alignment_qc
    plot_slit_alignment_qc(result.shifts, "Slit alignment (A frames)", "/path/to/plot.pdf")
    ```
    """
    import matplotlib.pyplot as plt

    frames = shifts.sort_values("mjd", kind="stable")["frame"].drop_duplicates().tolist()
    framePosition = {frame: index for index, frame in enumerate(frames)}
    fig, axis = plt.subplots(figsize=PLOT_SIZE_INCHES)
    try:
        for order, rows in shifts.groupby("order"):
            label = "frame" if order == WHOLE_FRAME else f"order {order}"
            rows = rows.assign(x=rows["frame"].map(framePosition)).sort_values("x")
            axis.plot(rows["x"], rows["shift_arcsec"], marker="o", markersize=4, label=label)
        fallback = shifts.loc[shifts["is_fallback"]]
        if len(fallback):
            axis.plot(
                fallback["frame"].map(framePosition),
                fallback["shift_arcsec"],
                linestyle="none",
                marker="x",
                markersize=10,
                color=FALLBACK_COLOUR,
                label="fallback",
            )
        axis.set_xticks(range(len(frames)))
        axis.set_xticklabels(frames, rotation=45, ha="right", fontsize=7)
        axis.set_ylabel("Slit shift from the first frame (arcsec)")
        axis.axhline(0.0, color="grey", linewidth=0.5)
        axis.legend(fontsize=7)
        axis.set_title(title, fontsize=10)
        fig.tight_layout()
        plt.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
    finally:
        plt.close(fig)
    return filePath
