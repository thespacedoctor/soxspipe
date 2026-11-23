"""*the A and B frames of a nodding SOF: the split, the balance check and the per-order rule*

The nod recipe and the stare recipe, when it reduces the frames of a nodding SOF (DY-2194), share these functions.
"""


def split_ab_frames(allObjectFrames, allFilenames, kw):
    """*split the science frames into A and B frames by the sign of their cumulative nodding offset*

    A frame with a positive offset is an A frame. Every other frame, a zero offset included, is a B frame.

    **Key Arguments:**

    - ``allObjectFrames`` -- the science frames
    - ``allFilenames`` -- the science frames' file names, in the same order
    - ``kw`` -- the keyword lookup

    **Return:**

    - ``allFrameA``, ``allFrameB``, ``allFrameAOffsets``, ``allFrameBOffsets``, ``allFrameANames``,
      ``allFrameBNames`` -- the A and B frames, their offsets and their file names
    """
    allFrameA, allFrameB = [], []
    allFrameAOffsets, allFrameBOffsets = [], []
    allFrameANames, allFrameBNames = [], []

    # CUMOFF Y IS THE OFFSET IN THE Y DIRECTION OF THE NODDING SEQUENCE. POSITIVE A, NEGATIVE B
    for frame, filename in zip(allObjectFrames, allFilenames, strict=False):
        offset = frame.header[kw("NOD_CUMULATIVE_OFFSETY")]
        if offset > 0:
            allFrameAOffsets.append(offset)
            allFrameA.append(frame)
            allFrameANames.append(filename)
        else:
            allFrameBOffsets.append(offset)
            allFrameB.append(frame)
            allFrameBNames.append(filename)

    return (
        allFrameA,
        allFrameB,
        allFrameAOffsets,
        allFrameBOffsets,
        allFrameANames,
        allFrameBNames,
    )


def check_ab_balance(log, allFrameAOffsets, allFrameBOffsets):
    """*refuse an unbalanced or A-less set of nodding frames and announce the cycles to process*

    **Key Arguments:**

    - ``log`` -- logger
    - ``allFrameAOffsets`` -- the cumulative offsets of the A frames
    - ``allFrameBOffsets`` -- the cumulative offsets of the B frames

    **Return:**

    - ``uniqueOffsets`` -- the distinct A-frame offsets

    **Raises:**

    - ``Exception`` -- when the numbers of A and B frames differ, or there is no A frame
    """
    uniqueOffsets = list(set(allFrameAOffsets))

    if len(allFrameAOffsets) != len(allFrameBOffsets):
        error = (
            f"Found {len(allFrameAOffsets)} A frames and {len(allFrameBOffsets)} B frames. "
            "The number of A and B frames must be the same for nodding reductions."
        )
        log.error(error)
        raise Exception(error)

    if len(uniqueOffsets) == 0:
        error = (
            "Did not find frames with a positive offset. "
            "Please check the `NOD_CUMULATIVE_OFFSETY` header keyword in the providing nodding frames."
        )
        log.error(error)
        raise Exception(error)

    s = "S" if len(uniqueOffsets) > 1 else ""
    log.print(
        f"# PROCESSING {len(allFrameAOffsets)} AB NODDING CYCLES WITH {len(uniqueOffsets)} UNIQUE "
        f"PAIR{s} OF OFFSET LOCATIONS"
    )

    return uniqueOffsets


def aligns_each_order(instrument, arm):
    """*decide whether each order of a frame gets its own slit shift*

    Only SOXS VIS does. SOXS NIR and every X-Shooter arm get one shift per frame.

    **Key Arguments:**

    - ``instrument`` -- the instrument name
    - ``arm`` -- the arm name

    **Return:**

    - ``isPerOrder`` -- True for SOXS VIS
    """
    return instrument.upper() == "SOXS" and arm.upper() == "VIS"
