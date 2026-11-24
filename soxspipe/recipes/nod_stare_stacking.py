"""
*Stack the frames of a nodding SOF that the stare recipe reduces in stare mode (DY-2194)*

Author
: David Young

Date Created
: November 16, 2026
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, NamedTuple

from soxspipe.commonutils.spectrum_names import SPECTRUM_TABLE_LABEL

if TYPE_CHECKING:
    import numpy as np

    from soxspipe.commonutils.slit_alignment import two_d_map_images


class reduction_inputs(NamedTuple):
    """The calibration inputs of a stare reduction: the calibration frames and the paths of the calibration tables.

    ``masterBias``, ``dark`` and ``masterFlat`` are a frame or False. ``responseFunctionPath`` is a path or False.
    """

    masterBias: object
    dark: object
    masterFlat: object
    orderTablePath: str
    dispMap: str
    twoDMap: str
    responseFunctionPath: object


class response_inputs(NamedTuple):
    """What the not-flattened extraction of a standard star needs, kept until the response curve is built."""

    unflattened: Any
    skyModel: Any
    twoDMap: str
    dispMap: str
    flatPass: Any


class stack_extraction(NamedTuple):
    """The flat-fielded extraction of one stack, before flux calibration or the response curve.

    ``responseInputs`` is set when the standard's not-flattened extraction is still to run, and ``notFlatPath`` and
    ``notFlatSpectrum`` when it already ran. All three are None for a science stack.
    """

    productPath: str | None
    combinedObject: Any
    mergedSpectrum: Any
    orderJoins: Any
    extractionPath: str | bool
    responseInputs: response_inputs | None
    notFlatPath: str | None
    notFlatSpectrum: Any
    traceFallback: bool
    slitEdgeSpill: bool
    extractionNotes: list[tuple[str, str]]


class position_stack(NamedTuple):
    """The stacks of the frames of one nod position, with the pixels of the slit that not every frame covers."""

    position: str
    objectFrames: list[Any]
    notFlattened: Any
    flattened: Any
    uncovered: np.ndarray


class nod_stare_stacking:
    """
    Mixin of the methods of `soxs_stare` that align and stack the A and B frames of a nodding SOF before extraction.

    The object sits at one slit position in the A frames and at another in the B frames, and no frame is subtracted
    from another. Each frame is detrended and the frames of each nod position are aligned to each other along the slit.
    The A frames and the B frames are then stacked apart, so each stack holds the object once, at its own slit position.
    The B frames are never shifted onto the A frames, which would leave a slit range that only the A frames cover.
    """

    # THE RAW FILE NAMES OF THE OBJECT FRAMES OF THE RUN, IN THE ORDER THE FRAMES WERE READ. THE SLIT ALIGNMENT NAMES
    # ITS FRAMES WITH THEM
    objectFrameNames = ()

    def _reduce_nod_frames_as_stare(self, allObjectFrames, rawFrames):
        """*align and stack the A and B frames apart, then reduce the stacks with the stare tail*

        **Key Arguments:**

        - ``allObjectFrames`` -- the A and B object frames, in observation order
        - ``rawFrames`` -- the raw file names the product header lists, or *False* to list every raw frame of the SOF

        **Return:**

        - ``productPath`` -- the path to the product
        - ``forceFailure`` -- the failure the tail reports
        """
        framesA, framesB, namesA, namesB = self._split_nod_frames(allObjectFrames)
        inputs = self._read_reduction_inputs(self.kw, self.arm)
        stacks = list(self._stack_nod_positions(framesA, framesB, namesA, namesB, inputs))
        return self._reduce_nod_positions(stacks, allObjectFrames, rawFrames, inputs)

    def _split_nod_frames(self, allObjectFrames):
        """*split the object frames into the A and B frames and refuse a set that is not balanced*

        **Key Arguments:**

        - ``allObjectFrames`` -- the object frames, in observation order

        **Return:**

        - ``framesA`` -- the A frames
        - ``framesB`` -- the B frames
        - ``namesA`` -- the raw file names of the A frames
        - ``namesB`` -- the raw file names of the B frames

        **Raises:**

        - ``ValueError`` -- when ``objectFrameNames`` does not hold one name for each frame
        """
        from soxspipe.commonutils.nod_frames import check_ab_balance, split_ab_frames

        names = list(self.objectFrameNames)
        if len(names) != len(allObjectFrames):
            raise ValueError(
                f"The recipe holds {len(allObjectFrames)} object frames but {len(names)} raw file names, so the "
                "frames cannot be named for the slit alignment."
            )
        framesA, framesB, offsetsA, offsetsB, namesA, namesB = split_ab_frames(allObjectFrames, names, self.kw)
        check_ab_balance(self.log, offsetsA, offsetsB)
        return framesA, framesB, namesA, namesB

    def _stack_nod_positions(
        self,
        framesA: list[Any],
        framesB: list[Any],
        namesA: list[str],
        namesB: list[str],
        inputs: reduction_inputs,
    ) -> tuple[position_stack, position_stack]:
        """*detrend and align the A and B frames, then stack the frames of each nod position on their own*

        Each frame is detrended with the master bias or dark and the scattered-light fit, but not flat-fielded. The A
        frames are aligned to each other and the B frames to each other, along the slit. Nothing is shifted between
        the two positions.

        **Key Arguments:**

        - ``framesA`` -- the A frames
        - ``framesB`` -- the B frames
        - ``namesA`` -- the raw file names of the A frames
        - ``namesB`` -- the raw file names of the B frames
        - ``inputs`` -- the ``reduction_inputs`` of `_read_reduction_inputs`

        **Return:**

        - ``stackA`` -- the ``position_stack`` of the A frames
        - ``stackB`` -- the ``position_stack`` of the B frames

        **Raises:**

        - ``ValueError`` -- when the alignment shifts of a position leave no slit range that every frame of it covers
          in any order
        """
        detrendedA = self._detrend_nod_frames(framesA, inputs.masterBias, inputs.dark, inputs.orderTablePath)
        detrendedB = self._detrend_nod_frames(framesB, inputs.masterBias, inputs.dark, inputs.orderTablePath)
        if len(detrendedA) + len(detrendedB) > 1:
            # EACH FRAME'S SCATTERED-LIGHT FIT ADDED A ROW FOR THE ONE PLOT FILE THE LAST FRAME LEFT
            self.products = self.products.drop_duplicates(subset=["product_label", "file_name"], keep="last")

        alignedA, alignedB, (shiftsA, shiftsB), mapImages = self._align_nod_frames(
            detrendedA, detrendedB, namesA, namesB, inputs.twoDMap, inputs.dispMap
        )
        return (
            self._stack_position("A", framesA, alignedA, shiftsA, mapImages, inputs.masterFlat),
            self._stack_position("B", framesB, alignedB, shiftsB, mapImages, inputs.masterFlat),
        )

    def _stack_position(
        self,
        position: str,
        objectFrames: list[Any],
        alignedFrames: list[Any],
        alignmentShifts: list[dict[int, float]],
        mapImages: two_d_map_images | None,
        masterFlat: Any,
    ) -> position_stack:
        """*stack the aligned frames of one nod position, not flat-fielded and flat-fielded*

        The aligned frames are released once they are stacked. A frame is flat-fielded before its stack, as the flat is
        fixed on the detector.

        **Key Arguments:**

        - ``position`` -- ``A`` or ``B``
        - ``objectFrames`` -- the raw object frames of the position
        - ``alignedFrames`` -- the detrended frames of the position, aligned to each other. The list is emptied
        - ``alignmentShifts`` -- for each frame, its ``{order: shift in arcsec}`` from the alignment
        - ``mapImages`` -- the ``two_d_map_images`` the alignment read, or *None* when it read none
        - ``masterFlat`` -- the master flat, or False when no flat is used

        **Return:**

        - ``stack`` -- the ``position_stack`` of the position. Its ``flattened`` is its ``notFlattened`` when no flat is
          used
        """
        uncovered = self._position_uncovered_mask(position, alignedFrames, alignmentShifts, mapImages)
        notFlattened = self._stack_aligned_frames(alignedFrames)
        if masterFlat is False:
            flattened = notFlattened
            alignedFrames.clear()
        else:
            flatFielded = self._flat_field_nod_frames(alignedFrames, masterFlat)
            alignedFrames.clear()
            flattened = self._stack_aligned_frames(flatFielded)
        return position_stack(position, objectFrames, notFlattened, flattened, uncovered)

    def _stack_aligned_frames(self, frames: list[Any]) -> Any:
        """*clip and stack the frames of a nod position, leaving each frame's masked pixels out*"""
        return self.clip_and_stack(
            frames=frames,
            recipe="soxs_stare",
            ignore_input_masks=False,
            post_stack_clipping=False,
        )

    def _position_uncovered_mask(
        self,
        position: str,
        frames: list[Any],
        alignmentShifts: list[dict[int, float]],
        mapImages: two_d_map_images | None,
    ) -> np.ndarray:
        """*the pixels of the slit that not every frame of a nod position covers (DY-2290)*

        The frames of a position are aligned to each other, so the shift of a frame can leave slit pixels that it does
        not reach. The sky changes between exposures, so the sky residual of the stack steps where a frame ends, and
        the step reads as a second source. The caller leaves the returned pixels out of the sky fit and the extraction
        and flags them in the products.

        **Key Arguments:**

        - ``position`` -- ``A`` or ``B``
        - ``frames`` -- the aligned frames of the position
        - ``alignmentShifts`` -- for each frame, its ``{order: shift in arcsec}`` from the alignment
        - ``mapImages`` -- the ``two_d_map_images`` the alignment read, or *None* when it read none, so no frame was
          shifted

        **Return:**

        - ``uncovered`` -- boolean image, True for a pixel that not every frame of the position covers

        **Raises:**

        - ``ValueError`` -- when the shifts leave no slit range that every frame covers in any order
        """
        import numpy as np

        from soxspipe.commonutils.slit_coverage import alignment_coverage_mask, orders_without_coverage

        if mapImages is None:
            return np.zeros(frames[0].data.shape, dtype=bool)
        uncovered = alignment_coverage_mask(alignmentShifts, mapImages)
        emptyOrders = orders_without_coverage(uncovered, mapImages)
        if emptyOrders and len(emptyOrders) == len(mapImages.orders):
            raise ValueError(
                f"The alignment shifts leave no slit range that every {position} frame covers in any order, so the "
                f"{position} frames cannot be stacked in stare mode."
            )
        if emptyOrders:
            orderList = ", ".join(str(order) for order in emptyOrders)
            self.log.warning(
                f"The alignment shifts leave no slit range that every {position} frame covers in order {orderList}, "
                f"so it is left out of the sky model and the extraction of the {position} frames."
            )
        self.log.info(
            f"{int(uncovered.sum())} pixels that not every {position} frame covers are left out of the sky model and "
            "the extraction."
        )
        return uncovered

    def _extract_position(
        self, stack: position_stack, inputs: reduction_inputs, rawFrames: list[str] | bool
    ) -> stack_extraction:
        """*sky-subtract and extract the stack of one nod position, under the names of that position*

        The SOF name carries the position (``<SOF>_A``) while the position is reduced, so the sky products, the sky QC
        plots and the extraction products of the A and B stacks do not overwrite each other. The sky fallback edits the
        recipe settings in place, so the position runs on a copy. The SOF name, filename template, recipe settings and
        sky flag are restored afterwards, whatever happens. The Horne extraction registers a merged-table product row
        that the combined spectrum registers again, so the position's own row is labelled with the position.

        **Key Arguments:**

        - ``stack`` -- the ``position_stack`` of the position
        - ``inputs`` -- the ``reduction_inputs`` of `_read_reduction_inputs`
        - ``rawFrames`` -- the raw file names the product header lists, or *False* to list every raw frame of the SOF

        **Return:**

        - ``extraction`` -- the ``stack_extraction`` of the position
        """
        import copy

        sofName, filenameTemplate = self.sofName, self.filenameTemplate
        recipeSettings, subtractSky = self.recipeSettings, self.subtractSky
        firstNewRow = len(self.products)
        try:
            self.sofName = f"{sofName}_{stack.position}"
            self.filenameTemplate = f"{self.productStem}.fits"
            self.recipeSettings = copy.deepcopy(recipeSettings)
            extraction = self._extract_stack(
                stack.notFlattened, stack.objectFrames, rawFrames, inputs, stack.flattened, stack.uncovered
            )
            extraction = self._extract_not_flattened_of_standard(extraction)
        finally:
            self.sofName, self.filenameTemplate = sofName, filenameTemplate
            self.recipeSettings, self.subtractSky = recipeSettings, subtractSky
        self._relabel_merged_table_rows_from(firstNewRow, stack.position)
        return extraction

    def _extract_not_flattened_of_standard(self, extraction: stack_extraction) -> stack_extraction:
        """*run the not-flattened extraction of a standard star, with the flat pass of the same position*

        Updates ``self.qc`` and ``self.products``. The extraction's inputs are released once it has run.

        **Key Arguments:**

        - ``extraction`` -- the ``stack_extraction`` of the flat-fielded pass

        **Return:**

        - ``extraction`` -- the same extraction with the path and spectrum of the not-flattened pass, or the extraction
          as it was when it has no ``responseInputs`` (a science frame, a standard with a response table, or one with
          no object trace)
        """
        if extraction.responseInputs is None:
            return extraction
        notFlatPath, notFlatSpectrum = self._extract_not_flattened(*extraction.responseInputs)
        return extraction._replace(responseInputs=None, notFlatPath=notFlatPath, notFlatSpectrum=notFlatSpectrum)

    def _relabel_merged_table_rows_from(self, firstRow: int, position: str) -> None:
        """*relabel every merged-table product row from ``firstRow`` on with a nod position*

        Updates ``self.products``.

        **Key Arguments:**

        - ``firstRow`` -- the position in ``self.products`` of the first row the extraction added
        - ``position`` -- ``A`` or ``B``
        """
        isMerged = (self.products["product_label"] == SPECTRUM_TABLE_LABEL) & (self.products.index >= firstRow)
        self.products.loc[isMerged, "product_label"] = f"{SPECTRUM_TABLE_LABEL}_{position}"

    def _reduce_nod_positions(
        self,
        stacks: list[position_stack],
        allObjectFrames: list[Any],
        rawFrames: list[str] | bool,
        inputs: reduction_inputs,
    ) -> tuple[str, bool | str]:
        """*extract the A and B stacks, combine the two extractions, then flux-calibrate and write the product*

        The A stack is extracted first, then the B stack. Each stack leaves ``stacks`` as it is extracted, so its data
        are released once its position is done. The two extractions are combined as the nod recipe combines its
        cycles, and flux calibration, the response curve of a standard, the ASCII spectrum and the plots take the
        combined spectrum.

        **Key Arguments:**

        - ``stacks`` -- the ``position_stack`` list of the A frames and of the B frames. It is emptied
        - ``allObjectFrames`` -- the A and B object frames, in observation order
        - ``rawFrames`` -- the raw file names the product header lists, or *False* to list every raw frame of the SOF
        - ``inputs`` -- the ``reduction_inputs`` of `_read_reduction_inputs`

        **Return:**

        - ``productPath`` -- the path to the combined spectrum
        - ``forceFailure`` -- True when the response curve or a slit-edge spill fails the recipe, or the failure message
          when a standard star has no object trace
        """
        extractions = {}
        while stacks:
            stack = stacks.pop(0)
            extractions[stack.position] = self._extract_position_keeping_its_qc_only_with_a_trace(
                stack, inputs, rawFrames
            )
        combined = self._combine_extractions(extractions)
        forceFailure = self._finish_extraction(combined, allObjectFrames, inputs.responseFunctionPath)
        return combined.productPath, forceFailure

    def _combine_not_flattened(self, extractions: dict[str, stack_extraction], orderJoins: Any) -> str | None:
        """*combine the not-flattened extractions of the nod positions of a standard star*

        **Key Arguments:**

        - ``extractions`` -- the ``stack_extraction`` of each nod position, by position
        - ``orderJoins`` -- the order joins of the combined spectrum

        **Return:**

        - ``notFlatPath`` -- the path to the combined not-flattened spectrum, or None when the positions made no
          not-flattened extraction
        """
        spectra = [extraction.notFlatSpectrum for extraction in extractions.values()]
        if all(spectrum is None for spectrum in spectra):
            return None
        _, notFlatPath = self.stack_extractions(spectra, notFlattened=True, orderJoins=orderJoins)
        return notFlatPath

    def _extract_position_keeping_its_qc_only_with_a_trace(
        self, stack: position_stack, inputs: reduction_inputs, rawFrames: list[str] | bool
    ) -> stack_extraction:
        """*extract a nod position, and drop the QC rows it added when it found no object trace*

        A position with no trace is left out of the combined spectrum, so its QC rows (the boxcar ``EXTRACT METHOD``,
        the sky and slit QCs of a frame with no object) must not describe that spectrum. Updates ``self.qc``.

        **Key Arguments:**

        - ``stack`` -- the ``position_stack`` of the position
        - ``inputs`` -- the ``reduction_inputs`` of `_read_reduction_inputs`
        - ``rawFrames`` -- the raw file names the product header lists, or *False* to list every raw frame of the SOF

        **Return:**

        - ``extraction`` -- the ``stack_extraction`` of the position
        """
        qcBefore = self.qc.copy()
        extraction = self._extract_position(stack, inputs, rawFrames)
        if extraction.traceFallback:
            self.qc = qcBefore
        return extraction

    def _combine_extractions(self, extractions: dict[str, stack_extraction]) -> stack_extraction:
        """*combine the extractions of the nod positions into one spectrum, as the nod recipe combines its cycles*

        A position that found no object trace gave a boxcar extraction, which is left out: the spectrum is the other
        position's alone, and a warning says so. When neither position found a trace the combination fails, as the nod
        recipe does. The header of the combined product comes from the first position that is kept. Updates
        ``self.qc`` and ``self.products``.

        **Key Arguments:**

        - ``extractions`` -- the ``stack_extraction`` of each nod position, by position, in extraction order

        **Return:**

        - ``combined`` -- the ``stack_extraction`` of the combined spectrum. Its notes carry the position of each note
        """
        kept = {position: extraction for position, extraction in extractions.items() if not extraction.traceFallback}
        notes = self._combined_notes(extractions, kept)
        reference = kept or extractions
        first, last = reference[next(iter(reference))], reference[next(reversed(reference))]
        self.masterHeaderFrame = first.combinedObject
        stackedSpectrum, extractionPath = self.stack_extractions(
            [
                None if position not in kept else extraction.mergedSpectrum
                for position, extraction in extractions.items()
            ],
            orderJoins=last.orderJoins,
        )
        notFlatPath = self._combine_not_flattened(extractions, last.orderJoins)
        return stack_extraction(
            productPath=extractionPath,
            combinedObject=first.combinedObject,
            mergedSpectrum=stackedSpectrum,
            orderJoins=last.orderJoins,
            extractionPath=extractionPath,
            responseInputs=None,
            notFlatPath=notFlatPath,
            notFlatSpectrum=None,
            traceFallback=False,
            slitEdgeSpill=any(extraction.slitEdgeSpill for extraction in reference.values()),
            extractionNotes=notes,
        )

    def _combined_notes(
        self, extractions: dict[str, stack_extraction], kept: dict[str, stack_extraction]
    ) -> list[tuple[str, str]]:
        """*the notes of the combined spectrum's plots: those of each position, tagged, and a warning for a dropped one*

        **Key Arguments:**

        - ``extractions`` -- the ``stack_extraction`` of each nod position, by position
        - ``kept`` -- the positions of ``extractions`` that found an object trace

        **Return:**

        - ``notes`` -- the ``(level, message)`` pairs, in position order
        """
        notes = [
            (level, f"[{position}] {message}")
            for position, extraction in extractions.items()
            for level, message in extraction.extractionNotes
        ]
        if not kept:
            return notes
        for position in extractions:
            if position in kept:
                continue
            other = next(iter(kept))
            warning = (
                f"No object trace was found in the {position} frames, so the stare-mode spectrum is the {other} "
                "position's extraction alone"
            )
            self.log.warning(warning)
            notes.append(("WARNING", warning))
        return notes

    def _detrend_nod_frames(self, frames, master_bias, dark, orderTablePath):
        """*detrend each frame with the master bias or dark and the scattered-light fit, but not the flat*

        **Key Arguments:**

        - ``frames`` -- the frames to detrend
        - ``master_bias`` -- the master bias, or False
        - ``dark`` -- the master dark, or False
        - ``orderTablePath`` -- the path to the order table

        **Return:**

        - ``detrended`` -- the detrended frames, in the order of ``frames``
        """
        return [
            self.detrend(
                inputFrame=frame, master_bias=master_bias, dark=dark, master_flat=False, order_table=orderTablePath
            )
            for frame in frames
        ]

    def _flat_field_nod_frames(self, frames, master_flat):
        """*divide each frame by the master flat*

        **Key Arguments:**

        - ``frames`` -- the frames to flat-field
        - ``master_flat`` -- the master flat

        **Return:**

        - ``flattened`` -- the flat-fielded frames, in the order of ``frames``
        """
        return [
            self.detrend(inputFrame=frame, master_bias=False, dark=False, master_flat=master_flat, order_table=False)
            for frame in frames
        ]

    def _align_nod_frames(
        self,
        framesA: list[Any],
        framesB: list[Any],
        namesA: list[str],
        namesB: list[str],
        twoDMap: str,
        dispMap: str,
    ) -> tuple[list[Any], list[Any], tuple[list[dict[int, float]], list[dict[int, float]]], two_d_map_images | None]:
        """*align the A frames to each other and the B frames to each other along the slit*

        The slit alignment also replaces the items of the ``framesA`` and ``framesB`` lists, to release the original
        frames early. The caller uses the lists this method returns. A position of one frame needs no alignment.

        **Key Arguments:**

        - ``framesA`` -- the A frames
        - ``framesB`` -- the B frames
        - ``namesA`` -- the raw file names of the A frames
        - ``namesB`` -- the raw file names of the B frames
        - ``twoDMap`` -- the path to the dispersion map image
        - ``dispMap`` -- the path to the dispersion map table

        **Return:**

        - ``framesA`` -- the aligned A frames
        - ``framesB`` -- the aligned B frames
        - ``alignmentShifts`` -- the alignment shift of each order of each A frame, and of each B frame, as two lists
          of ``{order: shift in arcsec}`` dictionaries
        - ``mapImages`` -- the ``two_d_map_images`` the alignment read, or *None* when no position needed aligning
        """
        from soxspipe.commonutils.nod_frames import aligns_each_order
        from soxspipe.commonutils.slit_alignment import align_frame_groups_along_slit, frame_group
        from soxspipe.commonutils.slit_coverage import frame_order_shifts

        results = align_frame_groups_along_slit(
            log=self.log,
            groups=[frame_group(framesA, framesB, namesA), frame_group(framesB, framesA, namesB)],
            twoDMapPath=twoDMap,
            dispersionMap=dispMap,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            recipeName=self.recipeName,
            perOrder=aligns_each_order(self.inst, self.arm),
            startNightDate=self.startNightDate,
        )
        mapImages = results[0].mapImages if results[0].mapImages is not None else results[1].mapImages
        alignmentShifts = (frame_order_shifts(results[0], namesA), frame_order_shifts(results[1], namesB))
        return list(results[0].frames), list(results[1].frames), alignmentShifts, mapImages
