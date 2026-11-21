"""
*Stack the frames of a nodding SOF that the stare recipe reduces in stare mode (DY-2194)*

Author
: David Young

Date Created
: November 16, 2026
"""

from typing import NamedTuple

from soxspipe.commonutils.toolkit import utcnow_string


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


class nod_stare_stacking:
    """
    Mixin of the methods of `soxs_stare` that align and stack the A and B frames of a nodding SOF before extraction.

    The object sits at one slit position in the A frames and at another in the B frames, and no frame is subtracted
    from another. Each frame is detrended, the frames of each nod position are aligned to each other along the slit and
    the B frames are shifted onto the A frames. The stack then holds the object once, at the A position.
    """

    # THE RAW FILE NAMES OF THE OBJECT FRAMES OF THE RUN, IN THE ORDER THE FRAMES WERE READ. THE SLIT ALIGNMENT NAMES
    # ITS FRAMES WITH THEM
    objectFrameNames = ()

    def _reduce_nod_frames_as_stare(self, allObjectFrames, rawFrames):
        """*align and stack the frames of a nodding SOF, then reduce the stack with the stare tail*

        **Key Arguments:**

        - ``allObjectFrames`` -- the A and B object frames, in observation order
        - ``rawFrames`` -- the raw file names the product header lists, or *False* to list every raw frame of the SOF

        **Return:**

        - ``productPath`` -- the path to the product
        - ``forceFailure`` -- the failure the tail reports
        """
        framesA, framesB, namesA, namesB = self._split_nod_frames(allObjectFrames)
        inputs = self._read_reduction_inputs(self.kw, self.arm)
        notFlattened, flattened, uncovered = self._stack_nod_frames(framesA, framesB, namesA, namesB, inputs)
        return self._reduce_stacked_object(notFlattened, allObjectFrames, rawFrames, inputs, flattened, uncovered)

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

    def _stack_nod_frames(self, framesA, framesB, namesA, namesB, inputs):
        """*detrend, align and stack the A and B frames*

        Each frame is detrended with the master bias or dark and the scattered-light fit, but not flat-fielded. The A
        frames are aligned to each other and the B frames to each other, and the B frames are shifted onto the A
        frames. A frame is flat-fielded before its shift, as the flat is fixed on the detector.

        **Key Arguments:**

        - ``framesA`` -- the A frames
        - ``framesB`` -- the B frames
        - ``namesA`` -- the raw file names of the A frames
        - ``namesB`` -- the raw file names of the B frames
        - ``inputs`` -- the ``reduction_inputs`` of `_read_reduction_inputs`

        **Return:**

        - ``notFlattened`` -- the stack of the detrended frames
        - ``flattened`` -- the stack of the flat-fielded frames. The same frame as ``notFlattened`` when no master flat
          is used
        - ``uncovered`` -- boolean image, True for a pixel of the slit that only the A frames cover (see
          `_uncovered_slit_mask`)

        **Raises:**

        - ``ValueError`` -- when the shift of the B frames leaves no slit range that every frame covers in any order
        """
        framesA = self._detrend_nod_frames(framesA, inputs.masterBias, inputs.dark, inputs.orderTablePath)
        framesB = self._detrend_nod_frames(framesB, inputs.masterBias, inputs.dark, inputs.orderTablePath)
        if len(framesA) + len(framesB) > 1:
            # EACH FRAME'S SCATTERED-LIGHT FIT ADDED A ROW FOR THE ONE PLOT FILE THE LAST FRAME LEFT
            self.products = self.products.drop_duplicates(subset=["product_label", "file_name"], keep="last")

        framesA, framesB, plan = self._align_nod_frames(
            framesA, framesB, namesA, namesB, inputs.twoDMap, inputs.dispMap
        )
        uncovered = self._uncovered_slit_mask(plan)
        notFlattened = self._stack_shifted_nod_frames(framesA, framesB, plan)
        if inputs.masterFlat is False:
            return notFlattened, notFlattened, uncovered

        framesA = self._flat_field_nod_frames(framesA, inputs.masterFlat)
        framesB = self._flat_field_nod_frames(framesB, inputs.masterFlat)
        return notFlattened, self._stack_shifted_nod_frames(framesA, framesB, plan), uncovered

    def _uncovered_slit_mask(self, plan):
        """*the pixels of the slit that the shifted B frames do not reach, which only the A frames cover (DY-2290)*

        The B frames are shifted onto the A frames, so the slit beyond the shift's reach holds A frames only. The sky
        changes between the A and the B exposures, so the sky residual of the stack steps where the B frames end, and
        the step reads as a second source. The stack keeps only the slit range every frame covers: the caller leaves
        the returned pixels out of the sky fit and the extraction and flags them in the products. Orders with no shift
        have none.

        **Key Arguments:**

        - ``plan`` -- the ``shift_and_add_plan`` whose ``orderShifts`` move a B frame onto the A frames

        **Return:**

        - ``uncovered`` -- boolean image, True for a pixel a shifted B frame does not cover

        **Raises:**

        - ``ValueError`` -- when the shift leaves no slit range that every frame covers in any shifted order
        """
        from soxspipe.commonutils.slit_alignment import shift_coverage_mask

        uncovered = shift_coverage_mask(plan.orderShifts, plan.mapImages)
        orderImage, orders = plan.mapImages.order, plan.mapImages.orders
        emptyOrders = [order for order in orders if uncovered[orderImage == order].all()]
        if emptyOrders and len(emptyOrders) == len(orders):
            raise ValueError(
                "The shifted B frames leave no slit range that every frame covers in any order, so the frames cannot "
                "be stacked in stare mode."
            )
        if emptyOrders:
            orderList = ", ".join(str(order) for order in emptyOrders)
            self.log.warning(
                f"The shifted B frames leave no slit range that every frame covers in order {orderList}, so it is "
                "left out of the sky model and the extraction."
            )
        self.log.info(
            f"{int(uncovered.sum())} pixels that only the A frames cover are left out of the sky model and the "
            "extraction."
        )
        return uncovered

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

    def _align_nod_frames(self, framesA, framesB, namesA, namesB, twoDMap, dispMap):
        """*align the A frames to each other and the B frames to each other, and plan the shift of B onto A*

        The slit alignment also replaces the items of the ``framesA`` and ``framesB`` lists, to release the original
        frames early. The caller uses the lists this method returns. The plan comes from the object position measured
        at both nod positions. When it cannot be measured there, the plan is the nominal nod throw of the headers.
        Updates ``self.qc``.

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
        - ``plan`` -- the ``shift_and_add_plan`` whose ``orderShifts`` move a B frame onto the A frames

        **Raises:**

        - ``ValueError`` -- when the alignment measured no group, so there is no 2D map for the nominal throw
        """
        from soxspipe.commonutils.nod_frames import aligns_each_order, nominal_throw_shift_plan
        from soxspipe.commonutils.nod_shift_and_add import build_shift_and_add_plan
        from soxspipe.commonutils.slit_alignment import align_frame_groups_along_slit, frame_group

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
            measureSingleFrames=True,
        )
        framesA, framesB = list(results[0].frames), list(results[1].frames)
        plan = build_shift_and_add_plan(self.log, *results)
        isNominal = plan is None
        if isNominal:
            mapImages = results[0].mapImages if results[0].mapImages is not None else results[1].mapImages
            if mapImages is None:
                raise ValueError("The slit alignment read no 2D map, so the nominal nod throw cannot be applied.")
            plan = nominal_throw_shift_plan(self.log, framesA, framesB, mapImages, self.kw)
        self._record_nod_shift_source(isNominal, framesA[0].header[self.kw("DATE_OBS")])
        return framesA, framesB, plan

    def _record_nod_shift_source(self, isNominal, obsDateUtc):
        """*record in the QC table whether the B frames were shifted by the nominal nod throw*

        The QC has no acceptable range, so it never fails the recipe. Updates ``self.qc``.

        **Key Arguments:**

        - ``isNominal`` -- True when the nominal nod throw of the headers was used, False when the object position was
          measured at both nod positions
        - ``obsDateUtc`` -- the observation date to record against the QC
        """
        self.add_qc(
            qcName="NOD STARE NOMINAL THROW",
            qcValue=int(isNominal),
            qcComment="1 if the B frames were shifted by the commanded nod throw, as no object position was measured",
            reductionDateUtc=utcnow_string(),
            qcUnit=None,
            toHeader=True,
            obsDateUtc=obsDateUtc,
        )

    def _stack_shifted_nod_frames(self, framesA, framesB, plan):
        """*shift the B frames onto the A frames and stack all of them, leaving the shift's masked pixels out*

        **Key Arguments:**

        - ``framesA`` -- the aligned A frames
        - ``framesB`` -- the aligned B frames
        - ``plan`` -- the ``shift_and_add_plan`` whose ``orderShifts`` move a B frame onto the A frames

        **Return:**

        - ``stacked`` -- the stack. A pixel the shift flagged as off the slit in a B frame is left out of that frame
        """
        from soxspipe.commonutils.slit_alignment import shift_frame_along_slit

        mjdKey = self.kw("MJDOBS")
        shiftedB = [shift_frame_along_slit(frame, plan.orderShifts, plan.mapImages) for frame in framesB]
        frames = sorted([*framesA, *shiftedB], key=lambda frame: frame.header[mjdKey])
        return self.clip_and_stack(
            frames=frames,
            recipe="soxs_stare",
            ignore_input_masks=False,
            post_stack_clipping=False,
        )
