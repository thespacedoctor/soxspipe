"""Strict-xfail pins for the resource leaks found by the DY-36 defect hunt.

Every test here pins a leak that has been reported but not fixed. Each test asserts the
CORRECT behaviour, so it fails today and is marked ``xfail(strict=True)``. When someone fixes
the leak the test XPASSes, strict mode turns that into a suite failure, and the fixer must
delete the ``xfail`` marker.

The ``reason`` of each marker reads ``<Linear issue> (DY-36 <ID>): <defect>``, naming the
issue that tracks the leak and the leak's ID in the DY-36 findings table.

Each marker names the exception the leak produces today (``AssertionError``), so a setup,
import or fixture error fails the suite instead of hiding as an xfail. Each pin also asserts
that the production call it measures really ran, so deleting the call cannot make it pass.

Several builders are imported from neighbouring test modules. That coupling is accepted for
now: if one of those helpers is renamed or moved, retarget the import here.
"""

from __future__ import annotations

import contextlib
import gc
import sys
import warnings
import weakref
from collections.abc import Iterator
from importlib import import_module
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import psutil
import pytest
from astropy.io import fits
from astropy.table import Table

import soxspipe.commonutils as commonutils
from soxspipe.commonutils import toolkit
from soxspipe.recipes import soxs_mbias
from tests.factories import product_table, qc_table
from tests.integration.test_recipe_from_frame_list import _bias_frames, _IsolatedOrganiser, _settings
from tests.integration.test_toolkit_split_characterization import (
    KEYWORD_MAP,
    _two_by_two_binned_2d_map,
    _vis_frame,
)
from tests.unit.test_cli_dispatch import RecordingOrganiser, _run_cli

pytestmark = pytest.mark.unit

# HOW LONG A POOL WORKER GETS TO EXIT ON ITS OWN BEFORE THE PIN COUNTS IT AS LEAKED
WORKER_EXIT_GRACE_SECONDS = 3.0
# THE ONE INPUT ON WHICH THE L3 POOL TASK RAISES, SO THE PARALLEL MAP FAILS PART-WAY THROUGH
FAILING_TASK_INPUT = 3


@pytest.fixture
def cyclic_gc_disabled() -> Iterator[None]:
    """Turn the cyclic garbage collector off, so only reference counting frees objects."""
    gc.collect()
    wasEnabled = gc.isenabled()
    gc.disable()
    try:
        yield
    finally:
        if wasEnabled:
            gc.enable()
        gc.collect()


def _reducer_module() -> Any:
    """Return the reducer module; the package re-exports a class of the same name over it."""
    import_module("soxspipe.commonutils.reducer")
    return sys.modules["soxspipe.commonutils.reducer"]


def _resource_warnings_for(path: Path, caught: list[warnings.WarningMessage]) -> list[str]:
    """Return the messages of the caught ResourceWarnings that name ``path``."""
    return [str(item.message) for item in caught if item.category is ResourceWarning and str(path) in str(item.message)]


# ---------------------------------------------------------------------------
# L1: A FAILED RECIPE LEAVES ITS SCRATCH DIRECTORY BEHIND
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "DY-1700 (DY-36 L1): a recipe that raises skips clean_up, so its workspace tmp/ scratch directory and the "
        "prepared frames in it are never removed on the serial `reduce` path; related to DY-37 I4 (kill path, "
        "tests/unit/test_state_and_time_defect_pins.py); file together, fix must clear both pins"
    ),
)
def test_a_recipe_that_raises_removes_its_scratch_directory(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A real master-bias run that fails part-way must not leave its scratch directory behind.

    Related to DY-37 I4 (kill path, tests/unit/test_state_and_time_defect_pins.py); file together, fix must
    clear both pins.
    """
    # ARRANGE
    monkeypatch.setattr(commonutils, "data_organiser", _IsolatedOrganiser)
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    framePaths = _bias_frames(tmp_path / "raw")
    seenAtFailure: dict[str, list[str]] = {}

    def failing_qc(self: Any, *args: object, **kwargs: object) -> None:
        # THE FAILURE LANDS AFTER THE FRAMES ARE PREPARED AND COMBINED, AS A REAL MID-RECIPE ERROR WOULD
        seenAtFailure[self.outDir] = sorted(path.name for path in Path(self.outDir).iterdir())
        raise RuntimeError("synthetic mid-recipe failure")

    monkeypatch.setattr(soxs_mbias, "qc_bias_structure", failing_qc)

    # ACT
    with pytest.raises(RuntimeError, match="synthetic mid-recipe failure"):
        _reducer_module().run_recipe(
            log,
            "mbias",
            framePaths,
            settings=_settings(workspace),
            overwrite=False,
            turnOffMP=True,
        )

    # ASSERT: THE PRODUCTION RUN REALLY PREPARED FRAMES INTO ITS OWN SCRATCH DIRECTORY
    assert len(seenAtFailure) == 1
    (scratchDirectory, preparedFrames) = next(iter(seenAtFailure.items()))
    assert workspace in Path(scratchDirectory).parents
    assert any(name.endswith("_pre.fits") for name in preparedFrames)

    # ASSERT: AND REMOVED IT DESPITE THE FAILURE
    assert not Path(scratchDirectory).exists()


# ---------------------------------------------------------------------------
# L2: A FINISHED RECIPE IS FREED ONLY BY THE CYCLIC GARBAGE COLLECTOR
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "DY-1707 (DY-36 L2): after run_recipe returns, the finished recipe and every array it holds stay alive "
        "until a cyclic gc pass, because the recipe is reachable from frames in reference cycles"
    ),
)
def test_a_finished_recipe_is_freed_when_run_recipe_returns(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    cyclic_gc_disabled: None,
) -> None:
    """`del soxs_recipe` in `run_recipe` must release the recipe without waiting for the cyclic collector."""
    # ARRANGE
    monkeypatch.setattr(commonutils, "data_organiser", _IsolatedOrganiser)
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    framePaths = _bias_frames(tmp_path / "raw")
    recipeRefs: list[weakref.ReferenceType[Any]] = []
    originalProduce = soxs_mbias.produce_product

    def tracked_produce(self: Any) -> Any:
        recipeRefs.append(weakref.ref(self))
        return originalProduce(self)

    monkeypatch.setattr(soxs_mbias, "produce_product", tracked_produce)

    # ACT
    productPath, _ = _reducer_module().run_recipe(
        log,
        "mbias",
        framePaths,
        settings=_settings(workspace),
        overwrite=False,
        turnOffMP=True,
    )

    # ASSERT: THE PRODUCTION RECIPE RAN AND WROTE ITS PRODUCT
    assert len(recipeRefs) == 1
    assert Path(productPath).exists()

    # ASSERT: AND WAS FREED BY REFERENCE COUNTING ALONE
    assert recipeRefs[0]() is None


# ---------------------------------------------------------------------------
# L3: fundamentals.fmultiprocess LEAVES ITS POOL WORKERS RUNNING WHEN A TASK RAISES
# ---------------------------------------------------------------------------


def _task_that_fails_once(value: int, log: object) -> int:
    if value == FAILING_TASK_INPUT:
        raise ValueError("synthetic task failure")
    return value


def _live_children() -> list[psutil.Process]:
    return [child for child in psutil.Process().children(recursive=True) if child.is_running()]


def _is_alive_and_not_zombie(process: psutil.Process) -> bool:
    """Return True when ``process`` still runs; a terminated but unreaped zombie counts as stopped."""
    try:
        return process.is_running() and process.status() != psutil.STATUS_ZOMBIE
    except psutil.NoSuchProcess:
        return False


def _stop_children_not_in(childrenBefore: set[int]) -> None:
    """Stop every child process this test started, letting the abandoned pool shut itself down first.

    Signalling pool workers directly can kill one while it holds the task-queue lock, which deadlocks
    the pool's own finalizer, and the pool's maintenance thread replaces dead workers. So a cyclic gc
    pass runs the pool's finalizer first, and only the processes that survive it are terminated.
    """
    gc.collect()
    started = [child for child in _live_children() if child.pid not in childrenBefore]
    for child in started:
        with contextlib.suppress(psutil.NoSuchProcess):
            child.terminate()
    _, stillAlive = psutil.wait_procs(started, timeout=WORKER_EXIT_GRACE_SECONDS)
    for child in stillAlive:
        with contextlib.suppress(psutil.NoSuchProcess):
            child.kill()
    psutil.wait_procs(stillAlive, timeout=WORKER_EXIT_GRACE_SECONDS)


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "DY-1701 (DY-36 L3): fundamentals.fmultiprocess skips Pool.close/join when a task raises, so the pool's "
        "worker processes stay alive until a cyclic gc pass (an upstream defect in fundamentals)"
    ),
)
@pytest.mark.filterwarnings("ignore:.*fork.*:DeprecationWarning")
def test_fmultiprocess_stops_its_workers_when_a_task_raises(log: Any, cyclic_gc_disabled: None) -> None:
    """A failed parallel map must not leave its worker processes running."""
    from fundamentals import fmultiprocess

    # ARRANGE
    childrenBefore = {child.pid for child in _live_children()}

    try:
        # ACT
        with pytest.raises(ValueError, match="synthetic task failure"):
            fmultiprocess(
                log=log,
                function=_task_that_fails_once,
                inputArray=list(range(8)),
                poolSize=2,
                turnOffMP=False,
                progressBar=False,
            )
        candidates = [child for child in _live_children() if child.pid not in childrenBefore]
        _, stillAlive = psutil.wait_procs(candidates, timeout=WORKER_EXIT_GRACE_SECONDS)
        leaked = [child for child in stillAlive if _is_alive_and_not_zombie(child)]

        # ASSERT
        assert leaked == []
    finally:
        # NEVER LET A WORKER OUTLIVE THE TEST, EVEN WHEN THE ASSERTION FAILS
        _stop_children_not_in(childrenBefore)


# ---------------------------------------------------------------------------
# L4: `soxspipe reduce` KEEPS THE SESSION-LOOKUP DATABASE CONNECTION OPEN FOR THE WHOLE RUN
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "DY-1708 (DY-36 L4): cl_utils.main opens a data_organiser (and its soxspipe.db connection) to look up the "
        "session and never closes it, so it stays open for the whole reduction"
    ),
)
def test_cli_reduce_closes_the_session_lookup_organiser_before_reducing(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Every organiser `main` opens before `reduce` starts must already be closed."""
    # ARRANGE
    (tmp_path / "soxspipe.db").touch()
    organiserCallsAtReduce: list[list[tuple[str, object]]] = []

    class SnapshotReducer:
        def __init__(self, **kwargs: object) -> None:
            return None

        def reduce(self, **kwargs: object) -> None:
            organiserCallsAtReduce.append(list(RecordingOrganiser.calls))

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "reducer", SnapshotReducer)

    # ACT
    _run_cli(monkeypatch, tmp_path, ["soxspipe", "reduce", "all"])

    # ASSERT: THE SESSION LOOKUP REALLY OPENED AN ORGANISER BEFORE THE REDUCTION
    assert len(organiserCallsAtReduce) == 1
    calls = [name for name, _ in organiserCallsAtReduce[0]]
    assert ("session_list", {"silent": True}) in organiserCallsAtReduce[0]

    # ASSERT: AND CLOSED IT BEFORE THE REDUCTION STARTED
    assert calls.count("close") == calls.count("init")


# ---------------------------------------------------------------------------
# L5: TWO READERS MUST CLOSE THE FITS FILE THEY OPEN (FIXED IN DY-1709; THESE ARE REGRESSION TESTS)
# ---------------------------------------------------------------------------


def test_disp_map_dataframe_reader_closes_the_map_file(log: Any, tmp_path: Path) -> None:
    """Reading a 2D dispersion map must close the map file before returning."""
    # ARRANGE
    mapPath = _two_by_two_binned_2d_map(tmp_path / "map.fits")
    frame = _vis_frame((2, 2), binx=4, biny=4)

    # ACT
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        result, _ = toolkit.twoD_disp_map_image_to_dataframe(
            log,
            slit_length=4,
            twoDMapPath=str(mapPath),
            kw=lambda key: KEYWORD_MAP[key],
            associatedFrame=frame,
            dispAxis="y",
        )
        gc.collect()

    # ASSERT: THE PRODUCTION READER RAN ON THE MAP
    assert len(result) == 4

    # ASSERT: WITHOUT LEAVING THE FILE FOR THE GARBAGE COLLECTOR TO CLOSE
    assert _resource_warnings_for(mapPath, caught) == []


def _standard_star_extraction(tmp_path: Path) -> tuple[Path, Path]:
    """Write a two-sample extracted standard and its catalogue table, as the constructor test does."""
    extractedPath = tmp_path / "standard.fits"
    standardsPath = tmp_path / "standards.fits"
    header = fits.Header(
        {
            "SEQ_ARM": "VIS",
            "WIN_BINY": 1,
            "DATE_OBS": "2024-01-02T03:04:05",
            "EXPTIME": 30.0,
            "OBJECT": "CD-30 17706",
            "HIERARCH ESO TEL AIRM START": 1.1,
            "HIERARCH ESO TEL AIRM END": 1.3,
        }
    )
    extractedTable = Table({"WAVE": [600.0, 500.0], "FLUX_COUNTS": [20.0, 10.0]})
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU(extractedTable)]).writeto(extractedPath)
    Table({"WAVE": [500.0, 600.0], "LTT7987": [1.0, 2.0]}).write(standardsPath)
    return extractedPath, standardsPath


def test_response_function_closes_the_standard_extraction_file(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Building a response function must close the extraction file it read the header from."""
    from soxspipe.commonutils.response_function import response_function

    # ARRANGE
    extractedPath, _ = _standard_star_extraction(tmp_path)
    monkeypatch.setattr(commonutils, "keyword_lookup", lambda **_: SimpleNamespace(get=lambda key: key))
    monkeypatch.setattr(
        commonutils,
        "detector_lookup",
        lambda **_: SimpleNamespace(get=lambda arm: {"flux-standards": "standards.fits"}),
    )
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: str(tmp_path))
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")))

    # ACT
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        response = response_function(
            log=log,
            stdExtractionPath=str(extractedPath),
            recipeName="soxs-response",
            settings={"instrument": "soxs", "soxs-response": {"vis": {}}},
            qcTable=qc_table().iloc[0:0],
            productsTable=product_table().iloc[0:0],
            sofName="SYNTHETIC",
            startNightDate="2024-01-02",
            orderJoins=["12"],
            airmass=1.2,
        )
        gc.collect()

    # ASSERT: THE CONSTRUCTOR REALLY READ THE EXTRACTION HEADER
    assert response.arm == "VIS"
    assert np.isclose(response.texp, 30.0)

    # ASSERT: WITHOUT LEAVING THE FILE FOR THE GARBAGE COLLECTOR TO CLOSE
    assert _resource_warnings_for(extractedPath, caught) == []
