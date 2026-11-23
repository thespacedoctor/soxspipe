"""*the hand-off of a nodding SOF to the stare recipe (DY-2194)*

The ``--stare`` command-line option and the ``reduce-as-stare`` settings key reduce the frames of a nodding SOF with the
stare algorithm, to avoid the A-B subtraction removing extended host emission. The nod and the reducer entry points ask
`nod_recipe` for a recipe object, and it returns either the nod recipe or the stare recipe in nod mode.
"""

REDUCE_AS_STARE_KEY = "reduce-as-stare"

# THE MARKER IN THE NAME OF EVERY PRODUCT OF A NODDING SOF REDUCED IN STARE MODE, SO NO PRODUCT OVERWRITES A NOD PRODUCT
NOD_STARE_MARKER = "_STAREMODE"


def is_standard_input(inputFrames):
    """*check whether the input frames are a standard-star SOF, by the ``_STD_`` rule the recipes use for their name*

    **Key Arguments:**

    - ``inputFrames`` -- the recipe's input: a SOF path, a directory or a list of frame paths

    **Return:**

    - ``isStandard`` -- True when ``inputFrames`` is a string with ``_STD_`` in it
    """
    return isinstance(inputFrames, str) and "_STD_" in inputFrames


def reduces_nod_as_stare(settings, isStandard, isRequested=False):
    """*decide whether a nodding SOF is reduced in stare mode*

    The request (the ``--stare`` option) wins. Without it the ``reduce-as-stare`` key of the ``soxs-nod-std`` block (a
    standard) or the ``soxs-nod`` block decides, and a missing block or key means nod mode.

    **Key Arguments:**

    - ``settings`` -- the settings dictionary
    - ``isStandard`` -- True for a standard-star SOF
    - ``isRequested`` -- True when the command line asked for stare mode. Default *False*

    **Return:**

    - ``isStare`` -- True to reduce the SOF in stare mode
    """
    if isRequested:
        return True
    block = settings.get("soxs-nod-std" if isStandard else "soxs-nod")
    if not isinstance(block, dict):
        return False
    return bool(block.get(REDUCE_AS_STARE_KEY, False))


def nod_recipe(log, settings, inputFrames, isRequested=False, **recipeKwargs):
    """*build the recipe that reduces a nodding SOF: the nod recipe, or the stare recipe in nod mode*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- the SOF the recipe reduces
    - ``isRequested`` -- True when the command line asked for stare mode. Default *False*
    - ``recipeKwargs`` -- the other arguments of the recipe constructors (``verbose``, ``overwrite``, ``command``, ...)

    **Return:**

    - ``recipe`` -- a ``soxs_stare`` built with ``reduceNodAsStare=True``, or a ``soxs_nod``

    **Usage:**

    ```python
    from soxspipe.commonutils.stare_mode import nod_recipe
    productPath, qcTable = nod_recipe(log=log, settings=settings, inputFrames=sofPath).produce_product()
    ```
    """
    if reduces_nod_as_stare(settings, is_standard_input(inputFrames), isRequested):
        from soxspipe.recipes import soxs_stare

        return soxs_stare(log=log, settings=settings, inputFrames=inputFrames, reduceNodAsStare=True, **recipeKwargs)

    from soxspipe.recipes import soxs_nod

    return soxs_nod(log=log, settings=settings, inputFrames=inputFrames, **recipeKwargs)


def nod_recipe_name_and_marker(settings, isStandard, isRequested=False):
    """*name the recipe that writes the products of a nodding SOF, and the marker its product file names carry*

    The nod mode and the stare mode write their products to different directories, and the stare mode names them with
    ``NOD_STARE_MARKER``.

    **Key Arguments:**

    - ``settings`` -- the settings dictionary
    - ``isStandard`` -- True for a standard-star SOF
    - ``isRequested`` -- True when the command line asked for stare mode. Default *False*

    **Return:**

    - ``recipeName`` -- ``soxs-nod`` or ``soxs-stare``, with ``-std`` added for a standard
    - ``marker`` -- ``NOD_STARE_MARKER`` in stare mode, else an empty string

    **Usage:**

    ```python
    from soxspipe.commonutils.stare_mode import nod_recipe_name_and_marker
    recipeName, marker = nod_recipe_name_and_marker(settings, isStandard=False)
    ```
    """
    isStare = reduces_nod_as_stare(settings, isStandard, isRequested)
    recipeName = ("soxs-stare" if isStare else "soxs-nod") + ("-std" if isStandard else "")
    return recipeName, NOD_STARE_MARKER if isStare else ""


def predicted_nod_product_path(settings, inputFrames, isRequested=False):
    """*predict the product path of a nodding SOF in the mode the recipe will run in*

    The nod mode and the stare mode write their products to different directories, and the stare mode names them with
    ``NOD_STARE_MARKER``.

    **Key Arguments:**

    - ``settings`` -- the settings dictionary. Its ``workspace-root-dir``, which ``-o`` sets, is where the recipe writes
      the product
    - ``inputFrames`` -- the SOF the recipe reduces
    - ``isRequested`` -- True when the command line asked for stare mode. Default *False*

    **Return:**

    - ``productPath`` -- the predicted path of the product of the recipe
    """
    from soxspipe.commonutils import toolkit

    recipeName, marker = nod_recipe_name_and_marker(settings, is_standard_input(inputFrames), isRequested)
    productPath, _ = toolkit.predict_product_path(
        inputFrames, recipeName, marker, outputDir=settings.get("workspace-root-dir")
    )
    return productPath
