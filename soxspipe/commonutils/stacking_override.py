"""*the settings override behind the `--no-stack` command-line option (DY-1723)*

The option turns off `stack_before_extraction` for one run. It also marks the override as asked for on the command
line, so that standard-star SOFs can be reduced frame by frame (or cycle by cycle), which a settings-file value alone
never allows.
"""

import copy

STACKING_KEY = "stack_before_extraction"
OVERRIDE_KEY = "no-stack"
STACKING_RECIPE_BLOCKS = (
    "soxs-stare",
    "soxs-stare-std",
    "soxs-nod",
    "soxs-nod-std",
    "soxs-offset",
    "soxs-offset-std",
)


def settings_without_stacking(settings):
    """*return a copy of the settings with stacking before extraction turned off for stare, nod and offset*

    **Key Arguments:**

    - ``settings`` -- the settings dictionary. It is not changed

    **Return:**

    - ``settings`` -- a deep copy with `stack_before_extraction` set to False in every stare, nod and offset recipe
      block that is present, and in every arm sub-block of those that carries the key. The top-level key `no-stack`
      is True

    **Usage:**

    ```python
    from soxspipe.commonutils.stacking_override import settings_without_stacking
    settings = settings_without_stacking(settings)
    ```
    """
    overridden = copy.deepcopy(settings)
    for blockName in STACKING_RECIPE_BLOCKS:
        block = overridden.get(blockName)
        if not isinstance(block, dict):
            continue
        block[STACKING_KEY] = False
        for armBlock in block.values():
            if isinstance(armBlock, dict) and STACKING_KEY in armBlock:
                armBlock[STACKING_KEY] = False
    overridden[OVERRIDE_KEY] = True
    return overridden
