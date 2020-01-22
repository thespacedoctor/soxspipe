"""
*The pipeline recipes*
"""


class baseRecipe():
    """
    The worker class for the baseRecipe module

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary

    **Usage:**

    ```python
    usage code 
    ```

    ---

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - when complete, clean baseRecipe class
    ```
    """
    # Initialisation
    # 1. @flagged: what are the unique attrributes for each object? Add them
    # to __init__

    def __init__(
            self,
            log,
            settings=False,

    ):
        self.log = log
        log.debug("instansiating a new '__init__' object")
        self.settings = settings
        # xt-self-arg-tmpx

        # 2. @flagged: what are the default attrributes each object could have? Add them to variable attribute set here
        # Variable Data Atrributes

        # 3. @flagged: what variable attrributes need overriden in any baseclass(es) used
        # Override Variable Data Atrributes

        # Initial Actions

        return None

    def close(self):
        del self
        return None

    # 4. @flagged: what actions does each object have to be able to perform? Add them here
    # Method Attributes
    def get(self):
        """get the baseRecipe object

        **Return:**
            - ``baseRecipe``

        ..todo::

            - @review: when complete, clean get method
            - @review: when complete add logging
        """
        self.debug.info('starting the ``get`` method')

        baseRecipe = None

        self.debug.info('completed the ``get`` method')
        return baseRecipe
    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
