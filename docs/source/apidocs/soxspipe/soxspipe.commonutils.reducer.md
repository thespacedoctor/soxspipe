# {py:mod}`soxspipe.commonutils.reducer`

```{py:module} soxspipe.commonutils.reducer
```

```{autodoc2-docstring} soxspipe.commonutils.reducer
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`reducer <soxspipe.commonutils.reducer.reducer>`
  - ```{autodoc2-docstring} soxspipe.commonutils.reducer.reducer
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`run_recipe <soxspipe.commonutils.reducer.run_recipe>`
  - ```{autodoc2-docstring} soxspipe.commonutils.reducer.run_recipe
    :summary:
    ```
* - {py:obj}`run_recipe_bulk <soxspipe.commonutils.reducer.run_recipe_bulk>`
  - ```{autodoc2-docstring} soxspipe.commonutils.reducer.run_recipe_bulk
    :summary:
    ```
````

### API

`````{py:class} reducer(log, workspaceDirectory, reductionTarget='all', settings=False, pathToSettings=False, quitOnFail=False, overwrite=False, daemon=False, verbose=False, refreshWorkspace=False)
:canonical: soxspipe.commonutils.reducer.reducer

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.reducer.reducer
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.reducer.reducer.__init__
```

````{py:method} reduce(batch=False, multiprocess=False)
:canonical: soxspipe.commonutils.reducer.reducer.reduce

```{autodoc2-docstring} soxspipe.commonutils.reducer.reducer.reduce
```

````

````{py:method} select_sof_files_to_process(recipe=False, reductionTarget=False, batch=False, arm=False)
:canonical: soxspipe.commonutils.reducer.reducer.select_sof_files_to_process

```{autodoc2-docstring} soxspipe.commonutils.reducer.reducer.select_sof_files_to_process
```

````

`````

````{py:function} run_recipe(log, recipe, sof, settings, overwrite, command=False, verbose=False, turnOffMP=False)
:canonical: soxspipe.commonutils.reducer.run_recipe

```{autodoc2-docstring} soxspipe.commonutils.reducer.run_recipe
```
````

````{py:function} run_recipe_bulk(log, recipe, sofList, commandList, settings, overwrite, workspaceDirectory, conn, sessionId)
:canonical: soxspipe.commonutils.reducer.run_recipe_bulk

```{autodoc2-docstring} soxspipe.commonutils.reducer.run_recipe_bulk
```
````
