# {py:mod}`soxspipe.commonutils.data_organiser`

```{py:module} soxspipe.commonutils.data_organiser
```

```{autodoc2-docstring} soxspipe.commonutils.data_organiser
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`data_organiser <soxspipe.commonutils.data_organiser.data_organiser>`
  - ```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser
    :summary:
    ```
````

### API

`````{py:class} data_organiser(log, rootDir)
:canonical: soxspipe.commonutils.data_organiser.data_organiser

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.__init__
```

````{py:method} categorise_frames(filteredFrames, verbose=False)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.categorise_frames

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.categorise_frames
```

````

````{py:method} create_directory_table(pathToDirectory, filterKeys)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.create_directory_table

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.create_directory_table
```

````

````{py:method} generate_sof_and_product_names(series, reductionOrder, rawFrames, calibrationFrames, calibrationTables)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.generate_sof_and_product_names

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.generate_sof_and_product_names
```

````

````{py:method} populate_products_table(series, reductionOrder)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.populate_products_table

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.populate_products_table
```

````

````{py:method} prepare()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.prepare

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.prepare
```

````

````{py:method} session_create(sessionId=False)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.session_create

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.session_create
```

````

````{py:method} session_list(silent=False)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.session_list

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.session_list
```

````

````{py:method} session_refresh()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.session_refresh

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.session_refresh
```

````

````{py:method} session_switch(sessionId)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.session_switch

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.session_switch
```

````

````{py:method} symlink_session_assets_to_workspace_root()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.symlink_session_assets_to_workspace_root

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.symlink_session_assets_to_workspace_root
```

````

````{py:method} sync_sql_table_to_directory(directory, tableName, recursive=False)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.sync_sql_table_to_directory

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.sync_sql_table_to_directory
```

````

`````
