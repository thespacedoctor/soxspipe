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

`````{py:class} data_organiser(log, rootDir, vlt=False, dbConnect=True)
:canonical: soxspipe.commonutils.data_organiser.data_organiser

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.__init__
```

````{py:method} build_sof_files()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.build_sof_files

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.build_sof_files
```

````

````{py:method} close()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.close

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.close
```

````

````{py:method} get_incomplete_raw_frames_set()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.get_incomplete_raw_frames_set

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.get_incomplete_raw_frames_set
```

````

````{py:method} get_raw_frames_and_groups(ttype=None, arm=None, tech=None, recipe=None, recipeOrder=None, filterName=None, unprocessedOnly=False)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.get_raw_frames_and_groups

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.get_raw_frames_and_groups
```

````

````{py:method} list_obs()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.list_obs

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.list_obs
```

````

````{py:method} list_raw(sofFile)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.list_raw

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.list_raw
```

````

````{py:method} list_sofs()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.list_sofs

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.list_sofs
```

````

````{py:method} predict_product_frames(productTypes, rawGroups, recipe)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.predict_product_frames

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.predict_product_frames
```

````

````{py:method} prepare(refresh=False, report=True)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.prepare

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.prepare
```

````

````{py:method} raw_frames_to_sof_map(rawGroups, containerSofs)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.raw_frames_to_sof_map

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.raw_frames_to_sof_map
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

````{py:method} session_refresh(silent=False, failure=True)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.session_refresh

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.session_refresh
```

````

````{py:method} session_switch(sessionId)
:canonical: soxspipe.commonutils.data_organiser.data_organiser.session_switch

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.session_switch
```

````

````{py:method} use_vlt_environment_folders()
:canonical: soxspipe.commonutils.data_organiser.data_organiser.use_vlt_environment_folders

```{autodoc2-docstring} soxspipe.commonutils.data_organiser.data_organiser.use_vlt_environment_folders
```

````

`````
