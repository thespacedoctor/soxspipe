# {py:mod}`soxspipe.commonutils.set_of_files`

```{py:module} soxspipe.commonutils.set_of_files
```

```{autodoc2-docstring} soxspipe.commonutils.set_of_files
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`ImageFileCollection <soxspipe.commonutils.set_of_files.ImageFileCollection>`
  -
* - {py:obj}`set_of_files <soxspipe.commonutils.set_of_files.set_of_files>`
  - ```{autodoc2-docstring} soxspipe.commonutils.set_of_files.set_of_files
    :summary:
    ```
````

### API

`````{py:class} ImageFileCollection(location=None, keywords=None, find_fits_by_reading=False, filenames=None, glob_include=None, glob_exclude=None, ext=0)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection

Bases: {py:obj}`soxspipe.commonutils.set_of_files.ImageFileCollection`

````{py:method} __repr__()
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.__repr__

````

````{py:method} ccds(ccd_kwargs=None, **kwd)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.ccds

````

````{py:method} data(do_not_scale_image_data=False, **kwd)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.data

````

````{py:property} ext
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.ext

````

````{py:property} files
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.files

````

````{py:method} files_filtered(**kwd)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.files_filtered

````

````{py:method} filter(**kwd)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.filter

````

````{py:property} glob_exclude
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.glob_exclude

````

````{py:property} glob_include
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.glob_include

````

````{py:method} hdus(do_not_scale_image_data=False, **kwd)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.hdus

````

````{py:method} headers(do_not_scale_image_data=True, **kwd)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.headers

````

````{py:property} keywords
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.keywords

````

````{py:property} location
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.location

````

````{py:method} refresh()
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.refresh

````

````{py:method} sort(keys)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.sort

````

````{py:property} summary
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.summary

````

````{py:method} values(keyword, unique=False)
:canonical: soxspipe.commonutils.set_of_files.ImageFileCollection.values

````

`````

`````{py:class} set_of_files(log, settings=False, inputFrames=[], verbose=True, recipeName=False, ext=0, session=None)
:canonical: soxspipe.commonutils.set_of_files.set_of_files

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.set_of_files.set_of_files
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.set_of_files.set_of_files.__init__
```

````{py:method} create_supplementary_file_dictionary(supplementaryFilepaths)
:canonical: soxspipe.commonutils.set_of_files.set_of_files.create_supplementary_file_dictionary

```{autodoc2-docstring} soxspipe.commonutils.set_of_files.set_of_files.create_supplementary_file_dictionary
```

````

````{py:method} get()
:canonical: soxspipe.commonutils.set_of_files.set_of_files.get

```{autodoc2-docstring} soxspipe.commonutils.set_of_files.set_of_files.get
```

````

`````
