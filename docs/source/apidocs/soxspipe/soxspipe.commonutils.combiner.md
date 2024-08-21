# {py:mod}`soxspipe.commonutils.combiner`

```{py:module} soxspipe.commonutils.combiner
```

```{autodoc2-docstring} soxspipe.commonutils.combiner
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`Combiner <soxspipe.commonutils.combiner.Combiner>`
  -
````

### API

`````{py:class} Combiner(ccd_iter, dtype=None)
:canonical: soxspipe.commonutils.combiner.Combiner

Bases: {py:obj}`ccdproc.Combiner`

````{py:method} average_combine()
:canonical: soxspipe.commonutils.combiner.Combiner.average_combine

```{autodoc2-docstring} soxspipe.commonutils.combiner.Combiner.average_combine
```

````

````{py:method} clip_extrema(nlow=0, nhigh=0)
:canonical: soxspipe.commonutils.combiner.Combiner.clip_extrema

````

````{py:property} dtype
:canonical: soxspipe.commonutils.combiner.Combiner.dtype

````

````{py:method} median_combine(median_func=None, scale_to=None, uncertainty_func=sigma_func)
:canonical: soxspipe.commonutils.combiner.Combiner.median_combine

````

````{py:method} minmax_clipping(min_clip=None, max_clip=None)
:canonical: soxspipe.commonutils.combiner.Combiner.minmax_clipping

````

````{py:property} scaling
:canonical: soxspipe.commonutils.combiner.Combiner.scaling

````

````{py:method} sigma_clipping(low_thresh=3, high_thresh=3, func='mean', dev_func='std', **kwd)
:canonical: soxspipe.commonutils.combiner.Combiner.sigma_clipping

````

````{py:method} sum_combine(sum_func=None, scale_to=None, uncertainty_func=_default_std())
:canonical: soxspipe.commonutils.combiner.Combiner.sum_combine

````

````{py:property} weights
:canonical: soxspipe.commonutils.combiner.Combiner.weights

````

`````
