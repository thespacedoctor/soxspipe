# Order Table

The order table gets built over [`soxs_order_centres`](../recipes/soxs_order_centres.md) and [`soxs_mflat`](../recipes/soxs_mflat.md) recipes. The final product contains polynomial fits for the order centres and the upper and lower edges of the order locations on the detector.

Here is an example of the content of the NIR order table from the [`soxs_order_centres`](../recipes/soxs_order_centres.md) recipes.

```text
order,degy,CENT_c0,CENT_c1,CENT_c2,CENT_c3
11,3,787.9616,0.387875383,-0.000175467653,7.67345562e-10
12,3,714.084224,0.392816612,-0.000177468713,1.24058095e-09
...
```

## Extensions

The order table is a FITS file with three binary-table extensions.

1. HDU 1 holds the polynomial coefficients, as shown above.
2. HDU 2 holds the order metadata: the order number and the fitted range of each order along both axes.
3. HDU 3, named `ORDER_TRACES`, lists the pixel coordinates of each order trace. Read it by name, not by number.

`ORDER_TRACES` is one long table with one row for each pixel along the dispersion axis, inside the fitted range of each order. The sampled range runs from the lower limit of the order up to, but not including, its upper limit. This is the same range that `unpack_order_table` returns when it samples every pixel. The coordinates are in unbinned detector pixels, whatever the binning of the frame that the order table was measured on. The columns are named after the detector axes. For the SOXS VIS arm and every X-Shooter arm, where the traces are solved as `x` for each `y`, the columns are:

| Column | Content |
| --- | --- |
| `ycoord` | the pixel position along the axis the traces are sampled on |
| `order` | the order number |
| `xcoord_centre` | the position of the order centre |
| `xcoord_edgeup` | the position of the upper order edge |
| `xcoord_edgelow` | the position of the lower order edge |

For the SOXS NIR arm, the traces are solved as `y` for each `x`, so the axes swap: the columns are `xcoord`, `order`, `ycoord_centre`, `ycoord_edgeup` and `ycoord_edgelow`.

The edge columns are only in the order table that [`soxs_mflat`](../recipes/soxs_mflat.md) writes. The table from [`soxs_order_centres`](../recipes/soxs_order_centres.md) has the centre column only. Order tables written before this extension existed have two extensions. The object trace tables (`OBJECT_TAB_*`) do not have an `ORDER_TRACES` extension.
