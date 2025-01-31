# Order Table

The order table gets built over [`soxs_order_centres`](../recipes/soxs_order_centres.md) and [`soxs_mflat`](../recipes/soxs_mflat.md) recipes. The final product contains polynomial fits for the order centres and the upper and lower edges of the order locations on the detector.

Here is an example of the content of the NIR order table from the [`soxs_order_centres`](../recipes/soxs_order_centres.md) recipes.

```text
order,degy,CENT_c0,CENT_c1,CENT_c2,CENT_c3
11,3,787.9616,0.387875383,-0.000175467653,7.67345562e-10
12,3,714.084224,0.392816612,-0.000177468713,1.24058095e-09
...
```

---

:::{todo}

- add order table example after soxs_mflat

:::
