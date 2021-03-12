# SOXSPIPE File Naming Convention

With the SOXSPIPE file-naming convention you should be able, at a glance, be able to tell the content of the frame without the need to dig into the file's metadata.

Here is the file-naming scheme for all SOXSPIPE generated frames: 

```text
<DATE>_<OBID>_<ARM>_<BIN>_<TYPE>[_<MASK|SLIT>_<OBSMODE>].fits
```

Those variables within square brackets are present only for the frames that require them.

`<DATE>` 
: UT date-time stamp reflecting the start of the observation to second precision (required)

`<OBID>`
: ID of the observation block used to acquire the frame (required)

`<ARM>`
: the instrument arm the frame originated from. `vis`, `nir` or `acam`  (required)

`<BIN>`
: pixel binning used along both axes

`<TYPE>`
: the type of frame. See [Frame Types](#frame-types) section (required)

`<MASK|SLIT>`
: the type of mask (`onepin` or `mulpin`) or width of slit (in arcsec) used (required for spectroscopic frames)

`<OBSMODE>`
: the observation mode used, `stare`, `nod` or `offset` (required for spectroscopic frames)

## Frame Types

| Type String | Description |
| ----------- | ----------- |
| `mbias`       | master bias |
| `mdark`       | master dark |
| `mskyflat`       | master sky flat |
| `mflat`       | master lamp flat |
| `arc`       | arc lamp frame |
| `flat`       | flat lamp frame |
| `tel_std`       | telluric standard |
| `spec_std`       | spectroscopic standard |
| `phot_std`       | photometric standard |
| `acqu`       | acquisition |
| `obj` | science object |
