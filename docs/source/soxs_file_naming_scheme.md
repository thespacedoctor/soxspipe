# SOXSPIPE File Naming Convention

With the SOXSPIPE file-naming convention you should be able, at a glance, be able to tell the content of the frame without the need to dig into the file's metadata.

Here is the file-naming scheme for all SOXSPIPE generated frames: 

```text
<DATE>_<ARM>_<BIN>_[<ROMODE>_]<TYPE>[_<MASK|SLIT>_<OBSMODE>].fits
```

Those variables within square brackets are present only for the frames that require them.

`<DATE>` 
: UT date-time stamp reflecting the start of the observation to second precision (required)

<!--
`<OBID>`
: ID of the observation block used to acquire the frame (required)
-->

`<ARM>`
: the instrument arm the frame originated from. `vis`, `nir` or `acam`  (required)

`<BIN>`
: pixel binning used along both axes

`<TYPE>`
: the type of frame. See [Frame Types](#frame-types) section (required)

`<ROMODE>`
: readout mode (required)

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

## Stare-mode products of nodding frames

You can reduce the frames of a nodding SOF in stare mode (the `--stare` option, or the `reduce-as-stare` setting). The recipe then writes the product files of the SOF with the marker `_STAREMODE` after the SOF name, before the rest of the file name. The marker keeps the stare-mode products apart from the nod-mode products of the same SOF. The recipe writes them to the `soxs-stare` or `soxs-stare-std` product directory, not to the `soxs-nod` directory.

| Product | Nod-mode name | Stare-mode name |
| ------- | ------------- | --------------- |
| Merged spectrum | `<SOF>_SPECTRUM.fits` | `<SOF>_STAREMODE_SPECTRUM.fits` |
| Response curve of a standard star | `<SOF>_RESP.fits` | `<SOF>_STAREMODE_RESP.fits` |
| Product of one frame (`--no-stack`) | - | `<SOF>_F<N>_STAREMODE_SPECTRUM.fits` |
| Sky-subtracted frame of one nod position | - | `<SOF>_A_STAREMODE_SKYSUB.fits` and `<SOF>_B_STAREMODE_SKYSUB.fits` |
| Extraction of one nod position | - | `<SOF>_A_STAREMODE_SPECTRUM_A.fits` and `<SOF>_B_STAREMODE_SPECTRUM_B.fits` |

The recipe stacks the A frames and the B frames apart, so the sky products and the extraction of each nod position carry the position before the marker. The combined spectrum has no position in its name. In `--no-stack` mode, the frame number comes before the marker. The names of the SOF in the database, in the QC metrics and in the status of the SOF do not carry the marker. The data organiser does not predict the stare-mode products, so a later SOF never uses them.
