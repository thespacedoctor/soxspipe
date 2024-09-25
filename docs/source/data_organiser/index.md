# Data Organiser

The data organiser (DO) is the backbone of the pipeline. It is required to sort and prepare data within a workspace, predict which products will be produced during data reduction, write out all of the [set-of-files (SOF)](../utils/set_of_files.md) files required by each of the soxspipe recipes and keep track of all data-products generated during the reduction cascade. The DO also provides functionality to rewrite SOF files on the fly if a recipe fails to produce a product required by a future recipe (e.g. a master flat frame), switching out the failed product for the next-best product (e.g. the next master flat frame generated closest in time to the recipe data). Finally, on subsequent executions of the pipeline, the organiser prevents data from being re-reduced if the products already exist (unless the user chooses to override this feature).

The algorithm the DO uses to prepare a workspace is shown in {numref}`data_organiser_util`.

:::{figure-md} data_organiser_util
![](soxs_data_organiser.png){width=600px}

The algorithm used by the soxspipe data-organiser to prepare a workspace for data-reduction.
:::

At the heart of the DO is a SQLite database called `soxspipe.db`. Here, the organiser's bookkeeping is performed, recorded, and maintained.

The ESO Science Archive Facility delivers FITS data in a `.Z` compressed format. When running `soxspipe prep`, the DO first finds and uncompresses for any `.Z` compressed FITS frames within the workspace root. The DO then reads the FITS headers of all of the FITS frames in the workspace root and selects out the raw (unreduced) frames, record one entry per raw frame in the `raw_frames` table of `soxspipe.db`. The DO then moves these raw frames to a `raw_frames` directory within the workspace. Any remaining files are moved out of the workspace root and into a `misc` directory.

A sanity check is performed to ensure that the data in the `raw_frames` database table matches the data in the `raw_frames` directory. If frames have been removed from the file system, the corresponding records in the database table are deleted. Also, frames within the `raw_frames` directory missing from the database table are added.

The next step is for the DO to define all sets of raw frames that can be used to produce next-stage products (master bias, master flat, order location tables, etc). These sets are recorded sets in the `raw_frame_sets` database table. These raw frame sets derive the raw frame content for all possible SOF files, which are recorded in the `sof_map` database table, assigning a human-readable 'tag' (e.g. BIAS_UVB) to individual frames and maps the frames to named sof files. 

The initial set of SOF files in the `sof_map` table is used to predict the product files written when soxspipe recipes are executed on the SOF files. The expected product information is written to a `product_frames` database table. From this `product_frames` table, the products are assigned to SOF files later in the reduction cascade (recorded, again, in the `sof_map` table). 

Finally, all SOF files from the `sof_map`  table are written in a `sof` directory in the workspace root and ready to be used by the various soxspipe recipes during a data-reduction session.


## Utility API

:::{autodoc2-object} soxspipe.commonutils.data_organiser.data_organiser
:::

