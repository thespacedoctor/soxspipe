

```bash
    
    Documentation for soxspipe can be found here: http://soxspipe.readthedocs.org
    
    Usage:
        soxspipe --version
        soxspipe prep [<workspaceDirectory> --vlt --refresh]
        soxspipe [-qwpVm] reduce all [<workspaceDirectory> -b <batchSize> -s <pathToSettingsFile>]
        soxspipe [-qxV] reduce sof <sofFile> [<workspaceDirectory> -s <pathToSettingsFile>]
        soxspipe session ((ls|new|<sessionId>)|new <sessionId>)
        soxspipe list (ob|sof) [<workspaceDirectory> -s <pathToSettingsFile>]
        soxspipe raw sof <sofFile> [<workspaceDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vxd] mdark <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vxd] mbias <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vxd] disp_solution <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<ooww>]
        soxspipe [-Vxd] order_centres <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<ooww>]
        soxspipe [-Vxd] mflat <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vxd] spat_solution <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<oowwss>]
        soxspipe [-Vxd] (stare|stare_std) <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vxd] (nod|nod_std) <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vxd] offset <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe watch (start|stop|status) [-s <pathToSettingsFile>]
    
    Options:
        list ob                                list all observations within the workspace
        list sof                               list all science object SOF files within the workspace
        prep                                   prepare a folder of raw data (workspace) for data reduction
        session ls                             list all available data-reduction sessions in the workspace
        session new [<sessionId>]              start a new session, name it (A-Z, a-z, 0-9 and/or _; 16 chars max)
        session <sessionId>                    use an existing data-reduction session (use `session ls` to see all IDs)
        reduce all                             reduce all of the data in a workspace.
        reduce sof                             reduce a single science object SOF file.
        raw sof                                export all the raw frames needed to reduce a science object SOF file to a directory called `exported` in the current working directory.
    
        mbias                                  the master bias recipe
        mdark                                  the master dark recipe
        mflat                                  the master flat recipe
        disp_solution                          the disp solution recipe
        order_centres                          the order centres recipe
        spat_solution                          the spatial solution recipe
        stare                                  reduce stare mode science frames
        nod                                    reduce nodding mode science frames
        stare_std                              reduce stare mode standard-star frames
        nod_std                                reduce nodding mode standard-star frames
        offset                                 reduce offset mode science frames
    
        start                                   start the watch daemon
        stop                                    stop the watch daemon
        status                                  print the status of the watch daemon
    
        inputFrames                            path to a directory of frames or a set-of-files file
    
        -b, --batch                            reduce data in batches of <batchSize> recipes (only when reducing all data)
        -d, --debug                            show debugging plots
        -h, --help                             show this help message
        -m, --multiprocess                     run reductions of recipe in parallel (experimental, use with caution and check your results carefully if using this flag)
        -o, --output <outputDirectory>         the output directory for the recipe product
        -p, --prep                             prepare a workspace before reducing data
        -q, --quitOnFail                       stop the pipeline if a recipe fails
        -r, --refresh                          trigger a complete refresh the workspace during preparation (back up the database to `backups/`, rebuild it keeping its QC history and the pass/fail status of each unchanged SOF, and do a complete prepare)
        -s, --settings <pathToSettingsFile>    the settings file
        -v, --version                          show version
        -V, --verbose                          more verbose output
        -w, --watch                            watch the workspace and reduce new raw data as it is added (similar to 'watch' mode but runs in the foreground)
        -x, --overwrite                        more verbose output
        --poly=<ORDERS>                        polynomial degrees (overrides parameters found in setting file). oowwss = order_x,order_y,wavelength_x,wavelength_y,slit_x,slit_y e.g. 345435. od = order,dispersion-axis
        --vlt                                  only use this flag if setting up a workspace on a VLT environment workstation
    

```

## Refresh a workspace

Use `soxspipe prep --refresh` to rebuild the workspace database from the raw frames. This section is for people who run `soxspipe prep` on an existing workspace.

The same rebuild starts by itself when `soxspipe prep` finds a workspace database that will not open.

### What is kept

Before the database is deleted, `soxspipe` saves it to the `backups/` directory in the workspace root. If the database cannot be saved, `soxspipe prep` leaves it in place, rebuilds nothing, and exits with status 1.

After the rebuild, `soxspipe` restores two things from the backup:

- The `quality_control` rows.
- The pass/fail status of each SOF (set-of-files) file, for every session directory in the workspace, not only the current session.

The QC acceptable-range checks run again after each restore. A restored `fail` status changes to `pass` if the QC rows for that product are now in range.

### When a status is restored

A status is restored only when both of these are true:

- The SOF name exists in the rebuilt database.
- The SOF holds the same set of files before and after the rebuild.

If the set of files is different, the status is not restored and the SOF is queued for reduction again. A SOF that exists only in the backup is dropped. The backup file keeps its status.

### Summary lines

`soxspipe prep` prints one line for each session. For example:

```text
session 'science': 212 statuses restored, 2 changed (frames differ, requeued), 3 dropped (product no longer exists)
```

When the changed count or the dropped count is more than zero, the line ends with `; the previous statuses are kept in <backup path>`.

### Limits

- The restore is best effort. If `soxspipe` cannot read the backup, it prints a warning that names the backup file, leaves the status columns empty, and continues. It does not ask you anything.
- For a session that is not the current session, `soxspipe` compares the SOF against the SOF map that it rebuilt for the current session.
- The `error_message` column is not restored.
