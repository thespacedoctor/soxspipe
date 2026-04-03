

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
        soxspipe [-Vxd] stare <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vxd] nod <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe watch (start|stop|status) [-s <pathToSettingsFile>]
    
    Options:
        list ob                                list all observations within the workspace
        list sof                               list all science object SOF files within the workspace
        prep                                   prepare a folder of raw data (workspace) for data reduction
        session ls                             list all available data-reduction sessions in the workspace
        session new [<sessionId>]              start a new data-reduction session, optionally give a name up to 16 characters A-Z, a-z, 0-9 and/or _-
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
    
        start                                   start the watch daemon
        stop                                    stop the watch daemon
        status                                  print the status of the watch daemon
    
        inputFrames                            path to a directory of frames or a set-of-files file
    
        -b, --batch                            reduce data in batches of <batchSize> recipes (only when reducing all data)
        -d, --debug                            show debugging plots
        -h, --help                             show this help message
        -m, --multiprocess                     run reductions of recipe in parallel (experimental, use with caution and check your results carefully if using this flag)
        -p, --prep                             prepare a workspace before reducing data
        -q, --quitOnFail                       stop the pipeline if a recipe fails
        -r, --refresh                          trigger a complete refresh the workspace during preparation (delete database and do a complete prepare)
        -s, --settings <pathToSettingsFile>    the settings file
        -v, --version                          show version
        -V, --verbose                          more verbose output
        -w, --watch                            watch the workspace and reduce new raw data as it is added (similar to 'watch' mode but runs in the foreground)
        -x, --overwrite                        more verbose output
        --poly=<ORDERS>                        polynomial degrees (overrides parameters found in setting file). oowwss = order_x,order_y,wavelength_x,wavelength_y,slit_x,slit_y e.g. 345435. od = order,dispersion-axis
        --vlt                                  only use this flag if setting up a workspace on a VLT environment workstation
    

```
