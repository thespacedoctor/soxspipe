

```bash
    
    Documentation for soxspipe can be found here: http://soxspipe.readthedocs.org
    
    Usage:
        soxspipe prep <workspaceDirectory>
        soxspipe [-qw] reduce all <workspaceDirectory> [-s <pathToSettingsFile>]
        soxspipe session ((ls|new|<sessionId>)|new <sessionId>)
        soxspipe [-Vx] mdark <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] mbias <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] disp_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<ooww>]
        soxspipe [-Vx] order_centres <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<ooww>]
        soxspipe [-Vx] mflat <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] spat_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<oowwss>]
        soxspipe [-Vx] stare <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] nod <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe watch (start|stop|status) [-s <pathToSettingsFile>]
    
    Options:
        prep                                   prepare a folder of raw data (workspace) for data reduction
        session ls                             list all available data-reduction sessions in the workspace
        session new [<sessionId>]              start a new data-reduction session, optionally give a name up to 16 characters A-Z, a-z, 0-9 and/or _-
        session <sessionId>                    use an existing data-reduction session (use `session ls` to see all IDs)
        reduce all                             reduce all of the data in a workspace.
    
        mbias                                  the master bias recipe
        mdark                                  the master dark recipe
        mflat                                  the master flat recipe
        disp_sol                               the disp solution recipe
        order_centres                          the order centres recipe
        spat_sol                               the spatial solution recipe
        stare                                  reduce stare mode science frames
        nod                                    reduce nodding mode science frames
    
        start                                   start the watch daemon
        stop                                    stop the watch daemon
        status                                  print the status of the watch daemon
    
        inputFrames                            path to a directory of frames or a set-of-files file
    
        -q, --quitOnFail                       stop the pipeline if a recipe fails
        -h, --help                             show this help message
        -v, --version                          show version
        -s, --settings <pathToSettingsFile>    the settings file
        -V, --verbose                          more verbose output
        -x, --overwrite                        more verbose output
        -w, --watch                            watch the workspace and reduce new raw data as it is added (similar to 'watch' mode but runs in the foreground)
        --poly=<ORDERS>                        polynomial degrees (overrides parameters found in setting file). oowwss = order_x,order_y,wavelength_x,wavelength_y,slit_x,slit_y e.g. 345435. od = order,dispersion-axis
    

```
