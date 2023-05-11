

```bash
    
    Documentation for soxspipe can be found here: http://soxspipe.readthedocs.org
    
    Usage:
        soxspipe init
        soxspipe prep <workspaceDirectory>
        soxspipe [-Vx] mbias <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>] 
        soxspipe [-Vx] mdark <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] disp_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] order_centres <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] mflat <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] spat_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe [-Vx] stare <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
    
    Options:
        init                                   setup the soxspipe settings file for the first time
        prep                                   prepare a folder of raw data (workspace) for data reduction
        mbias                                  the master bias recipe
        mdark                                  the master dark recipe
        mflat                                  the master flat recipe
        disp_sol                               the disp solution recipe
        order_centres                          the order centres recipe
        spat_sol                               the spatial solution recipe
        stare                                  reduce stare mode science frames
    
        inputFrames                            path to a directory of frames or a set-of-files file
    
        -h, --help                             show this help message
        -v, --version                          show version
        -s, --settings <pathToSettingsFile>    the settings file
        -V, --verbose                          more verbose output
        -x, --overwrite                        more verbose output
    

```
