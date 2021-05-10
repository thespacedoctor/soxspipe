

```bash
    
    Documentation for soxspipe can be found here: http://soxspipe.readthedocs.org
    
    Usage:
        soxspipe init
        soxspipe mbias <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>] 
        soxspipe mdark <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe disp_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe order_centres <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
        soxspipe spat_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
    
    Options:
        init                                   setup the soxspipe settings file for the first time
        mbias                                  the master bias recipe
        mdark                                  the master dark recipe
        disp_sol                               the disp solution recipe
        order_centres                          the order centres recipe
        spat_sol                               the spatial solution recipe
    
        inputFrames                            path to a directory of frames or a set-of-files file
    
        -h, --help                             show this help message
        -v, --version                          show version
        -s, --settings <pathToSettingsFile>    the settings file
    

```
