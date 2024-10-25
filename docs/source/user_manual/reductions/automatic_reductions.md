# Automatic Reductions Using a Watched Workspace

With soxspipe, it is possible to create a 'watched' workspace. When data is added to a watched workspace, soxspipe will attempt to reduce the data automatically.

First, create the directory you are going to designate as the watched workspace:

```bash
mkdir ~/soxs_workspace
```

Now, change into that directory and ask soxspipe to 'watch' the directory:

```bash
cd ~/soxs_workspace
soxspipe watch start
```

The daemon is now running in the background, and you can check its status by running:

```bash
soxspipe watch status
```

If you start moving FITS files into this workspace, soxspipe will attempt to reduce the data automatically. Note that the pipeline is running in daemon mode, so nothing will be output to the terminal. However, the daemon log file can be found at `~/.config/soxspipe/daemon.log`, and individual recipe logs can be found in their usual location beside the products in the workspace products directory. 

The daemon will continue running even if you close the terminal with which you initiated it. To stop the daemon, run:

```bash
soxspipe watch stop
```
