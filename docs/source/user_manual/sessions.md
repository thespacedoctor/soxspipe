## Data Reduction Sessions

soxspipe reduces data within a 'reduction session'. These sessions are designed to be self-contained and isolated from other sessions, allowing users to reduce the same set of raw data with multiple different pipeline settings. An initial ' base ' session is automatically created when a workspace is first prepared (using the `soxspipe prep` command). Sessions are stored in the `sessions` directory of the workspace and contain their own settings file, products and QC directories. soxspipe will remember and use the latest session you were working with unless you create a new session or switch to another one. 

To see which reduction session you are working within and all other sessions available, run the command:

```bash
soxspipe session ls
```

To create a new session, run the command:

```bash
soxspipe session new
```

The new session will be named with the date-time stamp of its creation date. Alternatively, you can also supply a memorable name for the session. The name can be up to 16 characters long and use alpha-numeric characters, \- or \_.

```bash
soxspipe session new my_supernova
```

Or to switch to an existing session, run the command:

```bash
soxspipe session <sessionId>
```

For user convenience, when you switch to a new session, the symbolic links found within the workspace root folder are automatically switched to point to the current session assets (`products`, `qc`, `sof`, `soxspipe.yaml`,`soxspipe.db` etc). If you run `ls -lrt` within your workspace root directory, you will see these symlinks reported:

```bash
products -> ./sessions/base/products
qc -> ./sessions/base/qc
soxspipe.db -> ./sessions/base/soxspipe.db
soxspipe.yaml -> ./sessions/base/soxspipe.yaml
sof -> ./sessions/base/sof
```

