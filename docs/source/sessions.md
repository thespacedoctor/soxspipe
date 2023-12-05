# Data Reduction Sessions

soxspipe reduces data within a 'reduction session'. These sessions are designed to be self-contained and isolated from other sessions, allowing a user to reduce the same set of raw data with multiple different pipeline settings. When a workspace is first prepared (using the `soxspipe prep` command), an initial `base` session is automatically created. Sessions are stored in the `sessions` directory of the workspace and contain their own settings file, products and QC directories. soxspipe will remember and use the latest session you were working with unless you create a new session or switch to another one. 

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
