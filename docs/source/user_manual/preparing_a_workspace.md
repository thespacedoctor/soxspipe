# Preparing a Data-Reduction Workspace

A 'workspace' is simply a directory you create and designate as the location on your machine where you will reduce your data. It is advisable to have multiple workspaces (folders) for multiple datasets. 

First, create your workspace (you can call it whatever you want):

```bash
cd ~/
mkdir soxs-workspace-one
cd soxs-workspace-one
```

Now download/move all of the SOXS<!--or Xshooter--> raw data from your data set into the root of this workspace. There is no need for nested folders; just fill the directory with FITS files.

 `soxspipe` has a built-in [data organiser](../data_organiser_and_reducer/data_organiser.md) that does the heavy lifting of sorting and preparing your raw data for data reduction. So to prepare your workspace, change into the workspace directory and run the `soxspipe prep` command:

```bash
cd ~/soxs-workspace-one
soxspipe prep .
```

Once the workspace has been prepared, you should find it contains the following files and folders:

- `misc/`: a lost-and-found archive of non-fits files
- `qc/`: a directory to host the quality control data and plots
- `raw/`: all raw frames to be reduced
- `reduced/`: a directory to host the reduced data
- `sessions/`: directory of data-reduction sessions
- `sof/`: the set-of-files (sof) files required for each reduction step
- `soxspipe.db`: an SQLite database needed by the data-organiser; please do not delete
- `soxspipe.yaml`: file containing the default settings for each pipeline recipe

Your workspace is now prepared, and you are now ready to reduce your data.

:::{include} sessions.inc
:::

## Update the settings of an existing workspace

Each time you run `soxspipe prep` on a workspace, it brings the `soxspipe.yaml` of the active session up to date with the default settings of the installed `soxspipe`. The file is written once, when the session is created, so without this step a change to a default never reaches an old workspace.

`soxspipe prep` changes the file before it reads the settings. It does two things:

1. Superseded values. The package keeps a list of default values that a later release replaced. If a setting in your file has exactly the old value of one of these, `soxspipe prep` writes the new default. For example, the NIR spatial solution `order-deg` of `[4, 4]` becomes `[5, 4]`. A value that differs from both the old and the new default is a choice you made, and it stays as it is.
2. Missing keys. If the default settings file has a key that your file does not have, at any depth, `soxspipe prep` adds it with its default value.

`soxspipe prep` does not remove keys. A key that the defaults no longer have stays in your file. `soxspipe prep` never adds a per-workspace value, such as the log file path, and never changes one.

`soxspipe prep` writes one log line for each change, with the setting path, the old value (or `missing`), the new value, and the ticket for a superseded value. It also prints one line with the counts, for example `14 superseded and 138 missing settings changed`. If nothing needs to change, it prints nothing (the counts go to the debug log only) and does not rewrite the file. If a setting in your file is empty or is not a group of settings where the defaults have a group, `soxspipe prep` leaves it alone and logs a warning with the setting path.

Two limits apply:

- The file is rewritten without its comments. `soxspipe prep` does not make a backup copy of the old file.
- Only the active session is changed. Other sessions keep their files until you switch to them and run `soxspipe prep`.

To list the superseded values, read `soxspipe/resources/settings_migrations.yaml` in the installed package. Each entry has the setting path, the old value, the new value, and the ticket. Workspace tools can call `migrate_session_settings(yamlPath, instrument, log)` from `soxspipe.commonutils.settings_migration` to do the same on any `soxspipe.yaml`.

## Listing The Content of a Workspace

### List OBs

The `soxspipe list` commands can be used to list the science content of the workspace. To list all of the science OBs contained within a workspace, run the command `soxspipe list obs <workspaceDirectory>`. 

:::{figure-md} soxspipe_list_obs
![image-20251119123700080](../_images/image-20251119123700080.png)

Listing the OBs contained within a workspace with the `soxspipe list obs .` command.
:::

### List SOF files

An OB can contain multiple telescope-instrument setups (slit-widths, observing modes etc), but each individual setup will be listed in a unique SOF file. To list all of the science set-of-files (SOF) files contained within a workspace, run the command `soxspipe list sof <workspaceDirectory>`. 

:::{figure-md} soxspipe_list_obs
![image-20251119124816171](../_images/image-20251119124816171.png)

Listing the science SOF files contained within a workspace with the `soxspipe list sof .` command.
:::

## Refresh a workspace

Use `soxspipe prep --refresh` to rebuild the workspace database from the raw frames and do a complete prepare. This section is for people who run `soxspipe prep` on an existing workspace.

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
