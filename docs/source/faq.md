# FAQs

- My conda install is out of date. How do I update it?

    You may see a warning like this. 

    ```bash
    ==> WARNING: A newer version of conda exists. <==
      current version: 4.9.2
      latest version: 23.3.1
    ```

    If so, run the following command:

    ```bash
    conda install -n base -c defaults conda=23.3.1
    ```

- Can I export all the raw frames needed to reduce a specific dataset?

    Yes. Within your workspace directory, first list all of the science SOF files:

    ```bash
    soxspipe list sof 
    ```

    Then select the file you are interested in from the list of sof files generated.

    Now run the `raw` command to export all of the associated raw frames to an `exported` folder in the workspace directory. For example

    ```bash
    soxspipe raw sof 20260128T075112_NIR_3_STARE-OBJ_SLIT1_0_600_0S_SOXS_PAQS135626+042348.sof
    ```

- My spectrum seems to be very noisy in the bluer region < 400 nm. Why is this?

    The QTH lamp SOXS currently uses to generate the flat-field frames produces little signal in this region. Therefore, when the data is divided by the master flat, this additional noise is introduced. This can also affect the quality of the response curve generated from the spectroscopic standard stars in the same region and, therefore, the final flux calibration of science sources.
