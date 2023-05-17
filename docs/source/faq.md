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
