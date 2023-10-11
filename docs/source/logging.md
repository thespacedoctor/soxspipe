# Logging

When running a recipe, `soxspipe` writes informative logs to the terminal (stdout), allowing the user to keep track of the reduction progress in real time. For the sake of provenance, this same information is written to a log file adjacent to the recipe's product file(s).

[![](https://live.staticflickr.com/65535/53246177840_c710373c18_z.png)](https://live.staticflickr.com/65535/53246177840_c710373c18_o.png)

If the recipe happens to fail, a separate error log is written to the directory the product file should have been written to had the recipe succeeded. Error logs are named with a *"_ERROR.log"* suffix.

[![](https://live.staticflickr.com/65535/53246182095_dae0648eb0_z.png)](https://live.staticflickr.com/65535/53246182095_dae0648eb0_o.png)
