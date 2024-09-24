## Terminology

`soxspipe` employs a 'recipes' concept to describe each discrete, modular data-reduction step. Each of these 'recipes' accepts a set of files (SOF) file containing a list of the input files, or ingredients, required by the recipe. The recipes are then strung together to form a complete data-reduction cascade (or cookbook if you like), which converts raw science frames into ESO Phase 3-compliant science data products. 

The term 'utility' defines the reusable functions called from multiple recipes. Recipes are named with the prefix 'soxs' followed by a succinct recipe description (e.g. ``soxs_mbias``).

