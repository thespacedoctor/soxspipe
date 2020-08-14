## vx.x.x - xxdatexx

* **FEATURE** detector lookup class added alongside yaml files to host detector specific parameters (rotation, science-pixels etc). Code has been updated to remove hard-wired detector values.
* inject a 'SOXSPIPE PRE' keyword with timestamp value into prepared frames
* check frames for 'SOXSPIPE PRE' keyword before preparing - raises exception if found
* moved basic input frame verifications to the `_base_recipe` - so not to repeat code
* added mixing of readout speeds to input frame verification checks
* added a cleanup method to remove intermediate file once recipe completes
* removed python 2.7 support - not feasible with CCDProc

## v0.2.0 - February 27, 2020

* **FEATURE** added keyword lookups - abstracting exact keyword names from code
