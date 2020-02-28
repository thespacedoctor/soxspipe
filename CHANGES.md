## vx.x.x - xxdatexx

* inject a 'SOXSPIPE PRE' keyword with timestamp value into prepared frames
* check frames for 'SOXSPIPE PRE' keyword before preparing - raises exception if found
* moved basic input frame verifications to the `_base_recipe` - so not to repeat code
* added mixing of readout speeds to input frame verification checks

## v0.2.0 - February 27, 2020

* added keyword lookups - abstracting exact keyword names from code
