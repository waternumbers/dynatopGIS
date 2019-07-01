# dynatopGIS 0.0.2

## New features
- Added ability to continue a run of sink_fill by starting with previous
  filled dem 
  
## Bug fixes
- Model generation - missing sbar added and class of data frame columns altered
- Linking Rcpp function - alter NAMESPACE so should be included

## Other changes
- Vignette know includes running of model and burning in a class

# dynatopGIS 0.0.1

## Other changes
- Removed vignette data from package to reduce size. Now on github
- Minor spelling and format changes to vignette

# dynatopGIS 0.0.0.9000

## Context
This package is the result of an almost complete rewrite of the dynatopmodel package
formally on CRAN and the associated development code (not in the public
domain).

This package contains the code for extracting teh GIS summaries required for
dynamic TOPMODEL and helpers for inserting default parameters. The package [dynatop]{http://bbc.co.uk} contains the
code for performing simulations.

## New Features
- Change in approach to wrting all intermediate files to
      - produce a record of analysis
	  - allow for future developments with larger size GIS data
- Porting of sink_fill and atb calculations from orginal topmodel CRAN package
  to Rcpp
- Improved method of determining river network intersection with DEM from shape files
- Changes to sink filling and atb calculations to recognise burnt in river network
  and correctly split flow by gradient
- Alterations to output mew style dynamic TOPMODEL object
- Vignettes added documeting function usage

## Regressions of note
- No processing of distance histogram for river routing
