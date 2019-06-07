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
