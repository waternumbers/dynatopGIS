# dynatopGIS

This R package contains the code to generate a dynamic TOPMODEL from GIS
data. A related package for computing dynamic TOPMODEL simulations can be found
[here](https://waternumbers.github.io/dynatop).

These two packages are intended as successors to the original Dynamic TOPMODEL implementation in R;
the [dynatopmodel](https://CRAN.R-project.org/package=dynatopmodel) package. Currently
(as a design choice) not all the features of the original package are
implemented. 

If there is a feature or improvement you want or you experience problems
please raise an [issue](https://github.com/waternumbers/dynatopGIS/issues)
(or better still submit a pull request with alterations or fixes!) rather than contact the
maintainers directly in the first instance.

Thanks goes to:
* Peter Metcalfe who developed the original port of Dynamic TOPMODEL to R
during his PhD sponsored by the [JBA Trust](https://www.jbatrust.org). 
* The original developers of Dynamic TOPMODEL, [Keith
Beven](https://www.lancaster.ac.uk/lec/about-us/people/keith-beven) and [Jim
Freer](http://www.bristol.ac.uk/geography/people/jim-e-freer/index.html), who
contributed to the original R port of Dynamic TOPMODEL.

## Using the code

Currently the packages are not available on
[CRAN](https://cran.r-project.org/). They can be installed from within R using
the devtools package: 

```
devtools::install_github("waternumbers/dynatopGIS")
```

## Future developments
Two key future developments are:
* To support the use of large (out of memory) GIS files
* To improve iterative algorithms (iterations should not be needed in many
cases!!)
