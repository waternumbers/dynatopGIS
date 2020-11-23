## script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- './dynatopGIS'
devtools::document(pacPath)
devtools::check(pacPath)

## check documentation build
pkgdown::clean_site(pacPath)
pkgdown::build_site(pacPath)
pkgdown::clean_site(pacPath)

## build, populate drat
## linux
dratPath <- "~/Documents/Software/drat"
tmp <- devtools::build(pacPath)
install.packages(tmp)
drat::insertPackage(tmp,dratPath)

## mac and windows
rhub::validate_email() # for first time that session
pkgName <- sub('\\.tar.gz$', '', basename(tmp)) 
## rhub::platforms()[,1] # lists platforms
mch <- rhub::check(path = tmp,
                   platform = c("macos-highsierra-release-cran","windows-x86_64-release"))

tmp <- paste0(pkgName,".tgz")
download.file(file.path(mch$urls()$artifacts[1],tmp),tmp)
drat::insertPackage(tmp,dratPath)

tmp <- paste0(pkgName,".zip")
download.file(file.path(mch$urls()$artifacts[2],tmp),tmp)
drat::insertPackage(tmp,dratPath)

## tidy up drat
drat::pruneRepo(dratPath,pkg=pkgName,remove="git")## this only does source files

############################################################################

## This code runs to generate the model used in the dynatop examples
## it also checks the verbose mode code
rm(list=ls())
#pacPath <- './dynatopGIS'
#devtools::load_all(pacPath)
library("dynatopGIS")
dem <- raster::raster(system.file("extdata", "SwindaleDTM4mFilled.tif", package="dynatopData"))
shp <- rgdal::readOGR(system.file("extdata", "SwindaleRiverNetwork.shp", package="dynatopData"))
property_names <- c(channel_id="identifier",
                     endNode="endNode",
                     startNode="startNode",
                     length="length")

profvis::profvis({
    c2 <- dynatopGIS$new(dem)
    c2$add_channel(shp,property_names)
    c2$sink_fill(verbose=TRUE)
    c2$compute_properties(verbose=TRUE)
    c2$classify(list(atanb=20))
    c2$create_model(verbose=TRUE)
    mdl <- c2$get_model()
    c3 <- borg(c2)
})

saveRDS(mdl,"./dynatop/inst/extdata/Swindale.rds")





