## script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- './dynatopGIS'
## Rcpp::compileAttributes(pacPath)
devtools::document(pacPath)
devtools::check(pacPath)
tmp <- devtools::build(pacPath)
install.packages(tmp)
pkgdown::clean_site(pacPath)
pkgdown::build_site(pacPath)


## This code runs to generate the model used in the dynatop examples
## it also checks the verbose mode code
rm(list=ls())
pacPath_ <- './dynatopGIS'
devtools::load_all(pacPath)
##library("dynatopGIS")
dem <- raster::raster(system.file("extdata", "SwindaleDTM4mFilled.tif", package="dynatopGIS"))
shp <- rgdal::readOGR(system.file("extdata", "SwindaleRiverNetwork.shp", package="dynatopGIS"))
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

saveRDS(mdl,"./dynatop/inst/extdata/Swindale_from_r6.rds")

