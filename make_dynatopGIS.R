## script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- './dynatopGIS'
Rcpp::compileAttributes(pacPath)
devtools::document(pacPath)
devtools::check(pacPath)
tmp <- devtools::build(pacPath)
install.packages(tmp)
pkgdown::clean_site(pacPath)
pkgdown::build_site(pacPath)


## This code runs to generate the model used in the dynatop examples
devtools::load_all(pacPath)
dem_file <- system.file("extdata", "SwindaleDTM4mFilled.tif", package="dynatopGIS")
channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp", package="dynatopGIS")
ctch <- create_catchment(dem_file)
shp <- rgdal::readOGR(channel_file)
property_names <- c(channel_id="identifier",
                    endNode="endNode",
                    startNode="startNode",
                    length="length")
chn <- create_channel(shp,property_names)
ctch <- add_channel(ctch,chn)

ctch <- sink_fill(ctch)
ctch <- compute_properties(ctch)
cuts <- list(atanb=20)
stck <- split_to_class(stck,'atb_split',cuts)
model <- create_model(stck,chn,'atb_split')
saveRDS(model,"./dynatop/inst/extdata/Swindale.rds")





