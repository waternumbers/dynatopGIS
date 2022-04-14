## worked example for changes to v03
rm(list=ls())
graphics.off()

devtools::load_all("../")
unlink( list.files(".",pattern="^demo"))


demo_file <- "demo.nc"

ctch <- dynatopGIS$new(demo_file)
dem_file <- system.file("extdata", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
ctch$add_dem(dem_file)
ctch$get_layer()

channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp", package="dynatopGIS", mustWork = TRUE)

chn <- convert_channel(channel_file,c(name="identifier",
                                      endNode="endNode",
                                      startNode="startNode",
                                      length="length"))
ctch$add_channel(chn)
ctch$get_layer()

ctch$sink_fill()
ctch$get_layer()

ctch$compute_properties()
ctch$get_layer()
