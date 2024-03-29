## worked example for changes to v03
rm(list=ls())
graphics.off()

devtools::load_all("../")
unlink( list.files("./demo",full.names=TRUE) )


demo_file <- "./demo"

ctch <- dynatopGIS$new(demo_file)
dem_file <- system.file("extdata", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
ctch$add_dem(dem_file)
ctch$get_layer()
#ctch$plot_layer("dem")

channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp", package="dynatopGIS", mustWork = TRUE)

chn <- convert_channel(channel_file,c(name="identifier",
                                      endNode="endNode",
                                      startNode="startNode",
                                      length="length"))

ctch$add_channel(chn)

ctch$get_layer()
ctch$plot_layer("dem")

ctch$sink_fill()
ctch$get_layer()

ctch$compute_properties()
ctch$get_layer()

ctch$compute_flow_lengths()
ctch$get_layer()

tmp <- ctch$get_layer("filled_dem")
tmp <- terra::classify( tmp,
                          matrix(c(0,500,NA,
                                   500,1000,-999),
                                 ncol=3,byrow=TRUE))
ctch$add_layer(tmp,"greater_500")
ctch$get_layer()




ctch$classify("atb_20","atb",cuts=20)
ctch$plot_layer("atb_20")
ctch$get_layer()

ctch$get_method("atb_20")

ctch$combine_classes("atb_20_band",c("atb_20","band_inc_chn"))
ctch$get_method("atb_20_band")

ctch$combine_classes("atb_20_band_500",pairs=c("atb_20_band"),burns="greater_500")
tmp <- ctch$get_method("atb_20_band_500")$groups

ctch$plot_layer("atb_20_band_500")


ctch$create_model("./demo/test_model","atb_20_band","band_inc_chn",verbose=TRUE)

#unlink( list.files("./demo",full.names=TRUE) )
devtools::load_all("../../dynatop")
tmp <- readRDS("./demo/test_model.rds")
dy <- dynatop$new(tmp$hru)
data(Swindale)
dy$add_data(Swindale$obs)
dy$initialise()
dy$sim(tmp$output_flux)
out <- dy$get_output()
plot(Swindale$obs[,"flow"]); lines(out,col="red")

