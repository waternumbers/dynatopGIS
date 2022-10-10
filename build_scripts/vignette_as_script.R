## ----load_library-------------------------------------------------------------
rm(list=ls())
devtools::load_all()
##library("dynatopGIS")
unlink(demo_dir)

## ----tempory_dir--------------------------------------------------------------
demo_dir <- tempfile("dygis")
dir.create(demo_dir)


## ---- initialization----------------------------------------------------------
ctch <- dynatopGIS$new(file.path(demo_dir,"meta.json"))


## ---- data_files--------------------------------------------------------------
dem_file <- system.file("extdata", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp", package="dynatopGIS", mustWork = TRUE)


## ---- add_dem-----------------------------------------------------------------
dem <- terra::rast(dem_file)
ctch$add_dem(dem)


## ---- channel_current---------------------------------------------------------
sp_lines <- terra::vect(channel_file)
head(sp_lines)


## ---- channel_properties------------------------------------------------------
property_names <- c(channel_id="identifier",
                    endNode="endNode",
                    startNode="startNode",
                    length="length")


## ---- add_channel-------------------------------------------------------------
ctch$add_channel(sp_lines,property_names)


## ----basic_properties---------------------------------------------------------
ctch$compute_areas()


## ---- list_layers-------------------------------------------------------------
ctch$get_layer()


## ---- plot--------------------------------------------------------------------
ctch$plot_layer("dem", add_channel=TRUE)


## ---- get_layer---------------------------------------------------------------
ctch$get_layer("dem")


## ---- get_meta----------------------------------------------------------------
tmp <- ctch$get_meta()


## ---- sink_fill---------------------------------------------------------------
ctch$sink_fill()

raster::plot( ctch$get_layer('filled_dem') - ctch$get_layer('dem'),
             main="Changes to height")


## ---- calc_atb----------------------------------------------------------------
ctch$compute_properties()


## ---- plot_atb----------------------------------------------------------------
## plot of topographic index (log(a/tan b))
ctch$plot_layer('atb')


## ----extract_filled-----------------------------------------------------------
tmp <- ctch$get_layer("filled_dem")


## ----height layer-------------------------------------------------------------
tmp <- raster::reclassify( tmp,
                          matrix(c(0,500,NA,
                                   500,1000,-999),
                                 byrow=TRUE))


## ---- write_height_layer------------------------------------------------------
raster::writeRaster(tmp,file.path(demo_dir,"greater_500.tif"))


## ---- add_height_layer--------------------------------------------------------
ctch$add_layer("greater_500",file.path(demo_dir,"greater_500.tif"))
ctch$get_layer()


## ---- flow_length-------------------------------------------------------------
ctch$compute_flow_lengths()


## ---- flow_length_plot--------------------------------------------------------
ctch$get_layer()
ctch$plot_layer("band")


## ---- atb_split---------------------------------------------------------------
ctch$classify("atb_20","atb",cuts=20)
ctch$plot_layer("atb_20")


## ----atb_splt_get_class-------------------------------------------------------
ctch$get_class_method("atb_20")


## ---- atb_20_band-------------------------------------------------------------
ctch$combine_classes("atb_20_band",c("atb_20","band"))
ctch$plot_layer("atb_20_band")


## ---- atb_20_band_burn--------------------------------------------------------
ctch$combine_classes("atb_20_band_500",pairs=c("atb_20","band"),burns="greater_500")
ctch$plot_layer("atb_20_band_500")


## ----see_class----------------------------------------------------------------
head( ctch$get_class_method("atb_20_band_500") )


## ---- model_atb_split---------------------------------------------------------
ctch$create_model("new_model","atb_20_band","band")


## ---- model files-------------------------------------------------------------
list.files(demo_dir,pattern="new_model*")


## ----plot_model---------------------------------------------------------------
ctch$plot_layer("new_model")

