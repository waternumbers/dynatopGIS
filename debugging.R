rm(list=ls())
pacPath <- "./dynatopGIS"
devtools::load_all(pacPath)
stck <- raster::stack("debugging.gri")
chn <- readRDS("chn_debug.rds")
model <- create_model(stck,chn,'atb_split')

#library.dynam.unload("dynatopGIS", libpath = file.path(pacPath,"src"))

