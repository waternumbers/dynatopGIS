rm(list=ls())
setwd('/home/paul/Documents/Software/dynamic_topmodel/')
pacPath <- "./dynatopGIS"
devtools::load_all(pacPath)
stck <- raster::stack("debugging.gri")
chn <- readRDS("chn_debug.rds")
model <- create_model(stck,chn,'atb_split')
saveRDS(model,"model_debug.rds")
tmp <- band_model(model,10)
saveRDS(tmp,"band_model_debug.rds")

#library.dynam.unload("dynatopGIS", libpath = file.path(pacPath,"src"))

