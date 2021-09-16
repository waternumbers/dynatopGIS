## script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- '.'#/dynatopGIS'
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
download.file(file.path(mch$urls()$artifacts[1],tmp),file.path("..",tmp))
drat::insertPackage(file.path("..",tmp),dratPath)

tmp <- paste0(pkgName,".zip")
download.file(file.path(mch$urls()$artifacts[2],tmp),file.path("..",tmp))
drat::insertPackage(file.path("..",tmp),dratPath)

## tidy up drat
drat::pruneRepo(dratPath,pkg=pkgName,remove="git")## this only does source files

############################################################################

## This code runs to generate the model used in the dynatop examples
## it also checks the verbose mode code
rm(list=ls())
pacPath <- '.'# ##/dynatopGIS'
devtools::load_all(pacPath)
#library("dynatopGIS")
dem <- raster::raster(system.file("extdata", "SwindaleDTM4mFilled.tif", package="dynatopData"))
shp <- rgdal::readOGR(system.file("extdata", "SwindaleRiverNetwork.shp", package="dynatopData"))
agg_dem <- raster::aggregate(dem,10)

property_names <- c(channel_id="identifier",
                     endNode="endNode",
                     startNode="startNode",
                     length="length")

profvis::profvis({
    demo_dir <- tempfile("Swindale")
    dir.create(demo_dir)
    c2 <- dynatopGIS$new(file.path(demo_dir,"meta.json"))
    c2$add_dem(agg_dem)
    c2$add_channel(shp,property_names)
    c2$compute_areas()
    c2$sink_fill(verbose=TRUE)
    c2$compute_properties(verbose=TRUE)
    c2$compute_flow_lengths(verbose=TRUE)
    c2$classify("atb_20","atb",20)
    c2$combine_classes("atb_20_band3",c("atb_20","band"))
    c2$create_model("Swindale_exp","atb_20_band","band",transmissivity="exp",verbose=TRUE)
    c2$create_model("Swindale_bexp","atb_20_band","band",transmissivity="bexp",verbose=TRUE)
    c2$create_model("Swindale_dexp","atb_20_band","band",transmissivity="dexp",verbose=TRUE)
    c2$create_model("Swindale_cnst","atb_20_band","band",transmissivity="cnst",verbose=TRUE)
})

file.copy(list.files(demo_dir,pattern="^Swindale.*\\.rds$", full.names = TRUE),
          "../dynatop/build_scripts/",TRUE)
file.copy(file.path(demo_dir,"Swindale.tif"),"../dynatop/inst/extdata/",TRUE)
file.copy(file.path(demo_dir,"channel_id.tif"),"../dynatop/inst/extdata/",TRUE)
file.copy(list.files(demo_dir, "^channel[.]", full.names = TRUE),
          "../dynatop/inst/extdata/",TRUE)



