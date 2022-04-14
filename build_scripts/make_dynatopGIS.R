## script for creating and building the dynatop R package
rm(list=ls())
graphics.off()

## path of the package
pacPath <- '../'
devtools::document(pacPath)
devtools::spell_check(pacPath,use_wordlist=TRUE)

devtools::check(pacPath)

## check documentation build
pkgdown::clean_site(pacPath)
pkgdown::build_site(pacPath)
pkgdown::clean_site(pacPath)

## build, populate drat
## linux
dratPath <- "~/Documents/Software/drat"
buildFile <- devtools::build(pacPath)
install.packages(buildFile)
drat::insertPackage(buildFile,dratPath)

## mac and windows
rhub::validate_email() # for first time that session
pkgName <- sub('\\.tar.gz$', '', basename(buildFile)) 
## rhub::platforms()[,1] # lists platforms
mch <- rhub::check(path = buildFile,
                   platform = c("macos-highsierra-release-cran","windows-x86_64-release",
                                "macos-m1-bigsur-release"))

ext <- c(".tgz",".zip",".tgz")
outPath <- dirname(buildFile)
for(ii in 1:2){ ## m1 not fixed yet in drat
    tmp <- paste0(pkgName,ext[ii])
    outFile <- file.path(outPath,tmp)
    #download.file(file.path(mch$urls()$artifacts[ii],tmp),outFile)
    drat::insertPackage(outFile,dratPath)
}

## tidy up drat
drat::pruneRepo(dratPath,pkg=pkgName,remove="git")## this only does source files

## prior to submission
mch <- rhub::check_for_cran(path = buildFile)
mch <- rhub::check_with_valgrind(path = buildFile)

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



############################################################################

## This code builds the inst extdata directory from the raw data files in dynatopData

rm(list=ls())
dem_file <- system.file("extdata", "SwindaleDTM4mFilled.tif", package="dynatopData", mustWork = TRUE)
raw_dem <- raster::raster(dem_file)
agg_dem <- raster::aggregate(raw_dem,10)
writeRaster(agg_dem,"./inst/extdata/SwindaleDTM40m.tif")

file.copy( list.files(dirname(dem_file),pattern="^SwindaleRiver",full.names=TRUE),
          "./inst/extdata/")

          
