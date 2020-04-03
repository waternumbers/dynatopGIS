#' Functions to help handling and checking catchments
#'
#' @description A mixture of provate and public functions that help in handling catchment objects
#'
#' @param ctch A catchment
#'
#' @details At some point in the futrue this will become a class...
#'
#' @name catchment
#'
#'
#' #' Function to check the stack
#' @rdname catchment
#' @export
check_catchment <- function(ctch,req_names=req_catchment_layers()){
    if(!is_glist(ctch)){
        stop("Unknown catchment format")
    }

    if(!all( req_names %in% names(ctch$layers) )){
        stop(paste(c("Catchment missing key layers:",
                     setdiff( req_names,  names(ctch$layers))),
                   collapse="\n"))
    }
    
    check_projection(ctch)
}

#' Function to check the projection of the stack is OK
#' @rdname catchment
#' @export
check_projection <- function(ctch,...){
    if(!is_glist(ctch)){
        stop("Check requires a glist object")
    }

    ## check it is  a square grid
    if( ctch$raster$res[1] != ctch$raster$res[2] ){
        stop("Processing currently only works on square gridded DEMs")
    }
    
    ## check it is projected
    if(raster::isLonLat(glist_to_raster(ctch,"dem"))){
        stop("Processing currently only works for projected DEM")
    }
    
}

#' returns the list of required layers in the catchment file
req_catchment_layers <- function(){
    c("dem","filled_dem","land_area","channel_area","channel_id",
      "atanb","gradient","upslope_area","band","flowDir")
}

#' plot a catchment
plot_catchment <- function(ctch,property="dem",add_channel=FALSE){
    property <- property[1]
    
    has_channel <- "channel" %in% names(ctch)
    all_prop <- c("channel"[has_channel], names(ctch$layers))

    if( property %in% all_prop ){
        if( property == "channel" ){
            raster::plot(ctch$channel)
        }else{
            raster::plot(glist_to_raster(ctch,property))
        }
    }else{
        stop("Selected property not available")
    }

    if( add_channel ){
        if( has_channel ){ raster::plot(ctch$channel,add=TRUE) }
        else{ warning("Channel is not present to be added") }
    }
}

## convert catchment to a gis layer
catchment_to_gis <- function(ctch,file_name,property="dem"){
    property <- property[1]
    
    has_channel <- "channel" %in% names(ctch)
    all_prop <- c("channel"[has_channel], names(ctch$layers))
    
    if( property %in% all_prop ){
        if( property == "channel" ){
            rgdal::writeOGR(ctch$channel,file_name)
        }else{
            raster::writeRaster(glist_to_raster(ctch,property),file_name)
        }
    }else{
        stop("Selected property not available")
    }
}

    
