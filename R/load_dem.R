#' Initialise the analysis by loading a dem
#'
#' @description Reads the DEM file, sets up the raster catchment description
#'
#' @param dem must be either a rasterLayer or a file name
#' @param catchment_file filename to save output to
#' @param fill_na should NA values in dem be filled. See details
#' @param ... passed to raster package functions to set options for reading files
#' @param stck a RasterStack or file name for testing to see if it is a catchment description
#'
#' @return A raster stack object with initialised layers or the name of file into which is is saved
#'
#' @details Reads in the dem using the raster package. If fill_na all NA values other then those that link to the edge of the dem are filled with -Inf, so they can be identified as sinks.
#' @name catchment_file

#' @rdname catchment_file
#' @export
create_catchment <- function(dem,fill_na=TRUE,...){

    ## read in dem
    if(!("RasterLayer" %in% class(dem))){
        if( is.character(dem) ){
            dem <- raster::raster(dem,...)
        }else{
            stop("Unknown dem format")
        }
    }


    ## fill dem so only NA values are on boundaries
    if(fill_na){
        ## so all catchment dem value a number
        na_clumps <- raster::clump(is.na(dem)) # clumps of na values
        edge_values <- setdiff( unique(c(na_clumps[1,],
                                         na_clumps[nrow(na_clumps),],
                                         na_clumps[,1],
                                         na_clumps[,ncol(na_clumps)])),
                               NA) #those clumps on the edge to be ignored
        na_clumps[na_clumps%in%edge_values] <- NA # set to NA to ignore
        dem[!is.na(na_clumps)] <- -Inf # set to low value to indicate missing
    }

    ## initialise the catchment object
    catchment <- raster_to_glist(dem,"dem")

    ## add land area
    catchment$layers$land_area <- rep(prod(catchment$res),length(catchment$layers$dem))*
        !is.na(catchment$layers$dem)

    ## check projection is plausible
    check_projection(catchment)

    return(catchment)
}


#' Function to check the stack
#' @rdname catchment_file
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
#' @rdname catchment_file
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
      "atanb","gradient","upslope_area","order","flowDir")
}
