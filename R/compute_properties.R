#' Function to compute statistics of catchment pixels
#'
#' @description Computes statistics e.g. log(a/tanb) for raster cells
#'
#' @param stck RasterStack or file containing the raster stack (as generated using create_catchment)
#' @param ... additional parameters passed to raster::stack if stck is a file name
#'
#' @details The algorithm works in two passes. The first computes the number of upstream pixels. The second sequences downslope to compute values.
#' @export
compute_properties <- function(stck,...){

    if(!("RasterStack" %in% class(stck))){
        if( is.character(stck) ){
            stck_file <- stck
            stck <- raster::stack(stck,...)
        }else{
            stop("Unknown format for input")
        }
    }else{
        stck_file <- character(0)
    }

    check_catchment(stck)

    ## check it is projected
    if(raster::isLonLat(stck)){
        stop("Currently only works for projected DEM")
    }

    ## create other things to go with the call
    offset <- c(-ncol(stck) + -1:1,-1,1,ncol(stck) + -1:1)

    ## distance in each direction
    dx <- rep(sqrt(raster::xres(stck)^2 + raster::yres(stck)^2),8)
    dx[c(2,7)] <- raster::yres(stck)
    dx[c(4,5)] <- raster::xres(stck)

    ## contour length
    if( raster::xres(stck) != raster::yres(stck) ){
        warning("Contour length presumes a square grid")
    }
    mres <- (raster::xres(stck) + raster::yres(stck))/2
    cl <- c(rep( mres /(1+sqrt(2)),8),mres) # TO DO this is based on a n octogan - but other papers return a different ratio

    ## call cpp code
    out <- rcpp_compute_properties(raster::getValues(stck[["filled_dem"]]),
                                   raster::getValues(stck[["channel_id"]]),
                                   raster::getValues(stck[["land_area"]]),
                                   offset,
                                   dx,
                                   cl)


    for(ii in names(out)){
        idx <- which( names(stck)==ii )
        stck <- raster::setValues(stck, out[[ii]], layer=idx)
    }

    if(length(stck_file)>0){
        raster::writeRaster(stck,stck_file)
        return(stck_file)
    }else{
        return(stck)
    }

}
