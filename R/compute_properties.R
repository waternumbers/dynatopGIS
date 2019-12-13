#' Function to compute statistics of brck pixels
#'
#' @description Computes statistics e.g. log(a/tanb) for raster cells
#' 
#' @param brck rasterBrick or file containing teh raster brick
#'
#' @details The algorithm works in two passes. the first computed the number fo upstream pixels. The second sequences downslope to compute values.
#' @export
compute_properties <- function(brck,...){

    if(!("RasterBrick" %in% class(brck))){
        if( is.character(brck) ){
            brck <- raster::brick(brck,...)
        }else{
            stop("Unknown brck format")
        }
    }

    check_brick(brck)
    
    ## check it is projected
    if(raster::isLonLat(brck)){
        stop("Currently only works for projected DEM")
    }

    ## create other things to go with the call
    offset <- c(-ncol(brck) + -1:1,-1,1,ncol(brck) + -1:1)

    ## distance in each direction
    dx <- rep(sqrt(raster::xres(brck)^2 + raster::yres(brck)^2),8)
    dx[c(2,7)] <- raster::yres(brck)
    dx[c(4,5)] <- raster::xres(brck)
    if( raster::xres(brck) != raster::yres(brck) ){
        warning("Contour length presumes a square grid")
    }

    ## contour length
    mres <- (raster::xres(brck) + raster::yres(brck))/2
    cl <- c(rep( mres /(1+sqrt(2)),8),mres) # TO DO this is based on a n octogan - but other papers return a different ratio

    ## call cpp code
    out <- fun_single_pass(as.vector(brck[["filled_dem"]]),
                              as.vector(brck[["channel_id"]]),
                              as.vector(brck[["land_area"]]),
                              offset,
                              dx,
                              cl)
    

    for(ii in names(out)){
        brck[[ii]][] <- out[[ii]]
    }
    
    return(brck)
}
