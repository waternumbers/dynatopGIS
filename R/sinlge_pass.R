#' Function to compute upslope area
#'
#' @description Computes a logical raster which identifies sinks
#' 
#' @param dem File or raster containing dem data
#'
#' @details The algorithm computes sinks at all locations where it is the
#' lowest cell
#' @export
single_pass <- function(brck,analysis_file=brck){

    if(!("RasterBrick" %in% class(brck))){
        brck <- raster::brick(brck)
    }else{
        analysis_file=character(0)
    }
    
    ## check it is projected
    if(isLonLat(brck)){
        stop("Currently only works for projected DEM")
    }

    ## create other things to go with the call
    offset <- c(-ncol(brck) + -1:1,-1,1,ncol(brck) + -1:1)
    dx <- rep(sqrt(raster::xres(brck)^2 + raster::yres(brck)^2),8)
    dx[c(2,7)] <- raster::yres(brck)
    dx[c(4,5)] <- raster::xres(brck)
    if( raster::xres(brck) != raster::yres(brck) ){
        warning("Contour length presumes a square grid")
    }
    mres <- (raster::xres(brck) + raster::yres(brck))/2
    cl <- c(rep( mres /(1+sqrt(2)),8),mres) # TO DO this is based on a n octogan - but other papers return a different ratio

    ## call cpp code
    out <- fun_single_pass(as.vector(brck[["dem"]]),
                              as.vector(brck[["channel_id"]]),
                              as.vector(brck[["land_area"]]),
                              offset,
                              dx,
                              cl)
    

    for(ii in names(out)){
        brck[[ii]] <- brck[['dem']];
        brck[[ii]][] <- out[[ii]]
    }
    
    if(length(analysis_file)>0){
        writeRaster(brck,paste0(analysis_file))
    }
    
    return(brck)
}
