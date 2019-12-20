#' Function to compute upslope area
#'
#' @description Computes a logical raster which identifies sinks
#' 
#' @param dem File or raster containing dem data
#'
#' @details The algorithm computes sinks at all locations where it is the
#' lowest cell
#' @export
downslope_pass <- function(brck,analysis_file=brck){

    if(!("RasterBrick" %in% class(brck))){
        brck <- raster::brick(brck)
    }else{
        analysis_file=""
    }
    
    ## check it is projected
    if(isLonLat(brck)){
        stop("Currently only works for projected DEM")
    }
    
    ## work out the sequence
    ord <- as.vector(brck[["order"]])
    sq <- (0:(length(ord)-1))[order(ord,decreasing=TRUE,na.last=NA)]
    
    ## create other things to go with the call
    offset <- c(-ncol(dem) + -1:1,-1,1,ncol(dem) + -1:1)
    dx <- rep(sqrt(raster::xres(brck)^2 + raster::yres(brck)^2),8)
    dx[c(2,7)] <- raster::yres(brck)
    dx[c(4,5)] <- raster::xres(brck)
    if( raster::xres(brck) != raster::yres(brck) ){
        warning("Contour length presumes a square grid")
    }
    mres <- (raster::xres(brck) + raster::yres(brck))/2
    cl <- c(rep( mres /(1+sqrt(2)),8),mres) # TO DO this is based on a n octogan - but other papers return a different ratio

    ## call cpp code
    out <- fun_downslope_pass(as.vector(brck[["dem"]]), as.vector(brck[["order"]]),
                              as.vector(sq),
                              offset,
                              as.vector(brck[['land_area']]),
                              dx,
                              cl)
    

    for(ii in c("upslope_area","contour_length","gradient","atanb")){
        brck[[ii]] <- brck[['dem']];
        brck[[ii]][] <- out[[ii]]
    }
    
    if(length(analysis_file)>0){
        writeRaster(brck,paste0(analysis_file))
    }
    
    return(brck)
}
