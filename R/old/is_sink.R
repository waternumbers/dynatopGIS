#' Function to compute upslope area
#'
#' @description Computes a logical raster which identifies sinks
#' 
#' @param dem File or raster containing dem data
#'
#' @details The algorithm computes sinks at all locations where it is the
#' lowest cell
#' @export
is_sink <- function(dem){

    if(!("RasterLayer" %in% class(dem))){
        dem <- raster::raster(dem)
    }
    
    ## check it is projected and square
    if(isLonLat(dem)){
        stop("Currently only works for projected DEM")
    }


    sflag <- dem
    sflag[] <- NA
    
    ## create blank output files
    sflag[] <- fun_is_sink(as.vector(dem),ncol(dem))
        
    return(sflag)
}

    
