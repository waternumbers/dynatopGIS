#' Function to compute upslope area
#'
#' @description Computes a logical raster which identifies sinks
#' 
#' @param dem File or raster containing dem data
#'
#' @details The algorithm computes sinks at all locations where it is the
#' lowest cell
#' @export
increase_order <- function(dem,order){
    if(!("RasterLayer" %in% class(dem))){
        dem <- raster::raster(dem)
    }
    if(!("RasterLayer" %in% class(dem))){
        order <- raster::raster(order)
    }
   
    ## check it is projected and square
    if(isLonLat(dem)){
        stop("Currently only works for projected DEM")
    }

    ## create blank output files
    order[] <- fun_increase_order(as.vector(dem),as.vector(order),ncol(dem))
        
    return(order)
}

    
