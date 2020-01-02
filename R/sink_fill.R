#' Fill the sinks with a dem
#'
#' @description TO DO
#'
#' @param stck a rasterBrick object or whatever
#' @param ... additional varaibles passed to raster::brick if stck is a file name
#'
#' @return a rasterBrick object
#'
#' @details TO DO
#'
#' @export
sink_fill <- function(stck,...){

    if(!("RasterStack" %in% class(stck))){
        if( is.character(stck) ){
            stck <- raster::brick(stck,...)
        }else{
            stop("Unknown stck format")
        }
    }

    warning("Sink filling is not implimented...this just copies the dem")
    stck[['filled_dem']] <- stck[['dem']]

    return(stck)
}
