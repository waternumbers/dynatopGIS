#' Fill the sinks with a dem
#'
#' @description TO DO
#'
#' @param brck a rasterBrick object or whatever
#' @param ... additional varaibles passed to raster::brick if brck is a file name
#'
#' @return a rasterBrick object
#'
#' @details TO DO
#' 
#' @export
sink_fill <- function(brck,...){

    if(!("RasterBrick" %in% class(brck))){
        if( is.character(brck) ){
            brck <- raster::brick(brck,...)
        }else{
            stop("Unknown brck format")
        }
    }

    warning("Sink filling is not implimented...this just copies the dem")
    brck[['filled_dem']] <- brck[['dem']]

    return(brck)
}
