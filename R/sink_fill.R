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
            stck_file <- stck
            stck <- raster::stack(stck,...)
        }else{
            stop("Unknown stck format")
        }
    }else{
        stck_file <- character(0)
    }
    
    offset <- c(-ncol(stck) + -1:1,-1,1,ncol(stck) + -1:1)
    out <- rcpp_sink_fill(raster::getValues(stck[["dem"]]),
                          raster::getValues(stck[["channel_id"]]),
                          offset)
    stck <- setValues(stck, out, layers==which(names(stck)=="filled_dem"))
    
    #warning("Sink filling is not implimented...this just copies the dem")
    #stck[['filled_dem']] <- stck[['dem']]

    if(length(stck_file)>0){
        raster::writeRaster(stck,stck_file)
        return(stck_file)
    }else{
        return(stck)
    }

}
