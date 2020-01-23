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

    check_catchment(stck)
    
    ## for testing
    #stck <- raster::stack("./Swindale_stck.gri")

    ## extract dem and channel_id values
    dem <- raster::getValues(stck[['dem']])
    ch_id <- raster::getValues(stck[['channel_id']])
    
    ## work out lowest neighbour
    nc <- ncol(stck)
    nr <- nrow(stck)
    min_neighbour <- rep(NA,length(dem))
    idx <- which(!is.na(dem))
    for(ii in idx){
        min_neighbour[ii] <- min(dem[fN(ii,nr,nc)])
    }
    
    ## determine sinks and set to Inf - this loop should only increase dem value
    ## sink should be lower then neighbours and not have a channel
    idx <- which( min_neighbour >= dem && !is.finite(ch_id) )
    dem[idx] <- Inf
    while(length(idx)>0){
        io <- fN(idx[1],nr,nc) # neighbours of first point
        for(ii in io){
            ## no point checking if already NA or Inf or a channel
            if(is.finite(dem[ii]) & !is.finite(ch_id[ii]) ){
                min_neighbour[ii] <- min(dem[fN(ii,nr,nc)])
                if(!is.na(min_neighbour[ii]) && (min_neighbour[ii] >= dem[ii])){
                    idx <- c(idx,ii)
                    dem[ii] <- Inf
                }
            }
        }
        idx <- idx[-1]
    }

    ## populate starting with lowest
    idx <- which(dem==Inf)
    idx <- idx[order(min_neighbour[idx])]
    
    for(ii in idx){
        dd <- dem[fN(ii,nr,nc)]
        dem[ii] <- mean(dd[is.finite(dd)])
    }
    
    ## copy filled dem back into stack
    stck <- raster::setValues(stck, dem, layer=which(names(stck)=="filled_dem"))
    
    if(length(stck_file)>0){
        raster::writeRaster(stck,stck_file)
        return(stck_file)
    }else{
        return(stck)
    }

}
