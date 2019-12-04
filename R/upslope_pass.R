#' Function to compute upslope area
#'
#' @description Computes a logical raster which identifies sinks
#' 
#' @param dem File or raster containing dem data
#'
#' @details The algorithm computes sinks at all locations where it is the
#' lowest cell
#' @export
upslope_pass <- function(brck,analysis_file=brck){

    if(!("RasterBrick" %in% class(brck))){
        brck <- raster::brick(brck)
    }else{
        analysis_file=""
    }
    

    ## initialise the order raster
    order <- brck[["channel_id"]]
    order[] <- NA
    order[!is.na(brck[["channel_id"]])&!is.na(brck[["dem"]])] <- 1

    ## create other dist to go to call
    offset <- c(-ncol(dem) + -1:1,-1,1,ncol(dem) + -1:1)
    
    out <- fun_upslope_pass(as.vector(dem), as.vector(order), offset)

    for(ii in c("dem_filled","order")){
        brck[[ii]] <- brck[['dem']];
        brck[[ii]][] <- out[[ii]]
    }
    
    if(length(analysis_file)>0){
        writeRaster(brck,paste0(analysis_file))
    }
    
    return(brck)
}
