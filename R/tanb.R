#' Function to compute tanb
#'
#' @description Computes the average downslope gradient for use in topographic index calculations
#' 
#' @param project_path folder which is being used for the analysis
#' @param use_filled_dem logical, shouldthe filled DEM be used?
#'
#' @details The algorithm computes the area draining to a cell from those upstream presuming the flow is split between the downsteam (lower)  neighbouring cells proportional to the gradient.
#' @export
tanb <- function(project_path,use_filled_dem=TRUE){

    ## load dem and channel id
    file_list <- list()
    if(use_filled_dem){
        file_list[[1]] <- file.path(project_path,'filled_dem.tif')
    }else{
        file_list[[1]] <- file.path(project_path,'dem.tif')
    }
    file_list[[2]] <- file.path(project_path,'channel_id.tif')
    
    brck <- raster::brick(file_list)
    names(brck) <- c('dem','channel_id')

    ## check it is projected
    if(isLonLat(brck)){
        stop("Currently only works for projected DEM")
    }

    ## compute av_tanb
    dist <- matrix(sqrt(raster::xres(brck)^2 + raster::yres(brck)^2),3,3)
    dist[2,1] <- dist[2,3] <- raster::xres(brck)
    dist[1,2] <- dist[3,2] <- raster::yres(brck)
    dist[2,2] <- 0
 
    ## create blank output files
    mat_tanb <- fun_tanb(as.matrix(brck[['dem']]),
                         is.finite(as.matrix(brck[['channel_id']])),
                         dist)
   
    tanb <- raster::raster(brck,layer=0)
    tanb <- raster::setValues(tanb,mat_tanb)
    raster::writeRaster(tanb, file.path(project_path,'tanb.tif') )

    return(TRUE)
}

    
