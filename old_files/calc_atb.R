#' function to fill sinks
#'
#' @description Compute the Topographic Index
#' 
#' @param project_path
#' @param use_filled_dem logical, shouldthe filled DEM be used?
#' @param max_iterations Maximum number of iterations
#'
#' @details The algorithm computes the area draining to a cell from those upstream presuming the flow is split between the downsteam (lower)  neighbouring cells proportional to the gradient.
#' @export
calc_atb <- function(project_path,use_filled_dem=TRUE,max_iter=10000){

    ## load dem and channel id
    file_list <- list()
    if(use_filled_dem){
        file_list[[1]] <- file.path(project_path,'filled_dem.tif')
    }else{
        file_list[[1]] <- file.path(project_path,'dem.tif')
    }
    file_list[[2]] <- file.path(project_path,'channel_id.tif')
    file_list[[3]] <- file.path(project_path,'land_area.tif')

    brck <- raster::brick(file_list)
    names(brck) <- c('dem','channel_id','land_area')

    ## check it is projected and square
    if(isLonLat(brck)){
        stop("Currently only works for projected DEM")
    }

    ## compute av_tanb
    #av <- tanb <- focal(brck[['dem']]
    dist <- matrix(sqrt(raster::xres(brck)^2 + raster::yres(brck)^2),3,3)
    dist[2,1] <- dist[2,3] <- raster::xres(brck)
    dist[1,2] <- dist[3,2] <- raster::yres(brck)
    dist[2,2] <- 0
    #dist <- as.vector(dist)
 
    ## create blank output files
    mat_atb <- fun_atb(as.matrix(brck[['dem']]),
                       as.matrix(brck[['land_area']]),
                       is.finite(as.matrix(brck[['channel_id']])),
                       dist, max_iter)

    
    atb <- raster::raster(brck,layer=0)
    atb <- setValues(atb,mat_atb)
    raster::writeRaster( atb, file.path(project_path,'atb.tif') )

    return(TRUE)
}

    
