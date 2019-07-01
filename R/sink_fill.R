#' function to fill sinks
#'
#' @description Fills sinks so that they have a slope of at least the minimum_tangent to the downstream cell
#' 
#' @param project_path folder which is being used for the analysis
#' @param max_iter maximum number of iterations of the algorithm
#' @param minimum_tangent  minimum value of the gradient to the next cell downstream
#' @param use_filled_dem logical indicating if the analysis should start with the filled DEM. If set this option will overwrite the current filled dem with the final output
#'
#' @details the algotithm works by looping over all the cells in the DEM and is relavitley robust but inefficent. If the maximum number of iterations is exceeded and there are still sinks the \code{use_filled_dem} flag can be used to restart from the final state.
#' @export
sink_fill<- function(project_path,max_iter=1000,minimum_tangent=0.001,use_filled_dem=FALSE){

    
    ## load dem and channel id
    file_list <- list()
    if(use_filled_dem){
        file_list[[1]] <- file.path(project_path,'filled_dem.tif')
        overwrite_dem <- TRUE
    }else{
        file_list[[1]] <- file.path(project_path,'dem.tif')
        overwrite_dem <- FALSE
    }
    file_list[[2]] <- file.path(project_path,'channel_id.tif')

    brck <- raster::brick(file_list)
    names(brck) <- c('dem','channel_id')

    ## check it is projected and square
    if(isLonLat(brck)){
        stop("Currently only works for projected DEM")
    }
   
    
    ## convert to correct type
    mat_dem <- as.matrix(brck[['dem']])
    is_channel <- is.finite( as.matrix(brck[['channel_id']]) )

    ## compute distance
    dist <- matrix(sqrt(raster::xres(brck)^2 + raster::yres(brck)^2),3,3)
    dist[2,1] <- dist[2,3] <- raster::xres(brck)
    dist[1,2] <- dist[3,2] <- raster::yres(brck)
    dist[2,2] <- 0
    dist <- dist*minimum_tangent
    
    mat_dem <- fun_sink_fill(mat_dem,is_channel, dist, as.integer(max_iter) )
    dem <- raster::raster(brck,layer=0)
    dem <- setValues(dem,mat_dem)
    raster::writeRaster( dem, file.path(project_path,'filled_dem.tif'),
                        overwrite=overwrite_dem)

    return(TRUE)
}
    
