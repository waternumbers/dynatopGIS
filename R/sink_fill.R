#' function to fill sinks
#'
#' @description Fills sinks so that they have a slope of at least the minimum_tangent to the downstream cell
#' 
#' @param project_path folder whcih is being used for the analysis
#' @param max_iter maximum number of iterations of the algorithm
#' @param minimum_tangent  minimum value of teh gradient to the next cell downstream
#' @param use_filled_dem logical indicating if the analysis should start with the filled DEM
#'
#' @details the algotithm works by looping over all the cells in the DEM and is relavitley robust but inefficent. If the maximum number of iterations is exceeded and there are still sinks the \code{use_filled_dem} flag can be used to restart from the final state.
#' @export
sink_fill<- function(project_path,max_iter=10,minimum_tangent=0.001,use_filled_dem=FALSE){

    ## load dem and channel id
    if(use_filled_dem){
        dem <- raster::raster( file.path(project_path,'filled_dem.tif') )
    }else{
        dem <- raster::raster( file.path(project_path,'dem.tif') )
    }
    channel <- raster::raster( file.path(project_path,'channel_id.tif') )

    
    ## convert to correct type
    mat_dem <- as.matrix(dem)
    is_channel <- is.finite( as.matrix(channel) )

    ## compute distance
    dist <- matrix(sqrt(raster::xres(dem)^2 + raster::yres(dem)^2),3,3)
    dist[2,1] <- dist[2,3] <- raster::xres(dem)
    dist[1,2] <- dist[3,2] <- raster::yres(dem)
    dist[2,2] <- 0
    dist <- dist*minimum_tangent
    
    mat_dem <- fun_sink_fill(mat_dem,is_channel, dist, as.integer(max_iter) )
    dem <- raster::setValues(dem,mat_dem)
    raster::writeRaster( dem, file.path(project_path,'filled_dem.tif') )

    return(TRUE)
}
    
