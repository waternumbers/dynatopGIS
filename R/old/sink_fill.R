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

#rm(list=ls())
#library(raster)
    project_path <- "./dynatopGIS/vignette_data/processed"

    ## load dem
    dem_raster <- file.path(project_path,'dem.tif')
    dem <- raster(dem_raster)
    
    ## load channel
    channel_raster <- file.path(project_path,'channel_id.tif')
    channel <- raster(channel_raster)
        
    ## fill all NA clumps that aren't on the edge with -Inf
    ## so all catchment dem value a number
    na_clumps <- clump(is.na(dem)) # clumps of na values
    edge_values <- setdiff( unique(c(na_clumps[1,],
                                     na_clumps[nrow(na_clumps),],
                                     na_clumps[,1],
                                     na_clumps[,ncol(na_clumps)])),
                           NA) #those clumps on the edge to be ignored
    na_clumps[na_clumps%in%edge_values] <- NA # set to NA to ignore
    dem[!is.na(na_clumps)] <- -Inf
    
    ## work out sinks
    sink <- is_sink(dem)
    sink[!is.na(channel)] <- FALSE
    
    ## set up order raster
    depth <- channel
    depth[] <- NA
    depth[!is.na(channel)&!is.na(dem)] <- 1
    
    ## depth layer to look at
    for(dp in 1:1000){
        print(paste("Looking at depth order",dp))
        #browser()
        idx <- Which(depth==dp,cells=TRUE) # idx of all cells of order
        print(length(idx))
        if(length(idx)==0){
            break
        }
        
        #adj <- raster::adjacent(dem,idx,directions=8, pairs=FALSE) # adjacent cell indexs - includes current cells
        #adj <- adj[(adj%in%idx)]
        
        ## check no adjacent cells are sinks
        #if( any(sink[adj]>0) ){
        #    stop("Not written yet")
        #}
        
        ## find adjacent cells that are greater then and add increased order
        order <- increase_order(dem,depth)
        depth[!is.na(channel)&!is.na(dem)] <- 1
    }
}
