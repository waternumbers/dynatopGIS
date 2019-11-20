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

    ## make sure all values in the dem area are finite
    dem_raster <- file.path(project_path,'dem.tif')
    dem <- raster(dem_raster)
    na_clumps <- clump(is.na(dem)) # clumps of na values
    edge_values <- setdiff( unique(c(na_clumps[1,],
                                     na_clumps[nrow(tmp),],
                                     na_clumps[,1],
                                     na_clumps[,ncol(tmp)])),
                           NA) #those clumps on the edge to be ignored
    na_clumps[na_clumps%in%edge_values] <- NA # set to NA to ignore

    for(ii in setdiff( unique(tmp), NA)){
        print(paste("Fill clump",ii))
        print("whoops not implimented")
    }

    ## check all non river cells have at least one lower cell
    f_less <- function(x){ ifelse(is.na(x[5]),NA,sum(x[-5]<x[5],na.rm=TRUE)) }
    f_fill <- function(x){
        x <- x[-5]
        if(all(is.na(x))){
            out <- NA
        }else{
            x <- sort(x)
            out <- ifelse(length(x)==1,x+1e-3,0.5*(x[1]+x[2]))
        }
        return(out)
    }
    
    
    n_less <- focal(dem,w = matrix(1,3,3),fun=f_less)
    idx <- Which(n_less==0,cells=TRUE)
    df <- dem
    while(length(idx) > 0){
        print(paste("Filling",length(idx),"cells"))
        df[idx] <- NA
        df <- focal(df,w = matrix(1,3,3),fun=f_fill,NAonly=TRUE )
        n_less <- focal(df,w = matrix(1,3,3),fun=f_less)
        idx <- Which(n_less==0,cells=TRUE)
    }
    

    
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
   

    ## remove na value from dem
    
     values f
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
    
