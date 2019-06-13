#' function to fill sinks
#'
#' @description Fills sinks so that they have a slope of at least the minimum_tangent to the downstream cell
#' 
#' @param project_path folder whcih is being used for the analysis
#' @param max_iter maximum number of iterations of the algorithm
#' @param minimum_tangent  minimum value of the gradient to the next cell downstream
#' @param use_filled_dem logical indicating if the analysis should start with the filled DEM
#'
#' @details the algotithm works by looping over all the cells in the DEM and is relavitley robust but inefficent. If the maximum number of iterations is exceeded and there are still sinks the \code{use_filled_dem} flag can be used to restart from the final state.
#' @export
sink_fill_R<- function(project_path,max_iter=10,minimum_tangent=0.001,use_filled_dem=FALSE){

    ## load dem and channel id
    if(use_filled_dem){
        dem <- raster::raster( file.path(project_path,'filled_dem.tif') )
    }else{
        dem <- raster::raster( file.path(project_path,'dem.tif') )
    }
    channel <- raster::raster( file.path(project_path,'channel_id.tif') )

    ## check it is square
    if( raster::xres(dem)!=raster::yres(dem) ){
        stop("Sink fill only works with dem with equal x & y resolutuions")
    }

    ## standardise dem so don't have to pass distances into focal function
    dem <- dem / raster::xres(dem)

    ## test for if it is a sinkgradient of scaled DEM
    grad_calc <- function(x,...){
        s2 <- sqrt(2)
        max((x[5]-x[-5])/c(s2,1,s2,1,1,s2,1,s2))
    }
    ## fill value for DEM
    fill_calc <- function(x,...){
        s2 <- sqrt(2)
        min( x[-5] + c(s2,1,s2,1,1,s2,1,s2) )
    }
        
        
    it <- 1
    while(it <= max_iter){
        ## compute the max gradient from each pixel
        max_grad <- raster::focal(dem,w=matrix(rep(1,9),3),fun=grad_calc,
                                  pad=TRUE)
        ## find all sinks
        is_sink <- raster::Which( max_grad<minimum_tangent,
                                 na.rm=TRUE,cells=TRUE) # na.rm ensure we don't fill boundaries
        
        ## remove sinks pixels that are partly channel
        is_sink <- is_sink[is.na(channel[is_sink])]
        ## print out some crud
        print(paste0("Iteration ",it,": ",length(is_sink)," sinks to fill"))

        if(length(is_sink)==0){
            break
        }
        
        ## loop to fill sinks
        dem[is_sink] <- NA
        dem <- raster::focal(dem,w=matrix(rep(1,9),3),fun=fill_calc,
                             pad=TRUE, NAonly=TRUE)

        ## increment iteration counter
        it <- it + 1
    }

    ## scale back raster
    dem <- dem * raster::xres(dem)

    ## write output
    raster::writeRaster( dem, file.path(project_path,'filled_dem.tif') ,overwrite=TRUE)
    return(TRUE)
}

##         && is.na(channel) )
    
##     ## convert to correct type
##     mat_dem <- as.matrix(dem)
##     is_channel <- is.finite( as.matrix(channel) )

##     ## compute distance
##     dist <- matrix(sqrt(raster::xres(dem)^2 + raster::yres(dem)^2),3,3)
##     dist[2,1] <- dist[2,3] <- raster::xres(dem)
##     dist[1,2] <- dist[3,2] <- raster::yres(dem)
##     dist[2,2] <- 0
##     dist <- dist*minimum_tangent
    
##     mat_dem <- fun_sink_fill(mat_dem,is_channel, dist, as.integer(max_iter) )
##     dem <- raster::setValues(dem,mat_dem)
##     raster::writeRaster( dem, file.path(project_path,'filled_dem.tif') )

##     return(TRUE)
## }
    
