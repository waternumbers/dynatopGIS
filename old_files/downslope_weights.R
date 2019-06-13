#' function to compute the downslope weights
#'
#' @description Blah...
#' 
#' @param project_path folder whcih is being used for the analysis
#' @param use_filled_dem logical indicating if the analysis should start with the filled DEM
#'
#' @details the algotithm works by looping over all the cells in the DEM and is relavitley robust but inefficent. If the maximum number of iterations is exceeded and there are still sinks the \code{use_filled_dem} flag can be used to restart from the final state.
#' @export
downslope_weights <- function(project_path,use_filled_dem=TRUE){

    ## load dem and channel id
    if(use_filled_dem){
        dem <- raster::raster( file.path(project_path,'filled_dem.tif') )
    }else{
        dem <- raster::raster( file.path(project_path,'dem.tif') )
    }
    channel <- raster::raster( file.path(project_path,'channel_id.tif') )

    ## check it is square
    if( raster::xres(dem)!=raster::yres(dem) ){
        stop("Upstream area only works with dem with equal x & y resolutuions")
    }

    ## standardise dem so don't have to pass distances into gradient
    dem <- dem / raster::xres(dem)

    ## weight for each direction
    weight_calc <- function(x,...){
        s2 <- sqrt(2)
        pmax((x[5]-x)/c(s2,1,s2,1,1,s2,1,s2),0)
    }

    delta <- raster::brick(
                         (dem - raster::focal(dem,w=matrix(c(0,0,1,0,0,0,0,0,0),3),
                                             pad=TRUE))/sqrt(2),
                         dem - raster::focal(dem,w=matrix(c(0,0,0,0,0,1,0,0,0),3),
                                             pad=TRUE),
                         (dem - raster::focal(dem,w=matrix(c(0,0,0,0,0,0,0,0,1),3),
                                              pad=TRUE))/sqrt(2),
                         dem - raster::focal(dem,w=matrix(c(0,1,0,0,0,0,0,0,0),3),
                                             pad=TRUE),
                         dem-dem,
                         dem - raster::focal(dem,w=matrix(c(0,0,0,0,0,0,0,1,0),3),
                                             pad=TRUE),
                         (dem - raster::focal(dem,w=matrix(c(1,0,0,0,0,0,0,0,0),3),
                                             pad=TRUE))/sqrt(2),
                         dem - raster::focal(dem,w=matrix(c(0,0,0,1,0,0,0,0,0),3),
                                             pad=TRUE),
                         (dem - raster::focal(dem,w=matrix(c(0,0,0,0,0,0,1,0,0),3),
                                              pad=TRUE))/sqrt(2)
                         )

    f_weight <- function(x){
        x <- pmax(x,0)
        sx <- sum(x,na.rm=TRUE) # zero if all NA
        if( sx == 0 ){
            x
        }else{
            x/sx
        }
    }

    
    ## calculate the weights
    weight <- calc(delta,fun=f_weight)
    ## set weight to channel reaches to be NA
    weight[ Which(channel>0) ] <- NA

    ## write output
    raster::writeRaster(weight, file.path(project_path,'downslope_weight.tif') ,overwrite=TRUE)
    return(TRUE)
}
