#' Initialise the analysis by loading a dem
#'
#' @description Reads the DEM file, sets up the raster brick
#'
#' @param dem must be either a rasterLayer or a file name
#' @param fill_na should NA values in dem be filled. See details
#' @param ... passed to raster::raster
#' @param brck a RasterBrick for testing to match those created by create_brick
#'
#' @return raster brick object with initialised layers
#'
#' @details Reads in the dem using the raster package. Fi fill_na all NA values other then those that link to the edge of the dem are filled with -Inf, so they can be identified as sinks.
#' @name gis_brick

#' @rdname gis_brick
#' @export
create_brick <- function(dem,fill_na=TRUE,...){

    ## read in dem
    if(!("RasterLayer" %in% class(dem))){
        if( is.character(dem) ){
            dem <- raster::raster(dem,...)
        }else{
            stop("Unknown dem format")
        }
    }

    ## fill dem so only NA values are on boundaries
    if(fill_na){
        ## so all catchment dem value a number
        na_clumps <- raster::clump(is.na(dem)) # clumps of na values
        edge_values <- setdiff( unique(c(na_clumps[1,],
                                         na_clumps[nrow(na_clumps),],
                                         na_clumps[,1],
                                         na_clumps[,ncol(na_clumps)])),
                               NA) #those clumps on the edge to be ignored
        na_clumps[na_clumps%in%edge_values] <- NA # set to NA to ignore
        dem[!is.na(na_clumps)] <- -Inf
    }

    ## initialise the other raster outputs
    brck <- raster::brick(dem);
    names(brck) <- "dem"
    tmp <- dem
    tmp[] <- NA
    for(ii in c("filled_dem","land_area","channel_area","channel_id",
                "atanb","gradient","upslope_area","contour_length","order")){
        ##browser()
        ##print(ii)
        names(tmp) <- ii
        brck <- raster::addLayer(brck,tmp)
        ##dem
        ##brck[[ii]][] <- NA
    }
    brck <- raster::brick(brck) # since addLeyer makes it a stack...

    ## compute land area
    if( raster::isLonLat(brck) ){
        brck[['land_area']] <- raster::area(brck)
    }else{
        brck[['land_area']][] <- prod(raster::res(brck))
    }
    brck[['land_area']] <- raster::mask( brck[['land_area']], brck[['dem']] )
 
    return(brck)
}


#' Function to check the brick
#' @rdname gis_brick
#' @export
check_brick <- function(brck){
    req_names <- c("filled_dem","land_area","channel_area","channel_id",
                   "atanb","gradient","upslope_area","contour_length","order")
    if(!is(brck,"RasterBrick")){stop("Not a raster brick")}
    if(!all( req_names %in% names(brck) )){stop("Brick missing key variables")}
}
