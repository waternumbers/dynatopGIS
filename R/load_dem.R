#' Initialise the analysis by loading a dem
#'
#' @description Reads the DEM file, sets up the raster catchment description
#'
#' @param dem must be either a rasterLayer or a file name
#' @param catchment_file filename to save output to
#' @param fill_na should NA values in dem be filled. See details
#' @param ... passed to raster package functions to set options for reading files
#' @param stck a RasterStack or file name for testing to see if it is a catchment description
#'
#' @return A raster stack object with initialised layers or the name of file into which is is saved
#'
#' @details Reads in the dem using the raster package. If fill_na all NA values other then those that link to the edge of the dem are filled with -Inf, so they can be identified as sinks.
#' @name catchment_file

#' @rdname catchment_file
#' @export
create_catchment <- function(dem,catchment_file="",fill_na=TRUE,...){

    ## read in dem
    if(!("RasterLayer" %in% class(dem))){
        if( is.character(dem) ){
            dem <- raster::raster(dem,...)
        }else{
            stop("Unknown dem format")
        }
    }

    ## check projection is plausible
    check_projection(dem)

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
        dem[!is.na(na_clumps)] <- Inf # set to high value to indicate missing
    }

    ## initialise the catchment file
    stck_names <- req_catchment_layers()
    tmp <- dem
    values(tmp) <- NA
    stck <- raster::stack( rep(list(tmp),length(stck_names)) )
    names(stck) <- stck_names

    ## set dem
    idx <- which(stck_names=="dem")
    stck <- setValues(stck, getValues(dem), layer=idx)

    ## compute land area
    idx <- which(stck_names=="land_area")
    land_area <- tmp
    land_area <- setValues(land_area,prod(raster::res(stck)))
    land_area <- raster::mask( land_area, stck[['dem']] )
    stck <- setValues(stck, getValues(land_area), layer=idx)

    ## write out
    if(catchment_file!="" && length(catchment_file)>0){
        writeRaster(stck,catchment_file)
        return(catchment_file)
    }else{
        return(stck)
    }
}


#' Function to check the stack
#' @rdname catchment_file
#' @export
check_catchment <- function(stck,...){
    if(!is(stck,"RasterStack")){
        if( is.character(stck) ){
            stck <- raster::stack(stck,...)
        }else{
            stop("Unknown catchment format")
        }
    }
    req_names <- req_catchment_layers()
    if(!is(stck,"RasterStack")){stop("Not a raster brick")}
    if(!all( req_names %in% names(stck) )){stop("Catchment missing key variables")}

    check_projection(stck)
}

#' Function to check the projection of the stack is OK
#' @rdname catchment_file
#' @export
check_projection <- function(stck,...){
    if(!is(stck,"RasterStack")){
        if( is.character(stck) ){
            stck <- raster::stack(stck,...)
        }else{
            stop("Unknown catchment format")
        }
    }

    ## check it is projected
    if(raster::isLonLat(stck)){
        stop("Processing currently only works for projected DEM")
    }

    ## check it is  a square grid
    if( raster::xres(stck) != raster::yres(stck) ){
        stop("Processing currently only works on square gridded DEMs")
    }

}

#' returns the list of required layers in the catchment file
req_catchment_layers <- function(){
    c("dem","filled_dem","land_area","channel_area","channel_id",
      "atanb","gradient","upslope_area","contour_length","order")
}
