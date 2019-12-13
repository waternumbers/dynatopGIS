#' Fill the sinks with a dem
#'
#' @description TO DO
#'
#' @param brck a rasterBrick object or whatever
#'
#' @return a rasterBrick object
#'
#' @details TO DO
#' 
#' @export
sink_fill <- function(brck,...){

    if(!("RasterBrick" %in% class(brck))){
        if( is.character(brck) ){
            brck <- raster::brick(brck,...)
        }else{
            stop("Unknown brck format")
        }
    }

    warning("Sink filling is not implimented...this just copies the dem")
    brck[['filled_dem']] <- brck[['dem']]

    return(brck)
}


    

    ## read in shp if required
    if( is.character(shp) ){
        shp <- rgdal::readOGR(shp,stringsAsFactors=FALSE)
    }

    ## check shp is a SpatialLinesDataFrame
    if(!is(shp,"SpatialLinesDataFrame") && !is(shp,"SpatialPolygonsDataFrame")){
        stop("shp is not a spatial lines or polygon data frame object (even when read in)")
    }

      
    ## check property_names is of correct type
    if(!is.vector(property_names) ||
       length(names(property_names))!=length(property_names) ){
        stop("property_names of wrong type")
    }

    ## check we have correct names
    if( !all( c("length","startNode","endNode") %in% names(property_names)) ){
        stop("A required property name is not specified")
    }

    ## check if there is an id feild which will be overwritten
    if( ("id" %in% names(property_names)) ){
        warning("The name id is reserved and will be overwritten")
        property_names <- property_names[-"id"]
    }

    ## see if all names exist
    if(!all( property_names %in% names(shp))){
        stop("Not all feilds names in property_names are in shp")
    }
    
    ## create channel object
    chn <- shp
    chn <- chn[,property_names]
    names(chn) <- names(property_names)
    chanel[['id']] <- 1:length(chn)

    if(!is(chn,"SpatialPolygonsDataFrame")){
        if(!("width" %in% names(chn))){
            warning("Modifying to spatial polygons using default width")
            chn[['width']] <- default_width
        }else{
            warning("Modifying to spatial polygons using provided width")
        }
        
        ## buffer to convert to spatial polygons object
        chn <- rgeos::gBuffer(chn, byid=TRUE, width=chn[['width']])
    }

    return(chn)
}


check_channel <- function(chn){
    req_names <- c("id","startNode","endNode","length")
    if(!is(chn,"SpatialPolygonsDataFrame")){stop("Not a SpatialPolygonsDataFrame")}
    if(!all( req_names %in% names(chn) )){stop("Channel missing key variables")}
}

##     raster::shapefile(buffered_chanel, file.path(project_path,'channel'))
    
##     chanel[['id']] <- id[cid]
    
##     ## add blank values to missing feilds
##     missing <- setdiff( c('startNode','endNode','length'), names(chanel) )
##     if(length(missing)>0){
##         warning("The following feilds are missing and blank values will be added: ",
##                 paste(missing,collapse=" "))
##         for(ii in missing){ chanel[[ii]] <- NA }
##     }

##     ## add default width
##     if( !("width" %in% names(chanel)) ){
##         chanel[['width']] <- default_width
##     }

##     ## buffer to convert to spatial polygons object
##     buffered_chanel <- rgeos::gBuffer(chanel, byid=TRUE, width=chanel[['width']])
##     raster::shapefile(buffered_chanel, file.path(project_path,'channel'))
##     return(TRUE)
## }
