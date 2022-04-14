
convert_channel <- function(sp_object,property_names=c(name = "DRN_ID",
                                                       length = "length",
                                                       startNode = "startNode",
                                                       endNode = "endNode",
                                                       width = "width"),
                            default_width=2){
    
    ## read in sp object is a character sting
    if(is.character(sp_object)){
        if(file.exists(as.character(sp_object))){
            sp_object <- raster::shapefile(sp_object)
        }else{
            stop("sp_object is a character string but the file specified does not exist")
        }
    }

    ## find out about the sp_object
    if(!is(sp_object,"SpatialLinesDataFrame") && !is(sp_object,"SpatialPolygonsDataFrame")){
        stop("The channel network is not a spatial lines or polygon data frame object (even when read in)")
    }
    is_polygon <- is(sp_object,"SpatialPolygonsDataFrame")

    ## check all feilds in property_names exist
    if( !all( property_names %in% names(sp_object) ) ){
        stop("A field given in property_names does not exist")
    }
    
    ## mutate the names so that they match `those on the property_names
    nm <- names(sp_object)
    for(ii in names(property_names)){
        nm[nm == property_names[ii]] <- ii
    }
    names(sp_object) <- nm

    ## work out required names
    req_names <- c("name","length","startNode","endNode")
    if(!all(req_names %in% names(sp_object))){
        stop("Not all the required field names are present")
    }

    ## check if SpatialPolygon object - if not buffer using a width
    if(!is_polygon){
        if(!("width" %in% names(property_names))){
            warning("Modifying to spatial polygons using default width")
            sp_object <- rgeos::gBuffer(sp_object, byid=TRUE, width=default_width)
        }else{
            warning("Modifying to spatial polygons using specified width")
            sp_object <- rgeos::gBuffer(sp_object, byid=TRUE, width=sp_object[[ property_names['width'] ]])
        }
    }

    sp_object$area <- raster::area(sp_object)
        
    ## some further basic checks
    sp_object$name <- as.character(sp_object$name)
    sp_object$startNode <- as.character(sp_object$startNode)
    sp_object$endNode <- as.character(sp_object$endNode)
    sp_object$length <- as.numeric(sp_object$length)
    if(!all(is.finite(c(sp_object$length,sp_object$area))) ){
        stop("Some non-finite values of channel lengths or areas found!")
    }
    return( sp_object )
}

