#' Function for assisting in the conversion of object to be suitable channel inputs to a dynatopGIS object
#'
#' @description Converts SpatialLinesDataFrame or SpatialPolygonsDataFrame to the correct format of SpatialPolygonsDataFrame for dynatopGIS.
#' @param sp_object a SpatialLinesDataFrame or SpatialPolygonsDataFrame object or a file which can read by raster::shapefile to create one
#' @param property_names a named vector of containing the columns of existing data properties required in the final SpatialPolygonsDataFrame
#' @param default_width the width in m to be used for buffering lines to produce polygons
#'
#' @details If the property_names vector contains a width this is used for buffering lines to produce polygons, otherwise the default_width value is used.
#' 
#' @examples
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' sp_lines <- rgdal::readOGR(channel_file)
#' property_names <- c(name="identifier",endNode="endNode",startNode="startNode",length="length")
#' chn <- convert_channel(sp_lines,property_names)
#' @export
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

