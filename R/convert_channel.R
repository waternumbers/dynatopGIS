#' Function for assisting in the conversion of object to be suitable channel inputs to a dynatopGIS object
#'
#' @description Converts SpatialLinesDataFrame or SpatialPolygonsDataFrame to the correct format of SpatialPolygonsDataFrame for dynatopGIS.
#' @param vect_object a SpatVect object or a file which can read by terra::vect to create one
#' @param property_names a named vector of containing the columns of existing data properties required in the final SpatialPolygonsDataFrame
#' @param default_width the width in m to be used for buffering lines to produce polygons
#'
#' @details If the property_names vector contains a width this is used for buffering lines to produce polygons, otherwise the default_width value is used.
#' 
#' @examples
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' vect_lines <- terra::vect(channel_file)
#' property_names <- c(name="identifier",endNode="endNode",startNode="startNode",length="length")
#' chn <- convert_channel(vect_lines,property_names)
#' @export
convert_channel <- function(vect_object,property_names=c(name = "DRN_ID",
                                                       length = "length",
                                                       startNode = "startNode",
                                                       endNode = "endNode",
                                                       width = "width"),
                            default_width=2){
    
    ## read in sp object is a character sting
    if(is.character(vect_object)){
        if(file.exists(vect_object)){
            vect_object <- terra::vect(vect_object)
        }else{
            stop("vect_object is a character string but the file specified does not exist")
        }
    }

    if(!(terra::geomtype(vect_object) %in% c("lines","polygons"))){
        stop("The channel network does not have a line or polygon geometry (even when read in)")
    }
    is_polygon <- terra::geomtype(vect_object)=="polygons"

    ## check all feilds in property_names exist
    if( !all( property_names %in% names(vect_object) ) ){
        stop("A field given in property_names does not exist")
    }
    
    ## mutate the names so that they match those on the property_names
    nm <- names(vect_object)
    for(ii in names(property_names)){
        nm[nm == property_names[ii]] <- ii
    }
    names(vect_object) <- nm

    ## work out required names
    req_names <- c("name","length","startNode","endNode")
    if(!all(req_names %in% names(vect_object))){
        stop("Not all the required field names are present")
    }

    ## check if SpatialPolygon object - if not buffer using a width
    if(!is_polygon){
        if(!("width" %in% names(vect_object))){
            warning("Modifying to spatial polygons using default width")
            vect_object$width <- default_width
            vect_object <- terra::buffer(vect_object, width=default_width/2)
        }else{
            warning("Modifying to spatial polygons using specified width")
            vect_object <- terra::buffer(vect_object, width=vect_object$width/2)
        }
    }
    
    vect_object$area <- terra::expanse(vect_object)

    if(!("width" %in% names(vect_object))){
        warning("Computing width from area and length")
        vect_object$width <- vect_object$area / vect_object$length
    }

    ## some further basic checks
    vect_object$name <- as.character(vect_object$name)
    vect_object$startNode <- as.character(vect_object$startNode)
    vect_object$endNode <- as.character(vect_object$endNode)
    vect_object$length <- as.numeric(vect_object$length)
    if(!all(is.finite(c(vect_object$length,vect_object$width,vect_object$area))) ){
        stop("Some non-finite values of channel lengths, widths or areas found!")
    }
    return( vect_object )
}

