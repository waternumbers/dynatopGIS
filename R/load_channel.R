#' load a chanel file and ensure it has the correct properties
#'
#' @description Import channel data from an OGR file to the project directory
#'
#' @param file_name name of file to read in
#' @param property_names named vector of columns of the spatial data frame to use for channel properties
#' @param project_path folder which is being used for the analysis
#' @param default_width default width of a channel if not specified in property_names. Defaults to 2 metres.
#' @return TRUE if runs otherwise an error message.
#'
#' @details reading in the input file and covnerts it then saves a shape file of polygons with properties chanel_id,length,width,startNode,endNode in the project directory under the name \code{channel.shp}. Feild in the input file corresponding to the properties can be specified in the \code{property_names} vector.
#' 
#' @export
load_channel <- function(file_name,property_names,project_path='.',default_width=2){
    ## check property names
    if(!is.vector(property_names) ||
       length(names(property_names))!=length(property_names) ){
        stop("property_names of wrong type")
    }
    if( !("chanel_id" %in% names(property_names)) ){
        stop("Must specify an id feild in property_names")
    }
    if( ("id" %in% names(property_names)) ){
        stop("The name is reserved")
    }
    if( !("width" %in% names(property_names)) ){
        warning("Using default width of ",default_width,"m")
    }
    
    

    ## load and check chanel data
    chanel <- rgdal::readOGR( file_name ,stringsAsFactors=FALSE)
    if(class(chanel)!="SpatialLinesDataFrame"){
        stop("Chanel data is not a spatial lines data frame object")
    }

    ## see if all names exist
    if(!all( property_names %in% names(chanel))){
        stop("Not all feilds names in property_names are in the chanel data frame")
    }

    ## alter chanel
    chanel <- chanel[,property_names]
    names(chanel) <- names(property_names)

    ## add a unique id for each unique chanel_id
    cid <- as.character( chanel[['chanel_id']] )
    ucid <- unique(cid)
    id <- 1:length(ucid); names(id) <- ucid
    chanel[['id']] <- id[cid]
    
    ## add blank values to missing feilds
    missing <- setdiff( c('startNode','endNode','length'), names(chanel) )
    if(length(missing)>0){
        warning("The following feilds are missing and blank values will be added: ",
                paste(missing,collapse=" "))
        for(ii in missing){ chanel[[ii]] <- NA }
    }

    ## add default width
    if( !("width" %in% names(chanel)) ){
        chanel[['width']] <- default_width
    }

    ## buffer to convert to spatial polygons object
    buffered_chanel <- rgeos::gBuffer(chanel, byid=TRUE, width=chanel[['width']])
    raster::shapefile(buffered_chanel, file.path(project_path,'channel'))
    return(TRUE)
}
