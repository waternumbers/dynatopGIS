#' Burn in distinct classes of hydrological response units (HRUs)
#'
#' @description Burn in distinct regions to a classification map
#' 
#' @param project_path folder which is being used for the analysis
#' @param class_name Name of the existing classification onto which to burn the new classes
#' @param output_name name to give the output files
#' @param burns list  Named list of geometries (supplied as rasters) to burn into discretisation of HRUs that will be stamped onto the classification. Overrides any classification already defined
#'
#' @return Logical indicating it has run. Outputs a raster of classifications and a json of HRU properties to the project directory.
#'
#' @details Replaces exisiting classes with those from the specified classification maps to be burnt in. burns are applied in the order theya are input.
#' 
#' @export
burn_in_class <- function(brck,class_name,output_name,burns){
    
    ## check burns are valid
    burns <- as.vector(burns)
    burns <- as.character(burns)

    ## check layers exist
    if( !all(c(class_name,burns) %in% names(brck)) ){
        stop("Missing layers used to perform burn")
    }

    ## check output name isn;t already in use
    if( output_name %in% names(brck) ){
        stop("Output name already used")
    }

    ## process in new class
    nc <- brck[[class_name]]
    for(ii in burns){
        idx <- Which(is.finite(brck[[ii]]),cells=TRUE)
        nc[idx] <- brck[[ii]][idx]
    }

    nc <- raster::mask(nc,brck[['dem']])
    brck[[output_name]] <- nc
    
    ## make json
    tmp <- jsonlite::fromJSON( file.path(json_path,paste0(class_name,'.json')), simplifyVector=FALSE )
    tmp$burns <- burns
    writeLines( jsonlite::toJSON(tmp,pretty=TRUE), file.path(json_path,paste0(output_name,'.json')) )
    
    return(brck)
}
