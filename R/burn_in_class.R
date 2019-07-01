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
burn_in_class <- function(project_path,class_name,output_name,burns){
    
    ## check burns are valid
    burns <- as.vector(burns)
    burns <- as.character(burns)

    ## load original classification
    class_file <- file.path(project_path,paste0(class_name,'.tif'))
    if( !file.exists(class_file) ){
        stop("File ",class_file," is missing")
    }
    class <- raster::raster(class_file)

    ## check the files for basing the burns on exist
    burn_files <- file.path(project_path,paste0(burns,'.tif'))
    burn_files_exist <- file.exists(burn_files)
    if(!all(burn_files_exist)){
        stop("Missing files for: ",paste(burns[!burn_files_exist],collapse=" "))
    }
    burn_brck <- raster::brick(list(burn_files))

    ## process if new class
    for(ii in 1:raster::nlayers(burn_brck)){
        class <- overlay(class, burn_brck[[ii]], fun = function(x, y) {
            x[is.finite(y)] <- y[is.finite(y)]
            return(x)
        })
    }

    ## write out raster
    raster::writeRaster(class,file.path(project_path,paste0(output_name,'.tif')))

    ## make json
    tmp <- jsonlite::fromJSON( file.path(project_path,paste0(class_name,'.json')), simplifyVector=FALSE )
    tmp$burns <- burns
    writeLines( jsonlite::toJSON(tmp,pretty=TRUE), file.path(project_path,paste0(output_name,'.json')) )

    return(TRUE)
}
