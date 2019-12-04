#' Compute the topographic index
#'
#' @description Compute the Topographic Index \code{log(area / tanb)}
#' 
#' @param project_path folder which is being used for the analysis
#'
#' @details The algorithm computes the topographic index from the upslope_area and tanb raster files in the project directory.
#' @export
atb <- function(project_path){

    ## load dem and channel id
    file_list <- list(file.path(project_path,'upslope_area.tif'),
                      file.path(project_path,'tanb.tif'))

    brck <- raster::brick(file_list)
    names(brck) <- c('upslope_area','tanb')

    ## check it is projected and square
    if(isLonLat(brck)){
        stop("Currently only works for projected DEM")
    }

    ## compute atb
    atb <- log(brck[['upslope_area']]) - log(brck[['tanb']])

    raster::writeRaster( atb, file.path(project_path,'atb.tif') )

    return(TRUE)
}

    
