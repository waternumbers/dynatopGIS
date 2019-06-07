#' Load a dem and compute its area
#'
#' @description Reads the DEM file and writes it to the project directory
#'
#' @param dem_file name of file to load as the dem
#' @param project_path the path of the project to which to write dem file
#'
#' @return TRUE if successful
#'
#' @details Currently reads and then writes file using the raster package. No checks as to the projection are enforced. Output is project_path/dem.tif
#' @export
load_dem <- function(dem_file,project_path='.'){
    out_file <- file.path(project_path,'dem.tif')

    dem <- raster::raster(dem_file)
    writeRaster(dem,out_file,format='GTiff')
    return(TRUE)
}
