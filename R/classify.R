#' Split a catchment into hydrological response units (HRUs)
#'
#' @description Split a catchment into a set hydrological response units (HRUs) according to any number of landscape layers and cuts
#' 
#' @param project_path folder which is being used for the analysis
#' @param output_name name to give the output files
#' @param cuts A named list of cuts of make to form the HRU. Names should correspond to raster layers in the project directory. Values should be numeric and define either the number of bands (single value) or breaks between band (multiple values)
#'
#' @return Logical imdicating it has run. Outputs a raster of classifications and a json of HRU properties to the project directory.
#'
#' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. These are numbered using a cantor pairing scheme. In theory you could reverse out the class of each layer if required but this isn't implimented.
#' @export
split_to_class <- function(project_path,output_name,cuts){

    ## check cuts are valid
    cuts <- as.list(cuts)
    cuts <- lapply(cuts,FUN=function(x){as.numeric(x)})
    cut_names<- names(cuts)
    if(length(cut_names)!=length(unique(cut_names))){
        stop("Replicated names in the cuts list")
    }

    ## check the files for basing the cuts on exist
    cut_files <- file.path(project_path,paste0(cut_names,'.tif'))
    cut_files_exist <- file.exists(cut_files)
    if(!all(cut_files_exist)){
        stop("Missing files for: ",paste(cut_names[!cut_files_exist],collapse=" "))
    }

    
    ## load files for cuts and apply cuts
    cut_brck <- raster::brick(as.list(cut_files))
    names(cut_brck) <- cut_names
    for(ii in cut_names){
        if( length(cuts[[ii]]) > 1 ){
            cut_brck[[ii]] <- raster::cut(cut_brck[[ii]],breaks=cuts[[ii]])
        }else{
            cut_brck[[ii]] <- raster::cut(cut_brck[[ii]],cuts[[ii]])
        }
    }

    ## compute single raster of classes
    cantor_pairing <- function(x){
        if(any(!is.finite(x))){
            return(NA)
        }else{
            
            for(ii in 1:length(x)){
                if(ii==1){
                    o <- x[ii]
                }else{
                    o <- 0.5*(o+x[ii])*(o+x[ii]+1) + x[ii]
                }
            }
            return(o)
        }
    }
    ## calculate the HRU numbers
    hrus <- calc(cut_brck,fun=cantor_pairing)

    ## write out raster
    raster::writeRaster(hrus,file.path(project_path,paste0(output_name,'.tif')))

    ## make and write json
    tmp <- list(cuts=cuts)
    writeLines( jsonlite::toJSON(tmp,pretty=TRUE), file.path(project_path,paste0(output_name,'.json')) )

    return(TRUE)
}

