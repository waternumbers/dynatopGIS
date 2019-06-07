#' Discrete a catchment into hydrological response units (HRUs)
#'
#' @description Discretise a catchment into a set hydrological response units (HRUs) according to any number of landscape layers and cuts
#' 
#' @param project_dir Directory of the project
#' @param class_name Name fo the existing classification onto which to burn the new classes
#' @param output_name name to give the output files
#' @param cuts A named list of cuts of make to form the HRU. Names should correspond to raster layers in the project directory. Values should be numeric and define either the number of bands (single value) or breaks between band (multiple values)
#' @param burns list  Named list of geometries (supplied as rasters) to burn into discretisation of HRUs that will be stamped onto the classification. Overrides any classification already defined
#'
#' @return Logical imdicating it has run. Outputs a raster of classifications and a json of HRU properties to the project directory.
#'
#' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. These are numbered using a cantor pairing scheme. In theory you could reverse out the class of each layer if required but this isn't implimented.
#' @export
#'
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
    hrus <- calc(cut_brck,fun=cantor_pairing)
    raster::writeRaster(hrus,file.path(project_path,paste0(output_name,'.tif')))
}

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
    for(ii in 1:nlayer(burn_brck)){
        class <- overlay(class, burn_brck[[ii]], fun = function(x, y) {
            x[!is.na(y[])] <- y[]
            return(x)
        })
    }

    raster::writeRaster(hrus,file.path(project_path,paste0(output_name,'.tif')))

    return(TRUE)
}
