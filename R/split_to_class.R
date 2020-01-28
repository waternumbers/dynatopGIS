#' Split a catchment into hydrological response units (HRUs)
#'
#' @description Split a catchment into a set hydrological response units (HRUs) according to any number of landscape layers and cuts
#'
#' @param stck a RasterStack as created by create_catchment
#' @param split_name name to give the split
#' @param cuts A named list of cuts of make to form the HRU. Names should correspond to raster layers in the project directory. Values should be numeric and define either the number of bands (single value) or breaks between band (multiple values)
#' @param json_path the folder to write out json file containing the summary of splits to
#'
#' @return as for stck but containing an additional layer with the split classes
#'
#' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. These are numbered using a cantor pairing scheme. In theory you could reverse out the class of each layer if required but this isn't implimented.
#' @export
split_to_class <- function(stck,split_name,cuts,json_path=getwd()){

        if(!("RasterStack" %in% class(stck))){
        if( is.character(stck) ){
            stck_file <- stck
            stck <- raster::stack(stck,...)
        }else{
            stop("Unknown format for input")
        }
    }else{
        stck_file <- character(0)
    }
    check_catchment(stck)

    ## check cuts are valid
    
    cuts <- as.list(cuts)
    cuts <- lapply(cuts,FUN=function(x){as.numeric(x)})
    cut_names<- names(cuts)
    if( length(cut_names) != length(cuts) ){
        stop("Unnamed cuts")
    }
    if( length(cut_names)!=length(unique(cut_names)) ){
        stop("Replicated names in the cuts list")
    }

    ## check output name is valid
    if( split_name %in% names(stck) ){
        stop("Output name already used")
    }

    ## check the layers exist
    if( !all(cut_names %in% names(stck)) ){
        stop(paste( "Missing layers:",
                   paste(setdiff(cut_names,names(stck)),collapse=" ")))
    }

    ## work out new cuts by cantor_pairing
    #browser()
    for(ii in 1:length(cuts)){
        x <- raster::cut(stck[[cut_names[ii]]],breaks=as.numeric(cuts[[ii]]))
        if(ii == 1){
            cp <- x
        }else{
            cp <- 0.5*(cp+x)*(cp+x+1)+x
        }
    }

    ## renumber from cantour paring to get sequential values
    ucp <- unique(cp)
    ucp <- data.frame(ucp,cellStats(stck[['channel_id']],max) + (1:length(ucp)))
    cp <- raster::subs(cp,ucp)

    ## add back into brick
    cp <- raster::mask(cp,stck[['dem']])
    names(cp) <- split_name
    stck <- addLayer(stck,cp)

    ## make and write json
    tmp <- list(cuts=cuts)
    writeLines( jsonlite::toJSON(tmp,pretty=TRUE),
               file.path(json_path,paste0(split_name,'.json')) )

    if(length(stck_file)>0){
        writeRaster(stck,stck_file)
        return(stck_file)
    }else{
        return(stck)
    }

}

