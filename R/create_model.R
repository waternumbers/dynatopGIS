#' Create a dynamic TOPMODEL suitabel more model evaluation
#'
#' @description take a classification map and associatd GIS data then computes the GIS properties needed for dynamic TOPMODEL
#' 
#' @param project_path directory of the project being worked on
#' @param hillslope_class Name of the existing hillslope classification to use to generate the model
#'
#' @return Logical imdicating it has run. Outputs an rds file named after the classification in the project directory containing the model summary.
#' @export
create_model <- function(project_path,hillslope_class){

    ## check we have all the files we need
    layer_names <- c('dem','land_area','atb','channel_area','channel_id',hillslope_class)
    file_names <- file.path( project_path, paste0(layer_names,'.tif') )
    if( !all(file.exists(file_names)) ){
        stop("Missing files for: ",paste(layer_names[!file.exists(file_names)],collapse=" "))
    }
    
    ## load files
    brck <- raster::brick(as.list(file_names))
    ## alter names so valid R names which appear in brck
    names(brck) <- make.names(layer_names)
    hillslope_class <- make.names(hillslope_class)

    ## initialise the hillslope
    model <- list(hillslope = data.frame(
                      id = raster::unique(brck[[hillslope_class]]),
                      area = raster::zonal(brck[['land_area']],brck[[hillslope_class]],
                                           'sum',digits=4)[,2],
                      atb_bar = raster::zonal(brck[['atb']],brck[[hillslope_class]],
                                              'mean',digits=4)[,2],
                      precip_input="unknown",
                      pet_input="unknown",
                      srz_max="srz_max_default",
                      srz_0="srz_0_default",
                      ln_t0="ln_t0_default",
                      m="m_default",
                      td="td_default",
                      tex="tex_default"
                  ),
                  channel = data.frame(
                      id = raster::unique(brck[['channel_id']]),
                      area = raster::zonal(brck[['channel_area']],brck[['channel_id']],
                                           'sum',digits=4)[,2],
                      precip_input="unknown",
                      pet_input="unknown"
                  ),
                  param <- c(srz_max_default=0.05,
                             srz_0_default=0.99,
                             ln_t0_default=19,
                             m_default=0.004,
                             td_default=20,
                             tex_default=100)
                  )

    ## make distance for computing the weighting matrices
    dist <- matrix(sqrt(raster::xres(brck)^2 + raster::yres(brck)^2),3,3)
    dist[2,1] <- dist[2,3] <- raster::xres(brck)
    dist[1,2] <- dist[3,2] <- raster::yres(brck)
    dist[2,2] <- 0
    ## create sequential hillslope and channel numbers to avoid messing about in the cpp code
    seq_hillslope <- raster::subs(brck[[hillslope_class]],
                                  data.frame(model$hillslope$id,(1:length(model$hillslope$id))-1))
    seq_channel <- raster::subs(brck[['channel_id']],
                                data.frame(model$channel$id,(1:length(model$channel$id))-1))
    tmp <- fun_redistribution(as.matrix(brck[['dem']]),
                              as.matrix(brck[['land_area']]),
                              as.matrix(seq_hillslope),
                              as.matrix(seq_channel),
                              length(model$hillslope$id),
                              length(model$channel$id),
                              dist)
    
    ## name the matrices 
    rownames(tmp$hillslope) <- colnames(tmp$hillslope) <- model$hillslope$id
    colnames(tmp$channel) <- model$hillslope$id
    rownames(tmp$channel) <- model$channel$id
    ## standardise by total area
    for(ii in 1:ncol(tmp$hillslope)){
        tmp$hillslope[,ii] <- tmp$hillslope[,ii] / model$hillslope[ii,'area']
        tmp$channel[,ii] <- tmp$channel[,ii] / model$hillslope[ii,'area']
    }
    
    model$Wsat <- model$Wex <- tmp$hillslope
    model$Fsat <- model$Fex <- tmp$channel

    saveRDS(model,file=file.path(project_path,paste0(hillslope_class,'.rds')))

    return(TRUE)
}
        
