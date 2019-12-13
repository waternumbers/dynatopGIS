#' Create a dynamic TOPMODEL suitabel more model evaluation
#'
#' @description take a classification map and associatd GIS data then computes the GIS properties needed for dynamic TOPMODEL.
#' 
#' @param project_path folder which is being used for the analysis
#' @param hillslope_class Name of the existing hillslope classification to use to generate the model
#'
#' @return Logical imdicating it has run. Outputs an rds file named after the classification in the project directory containing the model summary.
#' @export
create_model <- function(brck,chn,hillslope_class,out_path=getwd()){

    ## ####################################
    ##  Check required inputs are available
    ## ####################################
    
    ## check we have all the raster files we need
    layer_names <- c('dem','land_area','channel_area','channel_id',hillslope_class)
    if( !all(layer_names %in% names(brck)) ){
        stop(paste( "Missing layers:",
                   paste(setdiff(layer_names,names(brck)),collapse=" ")))
    }

    
    ## ###############################################################
    ## Sanity checks that :
    ##   - Channel shape and raster file are consistent
    ##   - Shape file contains all information
    ##   - No river connectivity situations that aren't handled occur
    ## ###############################################################
    
    ## check channel shape file values are not missing
    for(ii in c('endNode','startNode','length')){
        nm <- sum( is.na(chn[[ii]]) )
        max_mn <- ifelse(ii=="endNode",1,0)
        if(nm > max_mn){
            warning(paste("To many missing values (",nm,") in channel", ii, "variable - channels and gauges may be wrong"))
        }
    }
    
    ## check there are no channel bifurcations..since these aren't handled
    if( length(unique(chn[['startNode']])) != nrow(chn) ){
        stop("There are channel bifurcations which are not handled in the code")
    }
    ## check all channels in raster and in chn
    if( !all( raster::unique(brck[['channel_id']]) %in% chn[['id']] ) ){
        stop("Channels in raster file that are not in the shape file")
    }

    ## check all hillslope id numbers are unique to channels
    uid_hillslope <- raster::unique(brck[[hillslope_class]])
    uid_channel <- raster::unique(brck[['channel_id']])
    if( any( uid_hillslope %in% uid_channel ) ){
        stop("Hillslope class has the same number as the channel id")
    }

    if( !all( c(uid_hillslope,uid_channel) %in% 1:(length(uid_hillslope)+length(uid_channel)) ) &&
        min(c(uid_hillslope,uid_channel)!=1 ) ){
        stop("Hillslope and channel id's chould be sequential numbers starting from 1")
    }
    
    ## #######################################
    ## Basic computations of properties 
    ## #######################################
    max_index <- max(c(uid_hillslope,uid_channel)) # since Cpp is zero index
    offset <- c(-ncol(brck) + -1:1,-1,1,ncol(brck) + -1:1)
    dx <- rep(sqrt(raster::xres(brck)^2 + raster::yres(brck)^2),8)
    dx[c(2,7)] <- raster::yres(brck)
    dx[c(4,5)] <- raster::xres(brck)
    if( raster::xres(brck) != raster::yres(brck) ){
        warning("Contour length presumes a square grid")
    }
    mres <- (raster::xres(brck) + raster::yres(brck))/2
    cl <- c(rep( mres /(1+sqrt(2)),8),mres) # TO DO this is based on a n octogan - but other papers return a different ratio
    out <- fun_hru(as.vector(brck[['dem']]),
                   as.vector(brck[['gradient']]),
                   as.vector(brck[['land_area']]),
                   as.vector(brck[['channel_area']]),
                   as.vector(brck[['channel_id']])-1,
                   as.vector(brck[[hillslope_class]])-1,
                   offset,
                   dx,
                   cl,
                   max_index)
    ## standardise W by total area
    for(ii in 1:ncol(out$W)){
        out$W[,ii] <- out$W[,ii]/out$area[ii]
    }

    
    ## create the model list to be populated
    model <- list()

    ## hillsloep elements
    model$hillslope <- data.frame(
        id = uid_hillslope,
        area = out$area[uid_hillslope],
        s_bar = out$av_grad[uid_hillslope]/out$area[uid_hillslope],
        precip_input="unknown",
        pet_input="unknown",
        srz_max="srz_max_default",
        srz_0="srz_0_default",
        ln_t0="ln_t0_default",
        m="m_default",
        td="td_default",
        tex="tex_default",
        stringsAsFactors=FALSE
    )

    ## channel elements
    model$channel <- data.frame(
        id = uid_channel,
        area = out$area[uid_channel],
        length=NA,
        next_id=NA,                
        precip_input="unknown",
        pet_input="unknown",
        v_ch="v_ch_default",
        stringsAsFactors=FALSE
    )

    ## parameter values
    model$param <- c(srz_max_default=0.05,
                     srz_0_default=0.99,
                     ln_t0_default=19,
                     m_default=0.004,
                     td_default=20,
                     tex_default=100,
                     v_ch_default=100)

    
    ## ########################
    ## Process the redistirbution matrices for the hillslope
    ## ########################
    model$Dex <- Matrix::Matrix(out$W)
    model$Wsz <- Matrix::Matrix(out$W[uid_hillslope,uid_hillslope])
    model$Fsz <- Matrix::Matrix(out$W[uid_channel,uid_hillslope])

    ## #######################################
    ## Add channel routing information
    ## #######################################
    si <- NULL
    sj <- NULL
    for(ii in uid_channel){
        idx <- which(chn[['id']]==ii)
        if(length(idx) != 1){
            stop("Inconsitency between channel raster and shapefile")
        }
        ## set length
        model$channel$length[ii] <- chn[['length']][idx]
        ## work out next id
        jdx <- which(chn[['startNode']]==chn[['endNode']][idx])
        if( length(jdx) ==1 ){ # case of jdx>1 handled in bifurcation check
            si <- c(si,which(uid_channel==chn[['id']][jdx]))
            sj <- c(sj,ii)
        }        
    }
    model$Wch <- Matrix::sparseMatrix(i=si,j=sj,x=1)

    ## ############################################
    ## Add gauges at all outlets from river network
    ## ############################################
    idx <- which(colSums(model$Wch)==0)
    model$gauge <- data.frame(
        name = paste("channel",uid_channel[idx],sep="_"),
        channel_id = idx,
        fraction = rep(1,length(idx)),
        stringsAsFactors=FALSE
    )

    ## ##################################
    ## Add point inflow table
    ## ##################################
    ## blank point inflow table
    model$point_inflow <- data.frame(
        name = character(0),
        channel_id = numeric(0),
        fraction = numeric(0),
        stringsAsFactors=FALSE
    )

    ## Save model
    ## saveRDS(model,file=file.path(out_path,paste0(hillslope_class,'.rds')))

    return(model)
}
        
