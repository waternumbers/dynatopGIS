#' Simplifies a dynamic TOPMODEL by merging computational bands
#'
#' @description take a model and tries to merge the computational bands to produce a simplier description
#'
#' @param model a RasterBrick as created by create_model
#' @param nband number of bands to merge
#'
#' @return a simplified model
#' @export
band_model <- function(model, cuts){

    ## see if the input mode is valid
    dynatop::check_model(model)

    model <- readRDS("./dynatop/data/Swindale_ordered.rds")
    cuts <- seq(0,2760,by=10)
    hs <- model$hillslope
    hs$class <- hs$split_id
    
    ## create new uids
    grp <- as.numeric( cut(hs$band,cuts) )
    if(any(!is.finite(grp))){
        stop("Cut does not assign each HRU to a new band")
    }
    grp <- 0.5*(hs$class+grp)*(hs$class+grp+1)+grp # unique new uid by cantor pair
    ugrp <- unique(grp)
    hs[,'new_id'] <- NA
    for(ii in 1:length(ugrp)){
        hs[grp==ugrp[ii],'new_id'] <- ii
    }
    hs[,'new_id']  <-  hs[,'new_id'] + max(model$channel$id)
    

    ## create the new table
    hs2 <- data.frame(
        id = unique(hs[,'new_id']),
        area = as.numeric(by(hs$area,hs$new_id,sum)),
        atb_bar = as.numeric(by(hs$atb_bar*hs$area,hs$new_id,sum)),
        s_bar = as.numeric(by(hs$s_bar*hs$area,hs$new_id,sum)),
        delta_x = as.numeric(by(hs$delta_x,hs$new_id,sum)),
        class = as.numeric(by(hs$class,hs$new_id,unique)),
        band = as.numeric(by(hs$delta_x,hs$new_id,min)),
        precip_series="unknown",
        pet_series="unknown",
        qex_max="qex_max_default",
        srz_max="srz_max_default",
        srz_0="srz_0_default",
        ln_t0="ln_t0_default",
        m="m_default",
        td="td_default",
        tex="tex_default",
        stringsAsFactors=FALSE
    )
    hs2$s_bar <- hs2$s_bar/hs2$area
    hs2$atb_bar <- hs2$atb_bar/hs2$area

    ## adapt the redistribution matrices
    ##create summarion matrices on both parts so K %*% A %*% D %*% t(K) ish..
    Dsz <- Matrix(0,


    
    ## find the area, length and d/s id for each split in each band
    new_band <- unique(hs$new_band)
    split <- unique(hs$split_id)
    area <- matrix(NA,length(new_band),length(split))
    length <- matrix(NA,length(new_band),length(split))
    us_id <- matrix(NA,length(new_band),length(split))

    for(ii in 1:length(new_band)){
        for(jj in 1:length(split)){
            idx <- (hs$new_band==new_band[ii]) & (hs$split_id==split[jj])
            if(sum(idx)>0){
                area[ii,jj] <- sum(hs$area[idx])
                length[ii,jj] <- sum(hs$delta_x[idx])
                min_bnd <- min(hs$band[idx])
                jdx <- (hs$band==min_bnd) & (hs$split_id==split[jj])
                ds_id[ii,jj] <- hs$id[jdx]
            }
        }
    }
    
    ## find the split_ids in each band
    split_by_new_band <- by(hs$split_id,hs$new_band,unique,FALSE)

    

    
    for each split_id in each

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

    if(!is(chn,"SpatialPolygonsDataFrame")){
        if(is.character(chn)){
            chn <- rgdal::readOGR(chn,...)
        }else{
            stop("Unknown channel format")
        }
    }

    check_catchment(stck)
    check_channel(chn)





    ## ####################################
    ##  Check required inputs are available
    ## ####################################

    if( !all(hillslope_class %in% names(stck)) ){
        stop("Missing hillslope classification layer")
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
    if( !all( raster::unique(stck[['channel_id']]) %in% chn[['id']] ) ){
        stop("Channels in raster file that are not in the shape file")
    }

    ## check all hillslope id numbers are unique to channels
    uid_hillslope <- raster::unique(stck[[hillslope_class]])
    uid_channel <- raster::unique(stck[['channel_id']])
    if( any( uid_hillslope %in% uid_channel ) ){
        stop("Hillslope class has the same number as the channel id")
    }

    if( !all( c(uid_hillslope,uid_channel) %in% 1:(length(uid_hillslope)+length(uid_channel)) ) &&
        min(c(uid_hillslope,uid_channel)!=1 ) ){
        stop("Hillslope and channel id's chould be sequential numbers starting from 1")
    }

    ## #######################################
    ## Add additional splits by order
    ## #######################################
    cp <- stck[[hillslope_class]]
    x <- stck[['order']]
    cp <- 0.5*(cp+x)*(cp+x+1)+x
    ucp <- unique(cp)
    ucp <- data.frame(ucp,cellStats(stck[['channel_id']],max) + (1:length(ucp)))
    hillslope_id <- raster::subs(cp,ucp) # this is the new HRU class number
    uid_hillslope <- raster::unique(hillslope_id)
    
    ## #######################################
    ## Basic computations of properties
    ## #######################################
    max_index <- max(c(uid_hillslope,uid_channel)) # since Cpp is zero index
    offset <- c(-ncol(stck) + -1:1,-1,1,ncol(stck) + -1:1)
    dx <- rep(sqrt(raster::xres(stck)^2 + raster::yres(stck)^2),8)
    dx[c(2,7)] <- raster::yres(stck)
    dx[c(4,5)] <- raster::xres(stck)
    if( raster::xres(stck) != raster::yres(stck) ){
        warning("Contour length presumes a square grid")
    }
    mres <- (raster::xres(stck) + raster::yres(stck))/2
    cl <- c(rep( mres /(1+sqrt(2)),8),mres) # TO DO this is based on a n octogan - but other papers return a different ratio
    out <- rcpp_hru(raster::getValues(stck[['filled_dem']]),
                    raster::getValues(stck[['gradient']]),
                    raster::getValues(stck[["atanb"]]),
                    raster::getValues(stck[['land_area']]),
                    raster::getValues(stck[['channel_area']]),
                    raster::getValues(stck[['channel_id']])-1,
                    raster::getValues(hillslope_id)-1,
                    offset,
                    dx,
                    cl,
                    max_index)
    ## standardise W by total area
    out$W <- out$W %*% Diagonal(length(out$area),1/out$area)
    
    ## create the model list to be populated
    model <- list()

    ## hillsloep elements
    model$hillslope <- data.frame(
        id = uid_hillslope,
        area = out$area[uid_hillslope],
        atb_bar = out$av_atanb[uid_hillslope]/out$area[uid_hillslope],
        s_bar = out$av_grad[uid_hillslope]/out$area[uid_hillslope],
        delta_x = 1,
        split_id = unname(tapply(raster::getValues(stck[[hillslope_class]]),
                                 raster::getValues(hillslope_id),unique)),
        band = unname(tapply(raster::getValues(stck[['order']]),
                             raster::getValues(hillslope_id),unique)),
        precip_series="unknown",
        pet_series="unknown",
        qex_max="qex_max_default",
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
        precip_series="unknown",
        pet_series="unknown",
        v_ch="v_ch_default",
        stringsAsFactors=FALSE
    )

    ## parameter values
    model$param <- c(qex_max_default=Inf,
                     srz_max_default=0.05,
                     srz_0_default=0.99,
                     ln_t0_default=19,
                     m_default=0.004,
                     td_default=20,
                     tex_default=100,
                     v_ch_default=100)


    ## ########################
    ## Process the redistirbution matrices for the hillslope
    ## ########################
    #browser()
    model$Dex <- model$Dsz <- out$W

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
    

    ## ############################################
    ## Add gauges at all outlets from river network
    ## ############################################
    idx <- which(!(model$channel$id %in% model$channel$next_id))
    model$gauge <- data.frame(
        name = paste("channel",uid_channel[idx],sep="_"),
        id = idx,
        fraction = rep(1,length(idx)),
        stringsAsFactors=FALSE
    )

    ## ##################################
    ## Add point inflow table
    ## ##################################
    ## blank point inflow table
    model$point_inflow <- data.frame(
        name = character(0),
        id = numeric(0),
        fraction = numeric(0),
        stringsAsFactors=FALSE
    )

    ## Save model
    ## saveRDS(model,file=file.path(out_path,paste0(hillslope_class,'.rds')))

    return(model)
}

