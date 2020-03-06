#' Create a dynamic TOPMODEL suitabel more model evaluation
#'
#' @description take a classification map and associatd GIS data then computes the GIS properties needed for dynamic TOPMODEL.
#'
#' @param stck a RasterBrick as created by create_brick
#' @param chn a channel object as created by create_channel
#' @param hillslope_class Name of the existing hillslope classification to use to generate the model - a layer in stck
#'
#' @return Logical imdicating it has run. Outputs an rds file named after the classification in the project directory containing the model summary.
#' @export
create_model <- function(stck,chn,hillslope_class){

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

    if( !(hillslope_class %in% names(stck)) ){
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
            warning(paste("To many missing values (",nm,") in channel", ii, "variable - channels and outlets may be wrong"))
        }
    }

    ## check there are no channel bifurcations..since these aren't handled
    if( length(unique(chn[['startNode']])) != nrow(chn) ){
        stop("There are channel bifurcations which are not handled in the simulation code")
    }
    
    ## check all channels in raster and in chn
    if( !all( raster::unique(stck[['channel_id']]) %in% chn[['id']] ) ){
        stop("There are channels in raster file that are not in the shape file")
    }

    ## check all hillslope id numbers are unique to channels
    uid_hillslope <- raster::unique(stck[[hillslope_class]])
    uid_channel <- raster::unique(stck[['channel_id']])
    if( any( uid_hillslope %in% uid_channel ) ){
        stop("Hillslope class id numbers and channel id numbers are not unique")
    }

    ## check that uids are sequential numbers from 1
    uid <- c(uid_hillslope,uid_channel)
    if( !all( uid %in% 1:length(uid) ) &&
        min( uid )!=1 ){
        stop("The hillslope and channel id's chould be sequential numbers starting from 1")
    }

    ## #######################################
    ## Add additional splits by order
    ## #######################################
    cp <- stck[[hillslope_class]]
    x <- stck[['order']]
    cp <- 0.5*(cp+x)*(cp+x+1)+x
    ucp <- unique(cp)
    ucp <- data.frame(ucp,max(uid_channel) + (1:length(ucp)))
    hillslope_id <- raster::subs(cp,ucp) # this is the new HRU class number
    uid_hillslope <- raster::unique(hillslope_id)
    uid <- c(uid_hillslope,uid_channel)
    
    ## #######################################
    ## Basic computations of hillslope properties
    ## #######################################

    ## initialise varaibles
    class <- band <- delta_x <- rep(NA,length(uid))    
    area <- rep(0,length(uid))
    s_bar <- rep(0,length(uid))
    atb_bar <-  rep(0,length(uid))
    W <- Matrix(0,length(uid),length(uid),sparse=TRUE)
    #W <- matrix(0,length(uid),length(uid))
    
    ## computations using raster package
    class[uid_hillslope] <- unname(tapply(raster::getValues(stck[[hillslope_class]]),
                                          raster::getValues(hillslope_id),unique))
    band[uid_hillslope] <- unname(tapply(raster::getValues(stck[["order"]]),
                                         raster::getValues(hillslope_id),unique))
    band[uid_channel] <- max(band,na.rm=TRUE)+1 
    delta_x[uid_hillslope] <- raster::xres(stck)

    ## computations by looping
    dem <- raster::getValues(stck[['filled_dem']])
    ch_id <- raster::getValues(stck[['channel_id']])
    lnd_area <- raster::getValues(stck[['land_area']])
    ch_area <- raster::getValues(stck[['channel_area']])
    hs_id <- raster::getValues(hillslope_id)
    grad <- raster::getValues(stck[['gradient']])
    atb <- raster::getValues(stck[['atanb']])

    #browser()
    ## call cpp version
    ## area atb_bar atb s_bar are altered, W explicitly returned since can't be referenced
    W <- rcpp_hru(dem,hs_id,ch_id,lnd_area,ch_area,grad,atb,W,
                         area,s_bar,atb_bar,
                         c(nrow(stck),ncol(stck),xres(stck)))
    
    W <- W %*% Diagonal(length(area),1/area)
    s_bar <- s_bar/area
    atb_bar <- atb_bar/area
    
    ## create the model list to be populated
    model <- list()

    ## hillslope elements
    model$hillslope <- data.frame(
        id = uid_hillslope,
        area = area[uid_hillslope],
        atb_bar = atb_bar[uid_hillslope],
        s_bar = s_bar[uid_hillslope],
        delta_x = delta_x[uid_hillslope],
        class = class[uid_hillslope],
        sz_band = band[uid_hillslope],
        sf_band = band[uid_hillslope],
        precip="unknown",
        pet="unknown",
        q_sfmax="q_sfmax_default",
        s_rzmax="s_rzmax_default",
        s_rz0="s_rz0_default",
        ln_t0="ln_t0_default",
        m="m_default",
        t_d="t_d_default",
        t_sf="t_sf_default",
        stringsAsFactors=FALSE
    )

    ## channel elements
    model$channel <- data.frame(
        id = uid_channel,
        area = area[uid_channel],
        length=NA,
        sz_band = band[uid_channel],
        sf_band = band[uid_channel],
        next_id=NA,
        precip="unknown",
        pet="unknown",
        v_ch="v_ch_default",
        stringsAsFactors=FALSE
    )

    ## parameter values
    model$param <- c(q_sfmax_default=Inf,
                     s_rzmax_default=0.05,
                     s_rz0_default=0.99,
                     ln_t0_default=19,
                     m_default=0.004,
                     t_d_default=20,
                     t_sf_default=100,
                     v_ch_default=100)


    ## ########################
    ## Process the redistirbution matrices for the hillslope
    ## ########################
    #browser()
    model$Fsf <- model$Fsz <- W

    ## #######################################
    ## Add channel routing information
    ## #######################################
    ##si <- NULL
    ##sj <- NULL
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
            model$channel$next_id[ii] <- chn[['id']][jdx]
        }
    }
    

    ## ############################################
    ## Add gauges at all outlets from river network
    ## ############################################
    idx <- which(is.na(model$channel$next_id))
    model$gauge <- data.frame(
        name = paste("channel",model$channel$id[idx],sep="_"),
        id = model$channel$id[idx],
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

    return(model)
}

