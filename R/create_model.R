#' Create a dynamic TOPMODEL suitabel more model evaluation
#'
#' @description take a classification map and associatd GIS data then computes the GIS properties needed for dynamic TOPMODEL.
#'
#' @param stck a RasterBrick as created by create_brick
#' @param hillslope_class Name of the existing hillslope classification to use to generate the model - a layer in stck
#' @param ... passed to \code{raster::stack} if loading a file
#'
#' @return Logical imdicating it has run. Outputs an rds file named after the classification in the project directory containing the model summary.
#' @export
create_model <- function(ctch,hillslope_class,band_cut=5,verbose=FALSE){

    check_catchment(ctch)
    check_channel(ctch$channel)

    ## ####################################
    ##  Check required inputs are available
    ## ####################################

    if( !(hillslope_class %in% names(ctch$layers)) ){
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
        nm <- sum( is.na(ctch$channel[[ii]]) )
        max_mn <- ifelse(ii=="endNode",1,0)
        if(nm > max_mn){
            warning(paste("To many missing values (",nm,") in channel", ii, "variable - channels and outlets may be wrong"))
        }
    }

    ## check there are no channel bifurcations..since these aren't handled
    if( length(unique(ctch$channel[['startNode']])) != nrow(ctch$channel) ){
        stop("There are channel bifurcations which are not handled in the simulation code")
    }

    ## get unique hillslope class and channel id numbers
    uid <- list(hillslope = unique(ctch$layers[[hillslope_class]]),
                channel = unique(ctch$layers[['channel_id']])
                )
    for(ii in names(uid)){
        uid[[ii]] <- uid[[ii]][!is.na(uid[[ii]])]
    }
    
    ## check all channels in raster are in channel
    if( !all( uid[['channel']] %in% ctch$channel[['id']] ) ){
        stop("There are channels in raster file that are not in the shape file")
    }
    
    ## check that uids are sequential numbers from 1
    for(ii in names(uid)){
        tmp <- range(uid[[ii]])
        if( tmp[1] !=1 | tmp[2] != length(uid[[ii]]) ){
            stop("The",ii,"id numbers should run sequentially from 1")
        }
    }
    

    ## cut up the bands
    cut_band <- as.numeric(cut(ctch$layers$band,breaks=as.numeric(band_cut)))

    ## create hsu ids
    cp <- 0.5*(ctch$layers[[hillslope_class]] + cut_band)*
        (ctch$layers[[hillslope_class]] + cut_band + 1) +
        cut_band
    ucp <- unique(cp[!is.na(cp)])
    ucp <- setNames(1:length(ucp),
                    paste(ucp))
    ## create vector of HSU numbers for the hillslopes
    hillslope_hsu_id <- unname(ucp[paste(cp)]) + max(uid[['channel']]) ## this is the new HRU class number
    uid[['hillslope_hsu']] <- unique(hillslope_hsu_id[!is.na(hillslope_hsu_id)])

    ## Initialise the model
    model <- list()
    
    ## create the model hillslope table and do basic computations
    model$hillslope <- data.frame(
        id = as.integer(uid[['hillslope_hsu']]),
        area = unname(tapply(ctch$layers$land_area,
                             hillslope_hsu_id,
                             sum)
                      [ paste(uid[['hillslope_hsu']]) ]),
        atb_bar = unname(tapply(ctch$layers$land_area*ctch$layers$atanb,
                                hillslope_hsu_id,
                                sum)
                         [ paste(uid[['hillslope_hsu']]) ]),
        s_bar = unname(tapply(ctch$layers$land_area*ctch$layers$gradient,
                              hillslope_hsu_id,
                              sum)
                       [ paste(uid[['hillslope_hsu']]) ]),
        delta_x = unname(tapply(ctch$layers$band,
                                hillslope_hsu_id,
                                function(z){diff(range(z))})
                         [ paste(uid[['hillslope_hsu']]) ]),
        class = unname(tapply(ctch$layers[[hillslope_class]],
                              hillslope_hsu_id,
                              unique)
                       [ paste(uid[['hillslope_hsu']]) ]),
        sz_band = unname(tapply(ctch$layers$band,
                                hillslope_hsu_id,
                                min)
                         [ paste(uid[['hillslope_hsu']])]),
        sf_band = unname(tapply(ctch$layers$band,
                                hillslope_hsu_id,
                                sum)
                         [ paste(uid[['hillslope_hsu']]) ]),
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
    model$hillslope$atb_bar <- model$hillslope$atb_bar / model$hillslope$area
    model$hillslope$s_bar <- model$hillslope$s_bar / model$hillslope$area
    model$hillslope$delta_x <- model$hillslope$delta_x * mean(ctch$raster$res)
    
    ## create the model channel table and do basic computations
    model$channel <- data.frame(
        id = uid[['channel']],
        area = unname(tapply(ctch$layers$channel_area,
                             ctch$layers$channel_id,
                             sum)
                      [ paste(uid[['channel']]) ]),
        length = as.numeric(ctch$channel$length),
        sz_band = max(model$hillslope$sz_band)+1,
        sf_band = max(model$hillslope$sf_band)+1,
        precip="unknown",
        pet="unknown",
        v_ch="v_ch_default",
        stringsAsFactors=FALSE
    )

    ## #############################################
    ## Compute the hillslope HSU flow redistribution
    ## #############################################

    ## initialise the output
    model$hillslope[["sz_flowDir"]] <- rep(list(NULL),nrow(model$hillslope))

    ## work out the band of each hsu
    hsu_bnd <- rep(NA,max(uid$hillslope_hsu))
    hsu_bnd[model$hillslope$id] <- model$hillslope$sz_band
    hsu_bnd[model$channel$id] <- model$channel$sz_band

    ## work out the order in which to pass through the hillslope hsus
    ## from low to high band
    idx <- (1:nrow(model$hillslope))[order(model$hillslope$sz_band)]

    if(verbose){
        verbose_cnt <- c(numeric(2),sum(model$hillslope$area))
    }
    
    ## pass through them computing downslope
    for(ii in idx){
        ## initialise fraction in each direction
        tmp <- rep(0,max(uid$hillslope_hsu))
        id <- model$hillslope$id[ii]

        if(verbose){
            verbose_cnt[1] <- verbose_cnt[1] + 1
            cat("Processing ",verbose_cnt[1]," of ",length(idx)," hillsope HSUs","\n")
        }
        
        bm <- list()
        cnt <- 1
        ## loop through raster cells
        for(jj in which(hillslope_hsu_id == id)){
            ## set to downslope hillslope HSU value
            if( length( ctch$layers$flowDir[[jj]]$idx )== 0 ){
                ds <- ctch$layers$channel_id[jj]
                frc <- 1
            }else{
                ds <- hillslope_hsu_id[ ctch$layers$flowDir[[jj]]$idx ]
                
                ## if this is NA set to downslope channel HSU value
                ch <- ctch$layers$channel_id[ ctch$layers$flowDir[[jj]]$idx ]
                ds[is.na(ds)] <- ch[is.na(ds)]
                frc <- ctch$layers$flowDir[[jj]]$frc
            }
            bm[[cnt]] <- frc
            cnt <- cnt+1
            ## spread area downslope
            ## Need to allow for repeat downslope classes - to see the problem
            ## try running a <- 1:10; b <- c(3,4,5,3); f <- 1:4; a[b] <- a[b]+f
            lmp <- by(frc,ds,sum)
            ds <- as.numeric(names(lmp))
            frc <- as.numeric(lmp)
            
            tmp[ds] <- tmp[ds] + ctch$layers$land_area[jj]*frc


            
        }

        #browser()
        ## Ensure that flow only goes to HSUs later in sequence
        ## these have a higher band number
        kk <- which(tmp>0 & hsu_bnd > model$hillslope$sz_band[ii]) # trim to posistive and lower
        if( length(kk) > 0 ){
            ## if there is one
            ## then the flow can be redistributed to higher bands (more downslope)
            model$hillslope[["sz_flowDir"]][[ii]]$idx <- kk
            model$hillslope[["sz_flowDir"]][[ii]]$frc <- tmp[kk]/sum(tmp[kk])
        }else{
            ## if not look for the HSUs of the same class which are downslope
            kk <- which(model$hillslope$class == model$hillslope$class[ii] &
                        model$hillslope$sz_band > model$hillslope$sz_band[ii])
            if( length(kk) > 0 ){
                ## there are some select the nearest and send all flow there
                kk <- kk[which.min(model$hillslope$sz_band[kk])]
                model$hillslope[["sz_flowDir"]][[ii]]$idx <- model$hillslope$id[kk]
                model$hillslope[["sz_flowDir"]][[ii]]$frc <- 1
            }else{
                ## Cna't send flow anywhere issue a warning
                warning("Some HSU have no downslope output")
            }
        }

        if(verbose){
            #browser()
            verbose_cnt[2] <- verbose_cnt[2] + sum(tmp)
            cat("\t","Processed ",round(100*verbose_cnt[2]/verbose_cnt[3]),"% of area","\n")
        }
    }

    ## set surface flow directions to match saturated zone
    model$hillslope$sf_flowDir <- model$hillslope$sz_flowDir

    ## compute the channel HSU redistributions
    model$channel[["flowDir"]] <- rep(list(NULL),nrow(model$channel))

    for(ii in 1:nrow(model$channel)){
        idx <- which(ctch$channel[['id']]==model$channel$id[ii])
        if(length(idx) != 1){
            stop("Inconsitency between channel raster and shapefile")
        }
        
        ## work out next id
        jdx <- which(ctch$channel[['startNode']]==ctch$channel[['endNode']][idx])
        if( length(jdx) > 0 ){
            model$channel$flowDir[[ii]] <- list(idx=jdx,
                                                frc= 1/length(jdx))
        }
    }
    
    ## parameter values
    model$param <- c(q_sfmax_default=Inf,
                     s_rzmax_default=0.05,
                     s_rz0_default=0.99,
                     ln_t0_default=19,
                     m_default=0.004,
                     t_d_default=20,
                     t_sf_default=100,
                     v_ch_default=100)

    ## ############################################
    ## Add gauges at all outlets from river network
    ## ############################################
    idx <- which(sapply(model$channel$flowDir,length)==0)
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
        id = integer(0),
        fraction = numeric(0),
        stringsAsFactors=FALSE
    )

    return(model)
}

