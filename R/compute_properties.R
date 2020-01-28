#' Function to compute statistics of catchment pixels
#'
#' @description Computes statistics e.g. log(a/tanb) for raster cells
#'
#' @param stck RasterStack or file containing the raster stack (as generated using create_catchment)
#' @param ... additional parameters passed to raster::stack if stck is a file name
#'
#' @details The algorithm works in two passes. The first computes the number of upstream pixels. The second sequences downslope to compute values.
#' @export
compute_properties <- function(stck,...){

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

    ## for testing
    ##rm(list=ls())
    ##stck <- raster::stack("./Swindale_stck.gri")
    ##source("./dynatopGIS/R/fN.R")
    

    ## extract dem and channel_id values and dem resolution
    dem <- raster::getValues(stck[['filled_dem']])
    ch_id <- raster::getValues(stck[['channel_id']])
    lnd_area <- raster::getValues(stck[['land_area']])

    ## get properties of required from stck
    dx <-  raster::xres(stck)
    nc <- ncol(stck)
    nr <- nrow(stck)
    
    ## first pass to compute the number of higher cells
    n_higher <- rep(NA,length(dem))
    idx <- which(!is.na(dem))
    for(ii in idx){
        jj <- fN(ii,nr,nc)
        is_higher <- dem[jj] > dem[ii] & !is.finite(ch_id[jj]) # higher cells that aren't channels
        n_higher[ii] <- sum(is_higher,na.rm=TRUE) # na.rm ensure cells on edge of are evaluated
    }

    #browser()
    ## initialise the output
    band <- gradient <- atb <- rep(NA,length(dem))
    upslope_area <- lnd_area
    
    ## work down list of higher cells
    idx <- which(n_higher==0)
    cnt <- 0
    
    while( length(idx)>0 ){
        cnt <- cnt + 1
        
        band[idx] <- cnt
        n_higher[idx] <- -1
        for(ii in idx){
            
            ## skip of no land area
            if( lnd_area[ii] <= 0 ){
                next
            }
            
            if( is.finite(ch_id[ii]) ){
                is_channel <- TRUE
                use_upslope <- TRUE
            }else{
                is_channel <- FALSE
                use_upslope <- FALSE
            }
            
            ## neighbours and grd +ve means lower
            jj <- fN(ii,nr,nc,dx)
            jj$grd <- (dem[ii] - dem[jj$idx])/jj$dx

            ## compute stuff
            if( !use_upslope ){
                kk <- which(jj$grd > 0)
                sgn <- 1
                if(length(kk)<1){
                    if( any(is.na(jj$grd)) ){
                        ## Bugga no lower neighbours use upslope
                        kk <- which(jj$grd < 0)
                        sgn <- -1
                    }else{
                        stop("None edge cell with no down slope neighbours")
                    }
                    
                }
            }else{
                kk <- which(jj$grd < 0)
                sgn <- -1
                if(length(kk)<1){
                    ## Bugga no higher neighbours use downslope
                    kk <- which(jj$grd > 0)
                    sgn <- 1
                }
            }
            grad_cl <- sgn*jj$grd[kk]*jj$cl[kk]
            gradient[ii] <- sum(grad_cl) / sum(jj$cl[kk])
            atb[ii] <- log( upslope_area[ii]/sum(grad_cl) )
            if(atb[ii]==Inf){
                warning("None finite topographic index values produced - this requires investigation")
            }
            
            ## cascade downslope if not channel
            if( !is_channel ){
                kk <- jj$idx[which(jj$grd > 0)]
                
                upslope_area[kk] <- upslope_area[kk] +
                    upslope_area[ii]*grad_cl/sum(grad_cl)
                
                n_higher[kk] <- n_higher[kk] - 1
            }
        }
                
        # find next set of cells to evaluate
        idx <- which(n_higher==0)
    }
    print(cnt)
    #browser()
    ## write back to stack
    stck <- raster::setValues(stck, band, layer=which(names(stck)=="order"))
    stck <- raster::setValues(stck, gradient, layer=which(names(stck)=="gradient"))
    stck <- raster::setValues(stck, upslope_area, layer=which(names(stck)=="upslope_area"))
    stck <- raster::setValues(stck, atb, layer=which(names(stck)=="atanb"))
 
    
    
    if(length(stck_file)>0){
        raster::writeRaster(stck,stck_file)
        return(stck_file)
    }else{
        return(stck)
    }

}
