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

    ## initialise the output
    band <- gradient <- atb <- rep(NA,length(dem))
    upslope_area <- lnd_area
    
    ## work down list of higher cells
    idx <- which(n_higher==0)
    cnt <- 0
    while( length(idx)>0 ){
        cnt <- cnt + 1
        print(cnt)
        band[idx] <- cnt
        n_higher[idx] <- -1
        for(ii in idx){
            if( is.finite(ch_id[ii]) ){
                ## cell contains a channel so no movement downslope
                ## estiamte gradient from upslope cells
                jj <- fN(ii,nr,nc,dx)
                jj$dy <- dem[jj$idx] - dem[ii]
                is_higher <- is.finite(jj$dy) & jj$dy > 0
                jj <- lapply(jj,function(x,y){x[y]},y=is_higher)
                grad_cl <- (jj$dy/jj$dx)*jj$cl # gradient time coutour length
                gradient[ii] <- sum(grad_cl) / sum(jj$cl)
                atb[ii] <- log( upslope_area[ii]/sum(grad_cl) )
            }else{
                ## find lower neighbours and properties
                jj <- fN(ii,nr,nc,dx)
                jj$dy <- dem[ii] - dem[jj$idx]
                is_lower <- is.finite(jj$dy) & jj$dy > 0
                jj <- lapply(jj,function(x,y){x[y]},y= is_lower)
                
                ## work out split
                grad_cl <- (jj$dy/jj$dx)*jj$cl # gradient time coutour length
                
                ## move upslope area down slope
                upslope_area[jj$idx] <- upslope_area[jj$idx] +
                    upslope_area[ii]*grad_cl/sum(grad_cl)
                
                ## compute gradient - weighted sum by contour length
                gradient[ii] <- sum(grad_cl) / sum(jj$cl)
                atb[ii] <- log( upslope_area[ii]/sum(grad_cl) )
                
                ## subtract from from n_higher downslope
                n_higher[jj$idx] <- n_higher[jj$idx] - 1
            }
        }
        # find next set of cells to evaluate
        idx <- which(n_higher==0)
    }

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
