#' Function to compute statistics of catchment pixels
#'
#' @description Computes statistics e.g. log(a/tanb) for raster cells
#'
#' @param stck RasterStack or file containing the raster stack (as generated using create_catchment)
#' @param verbose print out additional diagnostic information
#' @param ... additional parameters passed to raster::stack if stck is a file name
#'
#' @details The algorithm works in two passes. The first computes the number of upstream pixels. The second sequences downslope to compute values.
#' @export
compute_properties <- function(ctch,verbose=FALSE){

    check_catchment(ctch,c("filled_dem","channel_id",
                           "land_area","channel_area"))


    ## get properties of required from ctch
    dx <- mean(ctch$raster$res)
    nc <- ctch$raster$dim[2]
    nr <- ctch$raster$dim[1]
    #browser()
    ## First pass to compute the number of higher cells, flow directions
    ## and variables that are not order specific

    cat("Pass 1")
    
    n_higher <- rep(NA,prod(ctch$raster$dim)) ## number of higher cells
    grad_cl <- rep(NA,prod(ctch$raster$dim)) ## sum of gradient * contour length of flow paths
    
    ctch$layers[["flowDir"]] <- rep(list(NULL),prod(ctch$raster$dim))
    ctch$layers[["gradient"]] <- rep(NA,prod(ctch$raster$dim))
    
    idx <- which(!is.na(ctch$layers$filled_dem))
    for(ii in idx){

        ## skip of no land area
        if( ctch$layers$land_area[ii] <= 0 ){
            next
        }

        ## if cell has a channel_id then set flags
        if( is.finite(ctch$layers$channel_id[ii]) ){
            is_channel <- TRUE
        }else{
            is_channel <- FALSE
        }

        
        ## work out neighbours and gradient to them
        jj <- fN(ii,nr,nc,dx)
        jj$grd <- (ctch$layers$filled_dem[ii] - ctch$layers$filled_dem[jj$idx])/jj$dx

        
        ## work out higher cells (negative gradient and not channels)
        ## these need evaluating first
        is_higher <- which( (jj$grd < 0) & !is.finite(ctch$layers$channel_id[jj$idx]) )

        ## work out the lower cells (positive gradient)
        is_lower <- which(jj$grd > 0)

        ## number fo higher cells
        n_higher[ii] <- length(is_higher) # na.rm ensure cells on edge are evaluated

        ## work out the flow directions - based on lower cells
        if( is_channel ){
            ## Cells with a channel element
            if( n_higher[ii] > 0 ){
                ## use upslope cells for gradient if possible
                g <- -jj$grd[is_higher]
                cl <- jj$cl[is_higher]
            }else{
                ## else use lower cells with warning
                warning(paste("Cell",ii,"is a channel cell with no higher neighbours"))
                g <- jj$grd[is_lower]
                cl <- jj$cl[is_lower]
            }
            gcl <- g*cl
            grad_cl[ii] <- sum(gcl)
            ctch$layers$gradient[ii] <- grad_cl[ii] / sum(cl)
        }else{
            ## cells with have no channel element
            if( length(is_lower) > 0 ){
                ## use downslope cells for gradient if possible
                g <- jj$grd[is_lower]
                cl <- jj$cl[is_lower]
            }else{
                ## else use higher cells with warning
                warning(paste("Cell",ii,"is a hillslope cell with no lower neighbours"))
                g <- jj$grd[is_higher]
                cl <- jj$cl[is_higher]
            }
            gcl <- g*cl
            grad_cl[ii] <- sum(gcl)
            ctch$layers$gradient[ii] <- grad_cl[ii] / sum(cl)

            ## do flow direction calculation
            if( length(is_lower) > 0 ){
                ctch$layers$flowDir[[ii]] <- list(idx=jj$idx[is_lower],
                                                  frc = gcl/grad_cl[ii])
            }
        }
    }

    
    ## Second pass to compute ordered valraibles
    cat("Pass 2")
    
    ctch$layers[["band"]] <- ctch$layers[["atanb"]] <- rep(NA,prod(ctch$raster$dim))
    ctch$layers$upslope_area <- ctch$layers$land_area
Q
    
    ## work down list of higher cells
    idx <- which(n_higher==0)
    cnt <- 0
    while( length(idx)>0 ){
        cnt <- cnt + 1
        if( verbose ){
            print(paste("Band",cnt))
        }
        
        ctch$layers$band[idx] <- cnt
        n_higher[idx] <- -1
        for(ii in idx){

            ## skip if no downslope linkages
            if(length(ctch$layers$flowDir[[ii]]) == 0){
                next
            }

            ## move upslopre area downslope
            kk <- ctch$layers$flowDir[[ii]]$idx
            ff <- ctch$layers$flowDir[[ii]]$frc
            ctch$layers$upslope_area[kk] <- ctch$layers$upslope_area[kk] +
                ctch$layers$upslope_area[ii]*ff
            n_higher[kk] <- n_higher[kk] - 1
        }
        
        # find next set of cells to evaluate
        idx <- which(n_higher==0)
    }

    ## compute tographic from summaries
    browser()
    ctch$layers$atanb <- log( ctch$layers$upslope_area / grad_cl )
    if( any(ctch$layers$atanb==Inf) ){       
        warning("None finite topographic index values produced - this requires investigation")
    }

    return(ctch)
}
