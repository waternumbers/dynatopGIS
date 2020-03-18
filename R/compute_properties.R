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
compute_properties <- function(stck,verbose=FALSE,...){

    check_catchment(ctch)


    ## get properties of required from ctch
    dx <- mean(ctch$raster$res)
    nc <- ctch$raster$dim[2]
    nr <- ctch$raster$dim[1]
    
    ## first pass to compute the number of higher cells and flow directions
    n_higher <- rep(NA,prod(ctch$raster$dim))
    idx <- which(!is.na(dem))
    for(ii in idx){
        jj <- fN(ii,nr,nc)
        is_higher <- (catch$layers$filled_dem[jj] > catch$layers$filled_dem[ii]) &
            !is.finite(catch$layers$channel_id[jj]) # higher cells that aren't channels
        n_higher[ii] <- sum(is_higher,na.rm=TRUE) # na.rm ensure cells on edge of are evaluated
    }

    ## initialise the output
    catch$layers[c("band","gradient","atanb")] <- rep(NA,prod(ctch$raster$dim))
    catch$layers$upslope_area <- catch$layers$land_area
    
    ## work down list of higher cells
    idx <- which(n_higher==0)
    cnt <- 0
    #browser()
    while( length(idx)>0 ){
        cnt <- cnt + 1
        if( verbose ){
            print(paste("Band",cnt))
        }
        
        catch$layers$band[idx] <- cnt
        n_higher[idx] <- -1
        for(ii in idx){
            
            ## skip of no land area
            if( catch$layers$lnd_area[ii] <= 0 ){
                next
            }
            
            if( is.finite(catch$layers$channel_id[ii]) ){
                is_channel <- TRUE
                use_upslope <- TRUE
            }else{
                is_channel <- FALSE
                use_upslope <- FALSE
            }
            
            ## neighbours and grd +ve means lower
            jj <- fN(ii,nr,nc,dx)
            jj$grd <- (catch$layers$filled_dem[ii] - catch$layers$filled_dem[jj$idx])/jj$dx

            ## compute stuff
            if( !use_upslope ){
                ## use downslope gradient
                kk <- which(jj$grd > 0)
                sgn <- 1
                if(length(kk)<1){
                    browser()
                    stop("Cell with no down slope neighbours")
                }
            }else{
                kk <- which(jj$grd < 0)
                sgn <- -1
                if(length(kk)<1){
                    ## Odd case where there is a channel cell with no higher neighbours
                    ## Bugga no higher neighbours use downslope
                    kk <- which(jj$grd > 0)
                    sgn <- 1
                }
            }
            grad_cl <- sgn*jj$grd[kk]*jj$cl[kk]
            ctch$layers$gradient[ii] <- sum(grad_cl) / sum(jj$cl[kk])
            ctch$layers$atanb[ii] <- log( ctch$layers$upslope_area[ii]/sum(grad_cl) )
            if(is.na(ctch$layers$atanb[ii]) | ctch$layers$atanb[ii]==Inf){
                browser()
                warning("None finite topographic index values produced - this requires investigation")
            }
            
            ## cascade downslope if not channel
            if( !is_channel ){
                jdx <- which(jj$grd > 0)
                kk <- jj$idx[jdx]
                grad_cl <- jj$grd[jdx]*jj$cl[jdx]
                ctch$layers$upslope_area[kk] <- ctch$layers$upslope_area[kk] +
                    ctch$layers$upslope_area[ii]*grad_cl/sum(grad_cl)
                
                n_higher[kk] <- n_higher[kk] - 1
            }
        }
                
        # find next set of cells to evaluate
        idx <- which(n_higher==0)
    }
    print(cnt)

    return(ctch)
}
