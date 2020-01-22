#' Fill the sinks with a dem
#'
#' @description TO DO
#'
#' @param stck a rasterBrick object or whatever
#' @param ... additional varaibles passed to raster::brick if stck is a file name
#'
#' @return a rasterBrick object
#'
#' @details TO DO
#'
#' @export
sink_fill <- function(stck,...){

    if(!("RasterStack" %in% class(stck))){
        if( is.character(stck) ){
            stck_file <- stck
            stck <- raster::stack(stck,...)
        }else{
            stop("Unknown stck format")
        }
    }else{
        stck_file <- character(0)
    }
    
    offset <- c(-ncol(stck) + -1:1,-1,1,ncol(stck) + -1:1)
    ##browser()

    ## raster are indexed column dominated i.e. fill lef to right first not top to bottom
    ## this function is tested against raster::adjacent since this matches the
    ## extraction order of raster::getValues
    fN <- function(k,nr,nc){
        ## find i,j location
        j <- ((k-1)%/%nc) +1
        i <- k - (j-1)*nc
        m <- matrix(c(i-1,i,i+1,i-1,i+1,i-1,i,i+1,
                      j-1,j-1,j-1,j,j,j+1,j+1,j+1),8,2)
        m <- m[ m[,1]>0 & m[,1]<=nr & m[,2]>0 & m[,2]<=nc ,]
        
        return((m[,2]-1)*nc + m[,1])
    }
    
    ## work out lowest neighbour
    nc <- ncol(stck)
    nr <- nrow(stck)
    min_neighbour <- rep(NA,n)
    for(ii in 1:length(dem)){
        jj <- ii + offset
        min_neighbour[ii] <- min(dem[fN(ii,nr,nc)])
    }

    ## determine sinks and set to Inf - this loop should only increase dem value
    idx <- which(min_neighbour >= dem)
    dem[idx] <- Inf
    while(length(idx)>0){
        io <- idx[1]+offset
        io <- io[io>0 & io<=n]
        for(ii in io){
            if(is.finite(dem[ii])){
                jj <- ii + offset
                min_neighbour[ii] <- min(dem[jj[jj>0 & jj<n]])
                if(!is.na(min_neighbour[ii]) && (min_neighbour[ii] >= dem[ii])){
                    idx <- c(idx,ii)
                    dem[ii] <- Inf
                }
            }
        }
        idx <- idx[-1]
    }

    ## populate starting with lowest
    idx <- which(min_neighbour >= dem)
    while(length(idx)>0){
        print(idx)
        mi <- which.min(min_neighbour[idx])
        ii <- idx[mi] # current value to valuate
        jj <- ii + offset
        dd <- dem[jj[jj>0 & jj<n]]
        dem[ii] <- mean(dd[is.finite(dd)])
        
        for(kk in jj){
            jjj <- kk + offset
            min_neighbour[kk] <- min(dem[jjj[jjj>0 & jjj<n]])
        }
        idx <- which(min_neighbour >= dem)
        #idx <- idx[-mi]
    }
    

                           
    out <- rcpp_sink_fill(raster::getValues(stck[["dem"]]),
                          raster::getValues(stck[["channel_id"]]),
                          offset)
    #browser()
    stck <- raster::setValues(stck, out, layer=which(names(stck)=="filled_dem"))
    
    #warning("Sink filling is not implimented...this just copies the dem")
    #stck[['filled_dem']] <- stck[['dem']]

    if(length(stck_file)>0){
        raster::writeRaster(stck,stck_file)
        return(stck_file)
    }else{
        return(stck)
    }

}
