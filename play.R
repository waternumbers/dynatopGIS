rm(list=ls())

dem <- terra::rast("~/Documents/Development/york/model/GIS/filled_dem.tif")
chn <- terra::rast("~/Documents/Development/york/model/GIS/channel.tif")

## divide up dem by height
hbrk <- seq(0,720,by=5)
hclass <- terra::classify(dem,hbrk)
hc <- terra::as.matrix( hclass, wide=TRUE )

## divide up dem by destination - use wd8
d <- terra::as.matrix( dem, wide=TRUE )
ch <- terra::as.matrix( chn,  wide=TRUE )
            
## create destination matrix
dst <- ch

## distances and contour lengths
## distance between cell centres
rs <- terra::res( dem )
dxy <- rep(sqrt(sum(rs^2)),8)
dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)

            
## if we go up in height order then we must have looked at all lower
## cells first
idx <- order(d,na.last=NA)

n_to_eval <- length(idx)
it <- 1
verbose <- TRUE
if(verbose){
    print_step <- round(n_to_eval/20)
    next_print <- print_step
}else{
    next_print <- Inf
}

for(ii in idx){
    if( is.finite(dst[ii]) ){ next }
        
    ## it is not already determined
    jdx <- ii+delta
    gcl <- dcl*(d[jdx]-d[ii])/dxy
    kdx <- which.min(gcl)
    dst[ii] <- dst[ jdx[kdx] ]

    ## verbose output here
    if(it >= next_print){
        cat(round(100*it / n_to_eval,1),
            "% complete","\n")
        next_print <- next_print+print_step
    }
    it <- it+1
}
dest <- terra::rast( dem, names="dest", vals=dst)


## combine hclass and dest
cp <- as.vector(hc)
x <- as.vector(dst)
cp <- 0.5*(cp+x)*(cp+x+1)+x
cp[ is.finite(as.vector(ch)) ] <- NA
ucp <- sort(unique(as.vector(cp)))
cp <- findInterval(cp, (ucp[-1]+head(ucp,-1))/2)
cp <- matrix(cp,nrow(hc),ncol(hc))
#cp[is.finite(ch)] <- NA
cl <- terra::rast( dem, names="class", vals=cp)
terra::plot(cl)
