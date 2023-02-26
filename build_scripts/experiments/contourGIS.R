## test code for channel based sequences
rm(list=ls())

rst <- terra::rast(c("/home/smithpj1/Documents/local_share/eden/model/filled_dem.tif",
                     "/home/smithpj1/Documents/local_share/eden/model/channel.tif"))

## distances and contour lengths
## distance between cell centres
rs <- terra::res( rst )
dxy <- rep(sqrt(sum(rs^2)),8)
dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
nr <- nrow(rst); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)


## contour the dem
tmp <- terra::as.contour(rst[["filled_dem"]], maxcells=terra::ncell(rst[["filled_dem"]]), levels=seq(8,1000,by=10))
rst[["cntr"]] <- terra::rasterize(tmp,rst,field = "level",touches=TRUE)

## remove contour cells whose dominant flow path is to cells on the same level
ctr <- terra::as.matrix( rst[["cntr"]], wide=TRUE)
d <- terra::as.matrix( rst[["filled_dem"]], wide=TRUE)

idx <- which(is.finite(ctr))
for(ii in idx){
    jdx <- ii+delta
    gcl <- (d[jdx]-d[ii])*dcl/dxy
    kdx <- which.min(gcl)
    if( is.finite( ctr[ jdx[kdx] ] ) ){
        ctr[ii] <- NA
    }
}
rst[["cntr"]] <- ctr

## compute level bands by splitting topography

## single pass grouping... top to bottom
## - if NA
##   - if higher cell then take from one with largest flow direction
##   - if not start a new group
## - work out d/s cell
##   - if d/s cell river do nothing
##   - if d/s cell in different level start a new group
##   - if d/s cell aleady has a group then merge group numbers

## compute properties for each unit - determin classes of other bit by modal value

grp <- d*NA

idx <- order(d,decreasing=TRUE,na.last=NA)

for(ii in idx){
    


## work out cell that is d/s using dominant flow direction
dfd <- d*NA
ch <- terra::as.matrix( rst[["channel"]], wide=TRUE)

idx <- order(d,na.last=NA)

for(ii in idx){
    if( is.finite(ch[ii]) ){
        dfd[ii] <- ii
    }else{
        
        ## it is not already in a group
        jdx <- ii+delta
        gcl <- (d[jdx]-d[ii])*dcl/dxy
        kdx <- which.min(gcl)
        dfd[ii] <- jdx[kdx]
    }
}

rst[["dfd"]] <- dfd

## work out group
idx <- which(is.finite(ctr))
grp <- d*NA
grp[idx] <- (1:length(idx)) + max(ch,na.rm=TRUE)
idx <- which(is.finite(ch))
grp[idx] <- ch[idx]

idx <- order(d,na.last=NA)
for(ii in idx){
    if( is.finite(grp[ii]) ){ next }
            
    grp[ii] <- grp[ dfd[ii] ]
}
rst[["grp"]] <- grp

## make a table to populate in the next bit
tbl <- data.frame(grp=sort(unique(as.vector(grp))),
                  isSource = FALSE,
                  isConnected = FALSE)
tbl$isChannel <- tbl$grp <= max(ch,na.rm=TRUE)


## work out groups that are drained into from a contour line

idx <- which(is.finite(ctr))
