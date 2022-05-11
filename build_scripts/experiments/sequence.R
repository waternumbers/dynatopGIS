## test code for channel based sequences
rm(list=ls())
private <- list()
private$brk <- terra::rast(c("/home/smithpj1/Documents/local_share/eden/model/filled_dem.tif",
                      "/home/smithpj1/Documents/local_share/eden/model/channel.tif"))

private$shp <- terra::vect("/home/smithpj1/Documents/local_share/eden/model/channel.shp")

rq <- c("filled_dem","channel")
if(!all( rq %in% names( private$brk) )){
    stop("Not all required input layers have been generated \n",
         "Try running sink_fill first")
}


d <- terra::as.matrix( private$brk[["filled_dem"]], wide=TRUE )
ch <- terra::as.matrix( private$brk[["channel"]],  wide=TRUE )

sq <- d*NA

## do channel bands
cnt <- 1
private$shp <- private$shp[  order(private$shp$id)  , ]
sqv <- rep(NA,nrow(private$shp))

n_to_eval <- nrow(private$shp)
verbose=TRUE
it <- 1
if(verbose){
    print_step <- round(n_to_eval/20)
    next_print <- print_step
}else{
    next_print <- Inf
}



ii <- 1
sqv[ii] <- 1
for( ii in 2:nrow(private$shp) ){
    
    idx <- private$shp$startNode == private$shp$endNode[ii]
    sqv[ii] <- max(sqv[idx])+1
    ## verbose output here
    if(it >= next_print){
        cat(round(100*it / n_to_eval,1),
            "% complete","\n")
        next_print <- next_print+print_step
    }
    
    it <- it+1
    
}

## do hillslope bands

## distances and contour lengths
## distance between cell centres
rs <- terra::res( private$brk )
dxy <- rep(sqrt(sum(rs^2)),8)
dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)


## if we go up in height order then we must have looked at all lower
## cells first
idx <- order(d,na.last=NA)
            
n_to_eval <- length(idx)
verbose=TRUE
it <- 1
if(verbose){
    print_step <- round(n_to_eval/20)
    next_print <- print_step
}else{
    next_print <- Inf
}

for(ii in idx){
    if( is.finite( ch[ii] ) ){
        sq[ii] <- sqv[ ch[ii] ]

    }else{
        jdx <- ii+delta
        gcl <- (d[jdx]-d[ii])*dcl/dxy
        kdx <- jdx[is_lower]
        sq[ii] <- max(sq[kdx])+1
    }

    ## verbose output here
    if(it >= next_print){
        cat(round(100*it / n_to_eval,1),
            "% complete","\n")
        next_print <- next_print+print_step
    }
    
    it <- it+1
}




## going up leaves "hanging" groups which don;t meet the upper contour which need merging

## cells first
idx <- order(d,na.last=NA)
            
n_to_eval <- length(idx)
verbose=TRUE
it <- 1
if(verbose){
    print_step <- round(n_to_eval/20)
    next_print <- print_step
}else{
    next_print <- Inf
}

for(ii in idx){
    if(is.finite(grp[ii])){ next }

    ## it is not already in a group
    jdx <- ii+delta
    gcl <- (d[jdx]-d[ii])*dcl/dxy
    kdx <- which.min(gcl)
    grp[ii] <- grp[ jdx[kdx] ]
    
    ## is_lower <- is.finite(gcl) & gcl<0
    ## kdx <- jdx[is_lower]
    ##                 bnd[ii] <- max(bnd[kdx])+1
    ##                 sfl[ii] <- min(sfl[kdx]+dxy[is_lower])
    ##                 efl[ii] <- sum((efl[kdx]+dxy[is_lower])*gcl[is_lower])/sum(gcl[is_lower])
    ##                 kdx <- which.min(gcl)
    ##                 dfl[ii] <- dfl[jdx[kdx]]+dxy[kdx]
    ##             }

    ## verbose output here
    if(it >= next_print){
        cat(round(100*it / n_to_eval,1),
            "% complete","\n")
        next_print <- next_print+print_step
    }
    
    it <- it+1
}


out <- terra::rast( private$brk[["filled_dem"]], names="sequence", vals=grp )


tmp <- table(ch)
quantile(tmp)
















## ## create some distance matrices
## sq <- d*NA

## ## distances and contour lengths
## ## distance between cell centres
## rs <- terra::res( private$brk )
## dxy <- rep(sqrt(sum(rs^2)),8)
## dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
## dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
## nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)

## cnt <- 1
## for(jj in sort(unique(as.vector(ch)))){ ## loop river segments
##     idx <- which(ch==jj) ## areas in that channel
##     sq[idx] <- cnt
##     idx <- unique(as.vector(sapply(idx, function(i){i + delta})))

##     it <- 1
##     while(length(idx)>0){
##         #print(paste(it,length(idx)))
##         ii <- idx[1]
##         ## if it is a channel remove and don't look
##         if(is.finite(ch[ii]) | is.na(d[ii]) | is.finite(sq[ii])){
##             idx <- idx[-1]
##             next
##         }
##         ## find neighbours
##         jdx <- ii+delta
##         gcl <- (d[jdx]-d[ii])*dcl/dxy
##         is_lower <- is.finite(gcl) & gcl<0
##         kdx <- jdx[is_lower]
##         tmp <- max(sq[kdx])+1
##         if( !is.na(tmp) ){
##             sq[ii] <- tmp
##             cnt <- max(cnt,tmp) + 1
##             idx <- unique( c(idx[-1], jdx[!is_lower]) )
##         }else{
##             idx <- idx[-1]
##         }
##         it <- it+1
##     }
##     print(paste(jj, range(sq,na.rm=TRUE)))
## }

## out <- terra::rast( private$brk[["filled_dem"]], names="sequence", vals=sq )
