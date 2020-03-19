#' Function to convert to and from a catchment object

glist_to_raster <- function(glist,layer="dem"){
    raster(crs = glist$raster$crs,
           ext = glist$raster$ext,
           resolution = glist$raster$res,
           vals = glist$layers[[layer]])
}

raster_to_glist <- function(rst,variable_name){
    out <- list(raster=list(crs = crs(rst),
                            ext = extent(rst),
                            res = res(rst),
                            dim = dim(rst)[1:2]),
                layers=list()
                )
    out$layers[[variable_name]] <- getValues(rst)
    return(out)
}

merge_glist <- function(...){
    all <- alist(...)
    out <- all[[1]]
    for(ii in 2:length(all)){
        ## TODO - this need to check projection properties
        jj <- setdiff(names(all[[ii]]),names(out))
        out[jj] <- all[[ii]][jj]
    }
    return(out)
}

check_glist <- function(glist){
    ## Check project info
    ##TODO

    ## check all layers the same and correct length
    rql <- prod(glist$dim)
    lngth <- sapply(glist$layers,length)
    if( !all( lngth == rql ) ){
        stop(paste("The following layer have the wrong length data",names(glist$layers)[lngth!=rql]))
    }
}

is_glist <- function(glist){
    ## Check it is a glist
    ## TODO
    TRUE
}
