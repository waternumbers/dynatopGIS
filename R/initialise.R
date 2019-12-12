#' Initialise the analysis by loading a dem
#'
#' @description Reads the DEM file, sets up the raster brick and saves it
#'
#' @param dem_file name of file to load as the dem
#' @param analysis_file file to save the analysis to
#'
#' @return raster brick object pointing to the analysis file
#'
#' @details Currently reads and then writes file using the raster package. No checks as to the projection are enforced. Output is project_path/dem.tif
#' @export
initialise <- function(dem,channel,analysis_file=character(0),fill_na=TRUE,...){

    ## read in dem
    if(!("RasterLayer" %in% class(dem))){
        dem <- raster::raster(dem)
    }

    ## read in channel
    ## if(!("???" %in% class(channel))){
    ##     channel <- load_channel(channel,...)
    ## }

    ## fill dem so only NA values are on boundaries
    if(fill_na){
        ## so all catchment dem value a number
        na_clumps <- clump(is.na(dem)) # clumps of na values
        edge_values <- setdiff( unique(c(na_clumps[1,],
                                         na_clumps[nrow(na_clumps),],
                                         na_clumps[,1],
                                         na_clumps[,ncol(na_clumps)])),
                               NA) #those clumps on the edge to be ignored
        na_clumps[na_clumps%in%edge_values] <- NA # set to NA to ignore
        dem[!is.na(na_clumps)] <- -Inf
    }

    ## initialise the other raster outputs
    land_area <- dem; land_area[] <- 4 ##TO DO change back after debugging raster::area(dem)
    channel_area <- dem; channel_area[] <- 0;
    channel_id <- dem; channel_id[] <- NA

    ## extract cells index, and fraction of river area in cell
    ch_cell<- raster::extract(land_area,channel,weights=TRUE,cellnumbers=TRUE,na.rm=TRUE)

    ## Loop and add the chanel id and channel area
    ch_area <- raster::area(channel) # areas of river channels
    for(ii in 1:length(ch_cell)){
        ch_cell[[ii]] <- data.frame( id = channel[['id']][ii],
                                    land_area = ch_cell[[ii]][,'value'],
                                    channel_area = ch_area[ii]*ch_cell[[ii]][,'weight'],
                                    cell = ch_cell[[ii]][,'cell'],
                                    stringsAsFactors=FALSE)
    }

    ## merge the list
    ch_cell <- do.call(rbind,ch_cell)
    
    ## assign the cells with multiple chanel id to the
    ## one with the largest area in the cell (equals returns first in list)
    ch_cell <- do.call(rbind,
                   lapply(split(ch_cell,ch_cell$cell),
                          function(chunk) chunk[which.max(chunk$channel_area),]))

    ## make land an channel area in ch_cell consistent
    ch_cell[,'channel_area'] <- pmin( ch_cell[,'channel_area'], ch_cell[,'land_area'] )
    ch_cell[,'land_area'] <- ch_cell[,'land_area'] - ch_cell[,'channel_area']

    ## add to rasters and write out
    land_area[ch_cell$cell] <- ch_cell$land_area
    channel_area[ch_cell$cell] <- ch_cell$channel_area
    channel_id[ch_cell$cell]  <- ch_cell$id
    
    brck <- raster::brick(list(dem=dem,
                               land_area=raster::mask( land_area, dem ),
                               channel_area=raster::mask( channel_area,dem),
                               channel_id=raster::mask( channel_id,dem)
                               )
                          )
    if(length(analysis_file)>0){
        writeRaster(brck,paste0(analysis_file))
    }
    
    return(brck)
}
