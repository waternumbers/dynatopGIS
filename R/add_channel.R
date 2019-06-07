#' Add a channel to the DEM
#'
#' @description Adds a chanel to the DEM by working out the land and channela areas and channel id
#'
#' @param project_path the path to the folder used for the analysis
#'
#' @details intersects the DEM and river network to determine which of the raster DEM cells contain parts of the river network. From this three raster maps are produced containing:
#' * the area in the DEM cell covered by land
#' * the area in the DEM cell covered by channel
#' * the id of the channel
#' If multiple river lengths intersect a DEM cell the id of that with the largest intersection is used. 
add_channel <- function(project_path){
    
    ## read in dem
    dem <- raster::raster(file.path(project_path,'dem.tif'))

    ## compute tha area of each pixel
    land_area <- mask( raster::area(dem), dem )

    ## read in the channel network
    channel <- raster::shapefile(file.path(project_path,'channel'))
    
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
    ch_cell[,'channel_area'] <- pmax( ch_cell[,'channel_area'], ch_cell[,'land_area'] )
    ch_cell[,'land_area'] <- ch_cell[,'land_area'] - ch_cell[,'channel_area']

    ## add to rasters and write out
    land_area[ch_cell$cell] <- ch_cell$land_area
    writeRaster(land_area,file.path(project_path,'land_area.tif'),format='GTiff')

    tmp <- dem; tmp[] <- NA
    tmp[ch_cell$cell] <- ch_cell$channel_area
    writeRaster(tmp,file.path(project_path,'channel_area.tif'),format='GTiff')
    tmp[ch_cell$cell] <- ch_cell$id
    writeRaster(tmp,file.path(project_path,'channel_id.tif'),format='GTiff')

    return(TRUE)
}
