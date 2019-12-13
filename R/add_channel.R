#' Add a channel to the gis_brick
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
add_channel <- function(brck,channel){
    
    ## pass to checks
    check_brck(brck)
    check_chn(channel)
    
    ## extract cells index, and fraction of river area in cell
    ch_cell<- raster::extract(brck[['land_area']],channel,weights=TRUE,cellnumbers=TRUE,na.rm=TRUE)

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
    brck[['land_area']][ch_cell$cell] <- ch_cell$land_area
    brck[['channel_area']][ch_cell$cell] <- ch_cell$channel_area
    brck[['channel_id']][ch_cell$cell]  <- ch_cell$id
        
    return(brck)
}
