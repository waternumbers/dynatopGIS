#' R6 Class for processing a catchment to make a Dynamic TOPMODEL
dynatopGIS <- R6::R6Class(
    "dynatopGIS",
    public = list(
        #' @description Reads the DEM file, sets up the raster catchment description
        #'
        #' @param dem a RasterLayer containing the dem
        #' @param fill_na should NA values in dem be filled. See details
        #'
        #' @details Reads in the dem using the raster package. If fill_na all NA values other then those that link to the edge of the dem are filled with Inf, so they can be identified as sinks.
        #' @return A new `dynatopGIS` object
        initialize = function(dem, fill_na=TRUE){
            #browser()
            ## check the file is a raster layer
            if(!("RasterLayer" %in% class(dem))){
                stop("dem should be a RasterLayer")
            }
            
            ## check the projection is OK
            if( !all.equal(diff(raster::res(dem)),0) | raster::isLonLat(dem) ){
                stop("Processing currently only works on projected dem's with a square grid")
            }
            
            ## fill dem so only NA values are on boundaries
            if(fill_na){
                ## so all catchment dem value a number
                na_clumps <- raster::clump(is.na(dem)) # clumps of na values
                edge_values <- setdiff( unique(c(na_clumps[1,],
                                                 na_clumps[nrow(na_clumps),],
                                                 na_clumps[,1],
                                                 na_clumps[,ncol(na_clumps)])),
                                       NA) #those clumps on the edge to be ignored
                na_clumps[na_clumps%in%edge_values] <- NA # set to NA to ignore
                dem[!is.na(na_clumps)] <- -Inf # set to low value to indicate missing
            }
            
            ## populate the data storage
            private$scope <- list(crs = raster::crs(dem),
                                  ext = raster::extent(dem),
                                  res = raster::res(dem),
                                  dim = dim(dem)[1:2])
            
            private$layers=list() # this will clear other data
            private$layers$dem <- raster::getValues(dem)
            private$layers$land_area <- rep(prod(private$scope$res),length(private$layers$dem))*(private$layers$dem>-Inf)
            
            ## return so can be chained
            invisible(self)
        },

        #' @description List the available of geographical layers
        #' @param is_class should layers of partial classification be returned
        ## @return a character vector of names
        list_layers =function(is_class=FALSE){
            if( is_class ){
                names(provate$class$partial)
            }else{
                names(private$layers)
            }
        },

        #' @description Add a layer of geographical information
        #' @param raster_layer a RasterLayer object to add
        #' @param layer_name name to give to the layer
        #' @param check_name logical indicating if name should be checked against reserved names
        #' @details The default is to check the name so that names reserved for computed varaibles are not overwritten. Disabling to overwrite computed layer may have unintended consequences.
        #' @return suitable for chaining
        add_layer = function(raster_layer,layer_name,check_name=TRUE){
            ## check it is a RasterLayer
            if(!("RasterLayer" %in% class(raster_layer))){
                stop("raster_layer should be a RasterLayer")
            }
            ## check crs etc
            if( !compareCRS(crs(raster_layer),private$scope$crs) ){stop("Projection of Rasterlayer does not match")}
            if( raster::extent(raster_layer) != private$scope$ext ){stop("Extent of Rasterlayer does not match")}
            if( !all(raster::res(raster_layer) == private$scope$res) ){stop("Resolution of Rasterlayer does not match")}

            ## check name
            if(check_name){
                layer_name <- as.character(layer_name)
                if( layer_name %in% private$reserved_layer_names ){stop("Layer name is reserved")}
                if( layer_name %in% names(private$layers) ){stop("Layer name is already in use")}
            }
                
            ## add layer
            private$layers[[layer_name]] <- raster::getValues(raster_layer)
            
            invisible(self)
        },
 
        #' @description Get a layer of geographical information
        #' @param layer_name name of the layer give to the layer
        #' @param is_class logical if the class based on this layer should be returned
        #' @return a RasterLayer of the requested information
        #' @details Setting is_class to TRUE and not providing a layer_name returns the overall classification
        get_layer = function(layer_name=character(0), is_class=FALSE){
            if( is_class ){
                if(length(layer_name)==0){
                    x <- private$class$total
                }else{
                    if( !(layer_name %in% names(private$class$partial)) ){
                        stop("Classification does not exist")}
                    x <- private$class$partial[[layer_name]]
                }
            }else{
                if( !(layer_name %in% names(private$layers)) ){
                    stop("layer does not exist")}
                x <- private$layers[[layer_name]]
            }
            
            raster(crs = private$scope$crs,
                   ext = private$scope$ext,
                   resolution = private$scope$res,
                   vals = x)
        },
        #' @description Plot a layer
        #' @param layer_name the name of layer to plot
        #' @param add_channel shouel the channel be added to the plot
        #' @param is_class logical if the class based on this layer should be returned
        #' @details Setting is_class to TRUE and not providing a layer_name returns a plot of the overall classification
        #' @return a plot
        plot = function(layer_name=character(0),add_channel=TRUE,is_class=FALSE){
            raster::plot( self$get_layer(layer_name,is_class) )
            if( add_channel ){ raster::plot( self$get_channel(), add=TRUE ) }
        },
        #' @description Import channel data from an OGR file to the `dynatopGIS` object
        #'
        #' @param sp_object a spatialLinesDataFrame containing the channel information
        #' @param property_names named vector of columns of the spatial data frame to use for channel properties
        #' @param default_width default width of a channel if not specified in property_names. Defaults to 2 metres.
        #'
        #' @details Takes the input file and covnerts it to polygons with properties length,width,startNode,endNode. The variable names in the sp_object data frame which corresponding to these properties can be specified in the \code{property_names} vector.
        #' The function opulates the channel_id & channel_area layers and alters the land_area layer to correspond.
        #'
        #' @return suitable for chaining
        add_channel = function(sp_object,property_names=c(length="length",
                                                           startNode="startNode",
                                                           endNode="endNode",
                                                           width="width"),default_width=2){
            
            ## check sp_object is a SpatialLinesDataFrame
            if(!is(sp_object,"SpatialLinesDataFrame") && !is(sp_object,"SpatialPolygonsDataFrame")){
                stop("sp_object is not a spatial lines or polygon data frame object (even when read in)")
            }
            
            ## check project of the sp_object
            if( !compareCRS(crs(sp_object),private$scope$crs) ){
                stop("Projection of channel object does not match")
            }
            
            ## check property_names is of correct type
            if(!is.vector(property_names) ||
               length(names(property_names))!=length(property_names) ){
                stop("property_names of wrong type")
            }

            ## check we have correct names
            if( !all( c("length","startNode","endNode") %in% names(property_names)) ){
                stop("A required property name is not specified")
            }
            
            ## check if there is an id feild which will be overwritten
            if( ("id" %in% names(property_names)) ){
                warning("The name id is reserved and will be overwritten")
                property_names <- property_names[-"id"]
            }
            
            ## see if all names exist
            if(!all( property_names %in% names(sp_object))){
                stop("Not all feilds names in property_names are in sp_object")
            }
            
            ## create channel object
            chn <- sp_object
            chn <- chn[,property_names]
            names(chn) <- names(property_names)
            chn[['id']] <- 1:length(chn)
            
            ## convert factors to strings
            idx <- sapply(chn@data, is.factor)
            chn@data[idx] <- lapply(chn@data[idx], as.character)
            
            if(!is(chn,"SpatialPolygonsDataFrame")){
                if(!("width" %in% names(chn))){
                    warning("Modifying to spatial polygons using default width")
                    chn[['width']] <- default_width
                }else{
                    warning("Modifying to spatial polygons using provided width")
                }
                
                ## buffer to convert to spatial polygons object
                chn <- rgeos::gBuffer(chn, byid=TRUE, width=chn[['width']])
            }

            ## store channel
            private$channel <- chn

            private$apply_add_channel()
            ## return so can be chained
            invisible(self)
        },
        #' @description Get the channel description
        #' @return an SpatialPolygonDataFrame containing the channels 
        get_channel = function(){
            private$channel
        },
        #' @description The sink filling algorithm of Planchona and Darboux 2001
        #'
        #' @param min_grad Minimum gradient between cell centres
        #' @param max_it maximum number of replacement cycles
        #' @param verbose print out additional diagnostic information
        #' @param hot_start start from filled_dem if it exists
        #' @details The algorithm implimented in Planchona and Darboux, "A fast, simple and versatile algorithm to fill the depressions in digital elevation models" Catena 46 (2001). A pdf can be found at https://horizon.documentation.ird.fr/exl-doc/pleins_textes/pleins_textes_7/sous_copyright/010031925.pdf.
        #'
        sink_fill = function(min_grad = 1e-4,max_it=1e6,verbose=FALSE, hot_start=FALSE){
            ## check for a required layers
            rq <- c("dem","channel_id","land_area")
            if(hot_start){rq <- c(rq,"filled_dem")}
            if(!all(rq %in% names(private$layers))){
                stop(paste(c("Algorithm requires the following layers",rq),
                           collapse=" "))
            }
                
            private$apply_sink_fill(min_grad,max_it,verbose,hot_start)
            invisible(self)
        },
        #' @description Computes statistics e.g. log(a/tanb) for raster cells
        #'
        #' @param missing_grad gradient that can be assigned to a pixel if it can't be computed
        #' @param verbose print out additional diagnostic information
        #' @details The algorithm works in two passes. The first computes the number of upstream pixels and values that do not depend on ordering. The second sequences downslope to compute values. Missing gradient is only used for pixels which are partially channel but have no upslope neighbours.
        compute_properties = function(missing_grad = 1e-4,verbose=FALSE){
            rq <- c("filled_dem","channel_id",
                    "land_area","channel_area")
            if(!all(rq %in% names(private$layers))){
                stop(paste(c("Algorithm requires the following layers",rq),
                           collapse=" "))
            }
            private$apply_compute_properties(missing_grad,verbose)
            invisible(self)
        },
        #' @description Split a catchment into a set hydrological response units (HRUs) according to any number of landscape layer cuts or burns
        #' @param cuts A named list of cuts of make to form the HRU. Names should correspond to raster layers in the project directory. Values should be numeric and define either the number of bands (single value) or breaks between band (multiple values)
        #' @param burns a vector of layer names which are to be used as burns
        #' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. Burns are added directly in the order they are given.
        classify = function(cuts,burns=NULL){

            ## check cuts
            cuts <- as.list(cuts)
            cuts <- lapply(cuts,FUN=function(x){as.numeric(x)})
            cut_names<- names(cuts)
            if( length(cut_names) != length(cuts) ){
                stop("Unnamed layers in the cuts")
            }
            all_names <- c(burns,cut_names)
            if( length(all_names)!=length(unique(all_names)) ){
                stop("Layer names can not be repeated")
            }
            if( !all(all_names %in% names(private$layers)) ){
                stop(paste( "Missing layers:",
                           paste(setdiff(cut_names,names(ctch$layers)),collapse=" ")))
            }
            
            private$apply_classify(cuts,burns)
        },
        #' @description Compute a Dynamic TOPMODEL
        #'@return a Dynamic TOPMODELGIS description suitable for the \code{dyantop} package
        create_model = function(){
            if( length(private$class$total) != length(private$layers$dem) ){
                stop("There is no valid classification")
            }
            
            private$apply_create_model()
        },
        
        #' @description Return the index of neighbouring cells
        #' @param x index of cells for which to find neighbours
        #' @return a list of indexes of neighbours
        neighbour = function(x){
            ## check if data loaded
            lapply(x,function(i){private$fN(i)$idx})
        },
        #' @description get the cuts and burns used to classify
        #' @return a list with two elements, cuts and burns
        get_class = function(){
            private$class$method
        },
        
            
        #' @description get the version number
        #' @return a numeric version number
        #' @details teh version number indicates the version of the algorithms within the object
        get_version = function(){
            private$version
        }
               
    ),
    private = list(
        version = 0.1,
        reserved_layer_names=c("dem","filled_dem","land_area",
                               "channel_area","channel_id",
                               "atanb","gradient","upslope_area","band"),
        scope=list(),
        layers=list(),
        class=list(total=NULL,
                   method=list(),
                   partial=list()),
        channel=NULL,
        
        apply_add_channel = function(){
            ## create a land area raster
            land_area <- self$get_layer("land_area")
            
            ## extract cells index, and fraction of river area in cell
            ch_cell<- raster::extract(land_area,private$channel,weights=TRUE,cellnumbers=TRUE,na.rm=TRUE)

            ## Loop and add the chanel id and channel area
            ch_area <- raster::area(private$channel) # areas of river channels
            for(ii in 1:length(ch_cell)){
                ch_cell[[ii]] <- data.frame( id = private$channel[['id']][ii],
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
            
            ## add to the catchment
            private$layers$channel_area <- rep(0,prod(private$scope$dim))
            private$layers$channel_id <- rep(NA,prod(private$scope$dim))
            private$layers$land_area[ch_cell$cell] <- ch_cell$land_area
            private$layers$channel_area[ch_cell$cell] <- ch_cell$channel_area
            private$layers$channel_id[ch_cell$cell]  <- ch_cell$id
        },
        apply_sink_fill = function(min_grad,max_it,verbose,hot_start){
    
            ## extract dem and channel_id values
            if(hot_start){
                w <- provate$layer$filled_dem
            }else{
                w <- private$layers$dem
                w[!is.na(w)] <- Inf
                w[is.finite(private$layers$channel_id)] <- private$layers$dem[is.finite(private$layers$channel_id)]
            }
            
            finite_neighbour <- rep(FALSE,length(w))
            for(ii in which(is.finite(w))){
                finite_neighbour[private$fN(ii)$idx] <- TRUE
            }
            
                                        # idx <- which(!is.na(w))
            
            something_done <- TRUE
            it <- 1
            eflag <- FALSE
            while( something_done ){
                something_done <- FALSE
                idx <- which( (w > private$layers$dem) & finite_neighbour )

                if( verbose ){
                    cat("Iteration",it,"\n")
                    cat("\t","Cells to evaluate:",length(idx),"\n")
                    cat("\t","Percentage Complete:",
                        round(100*sum(is.finite(w))/sum(is.finite(private$layers$dem)),1),"\n")
                }
                
                for(ii in idx){
                    ## compute neighbours
                    jj <- private$fN(ii)
                    w_min <- min(w[jj$idx]+jj$dx*min_grad,na.rm=TRUE)
            
                    if( private$layers$dem[ii] > w_min ){
                        ## set to private$layers$dem and leave
                        w[ii] <- private$layers$dem[ii]
                        something_done <- TRUE
                        finite_neighbour[jj$idx] <- TRUE
                    }else{
                        if( w[ii] > w_min ){
                            w[ii] <- w_min
                            finite_neighbour[jj$idx] <- TRUE
                            something_done <- TRUE
                        }
                    }
                }
                if( it > max_it ){
                    eflag <- TRUE
                    something_done <- FALSE # cause end of while loop
                }
                it <- it+1                
            }
            
            ## copy back into storage
            private$layers$filled_dem <- w
            
            if(eflag){ stop("Maximum number of iterations reached, sink filling not complete") }

        },

        ## Function to compute the properties
        apply_compute_properties = function(missing_grad,verbose){

            ## initialise the new layers to be generated
            new_layers <- list()
            new_layers$gradient <- new_layers$band <- new_layers$atanb <- rep(NA,prod(private$scope$dim))
            new_layers$upslope_area <- private$layers$land_area

            ## initialise the local storage vectors
            n_higher <- rep(NA,prod(private$scope$dim)) ## number of higher cells
            grad_cl <- rep(NA,prod(private$scope$dim)) ## sum of gradient * contour length of flow paths
            
            ## First pass to compute the number of higher cells, flow directions
            ## and variables that are not order specific
            cat("Pass 1","\n")

            idx <- which(!is.na(private$layers$filled_dem))
            for(ii in idx){
                
                ## cells with no land area shouldn't have values
                if( private$layers$land_area[ii] <= 0 ){
                    next
                }
                
                ## if cell has a channel_id then set flag
                is_channel <- is.finite(private$layers$channel_id[ii])

                ## direction properties
                jj <- private$fGCL(ii)

                ## number of higher cells
                n_higher[ii] <- length(jj$is_higher) # na.rm ensure cells on edge are evaluated

                ## work out the flow directions
                if( is_channel ){
                    ## use upslope cells
                    if( n_higher[ii] > 0 ){
                        ## use upslope cells for gradient if possible
                        g <- -jj$grd[jj$is_higher]
                        cl <- jj$cl[jj$is_higher]
                    }else{
                        ## else use mingrad with warning
                        g <- missing_grad
                        cl <- private$scope$res[1]
                        #warning(paste("Cell",ii,"is a channel cell with no higher neighbours"))
                    }
                }else{
                    ## cells with have no channel element
                    if( length(jj$is_lower) > 0 ){
                        ## use downslope cells for gradient if possible
                        g <- jj$grd[jj$is_lower]
                        cl <- jj$cl[jj$is_lower]
                    }else{
                        ## else use higher cells with warning
                        stop(paste("Cell",ii,"is a hillslope cell with no lower neighbours"))
                    }
                }
                gcl <- g*cl
                grad_cl[ii] <- sum(gcl)
                new_layers$gradient[ii] <- grad_cl[ii] / sum(cl)
                
            }

            
            
            
            ## Second pass to compute ordered valraibles
            cat("Pass 2","\n")         
            
            ## work down list of higher cells
            idx <- which(n_higher==0)
            cnt <- 0
            while( length(idx)>0 ){
                cnt <- cnt + 1
                if( verbose ){
                    print(paste("Band",cnt))
                }
            
                new_layers$band[idx] <- cnt
                n_higher[idx] <- -1
                idx_list <- list()
                icnt <- 0
                for(ii in idx){
                                        
                    ## move upslopre area downslope
                    fd <- private$fFD(ii)
                    kk <- fd$idx # private$layers$flowDir[[ii]]$idx
                    ff <- fd$frc # private$layers$flowDir[[ii]]$frc
                    new_layers$upslope_area[kk] <- new_layers$upslope_area[kk] +
                        new_layers$upslope_area[ii]*ff
                    n_higher[kk] <- n_higher[kk] - 1

                    ## since these have no known n_higher
                    nk <- kk[private$layers$land_area[kk]>0]
                    nk <- nk[n_higher[nk]==0]
                    
                    if(length(nk)>0){
                        icnt <- icnt+1
                        idx_list[[icnt]] <- nk
                    }
                    
                }

                
                ## find next set of cells to evaluate
                idx <- unique(do.call(c,idx_list))
                ## idx2 <- which(n_higher==0)
                ## if( !all(idx %in% idx2) | length(idx)!=length(idx2) ){
                ##     browser()
                ## }
                
            }
            
            ## compute tographic index fromm summaries
            
            new_layers$atanb <- log( new_layers$upslope_area / grad_cl )
            if( any(new_layers$atanb==Inf,na.rm=TRUE) ){       
                warning("None finite topographic index values produced - this requires investigation")
            }

            ## copy back into private store
            for(ii in names(new_layers)){
                private$layers[[ii]] <- new_layers[[ii]]
            }
            
            
        },

        ## split_to_class
        apply_classify = function(cuts,burns){
            
            ## work out new cuts by cantor_pairing
            nms <- names(cuts)
            for(ii in nms){
                brk <- as.numeric(cuts[[ii]])
                if( length(brk)==1 & is.na(brk[1]) ){
                    x <- private$layers[[ii]]
                }else{
                    x <- cut(private$layers[[ii]],breaks=as.numeric(cuts[[ii]]),
                             labels = FALSE)
                }
                private$class$partial[[ii]] <- x
                if(ii == nms[1]){
                    cp <- x
                }else{
                    cp <- 0.5*(cp+x)*(cp+x+1)+x
                }
            }

            ## add in burns
            for(ii in burns){
                idx <- is.finite(private$layers[[ii]])
                cp[idx] <- private$layers[[ii]][idx]
            }
            
            ## renumber from cantour paring to get sequential values
            ucp <- unique(cp,na.rm=TRUE)
            ucp <- ucp[!is.na(ucp)]
            ucp <- setNames( (1:length(ucp)) + max(private$layers$channel_id,na.rm=TRUE),
                            paste(ucp) ) # only works since cp is an integer
            cp <- ucp[paste(cp)]
            
            ## add back into layer
            private$class$total <- cp
            
            ## record in splits
            private$class$method <- list(cuts = cuts,burns=burns)
        },
    
        apply_create_model = function(){
            
            ## Initialise the model
            model <- list()
            
            ## create the model hillslope table and populate computations
            uid <- unique(private$class$total)
            uid <- uid[!is.na(uid)]
            model$hillslope <- data.frame(
                id = as.integer(uid),
                area = rep(NA,length(uid)),
                atb_bar = rep(NA,length(uid)),
                s_bar = rep(NA,length(uid)),
                delta_x = rep(NA,length(uid)),
                ## class = rep(NA,length(uid)),
                sz_band = rep(NA,length(uid)),
                precip="unknown",
                pet="unknown",
                q_sfmax="q_sfmax_default",
                s_rzmax="s_rzmax_default",
                s_rz0="s_rz0_default",
                ln_t0="ln_t0_default",
                m="m_default",
                t_d="t_d_default",
                t_sf="t_sf_default",
                stringsAsFactors=FALSE
            )
            
            ## populate the general properties
            for(id in model$hillslope$id){
                ## row of data.frame
                rw <- which(model$hillslope$id==id)
                
                ## index of hillslope points
                idx <- which(private$class$total==id)
                
                area <- sum( private$layers$land_area[idx] )
                model$hillslope$area[rw] <- area
                
                model$hillslope$atb_bar[rw] <- sum(
                    private$layers$atanb[idx]*private$layers$land_area[idx] )/area
                
                model$hillslope$s_bar[rw] <- sum(
                    private$layers$gradient[idx]*private$layers$land_area[idx] )/area
                
                model$hillslope$delta_x[rw] <- (diff(range(
                    private$layers$band[idx]))+1) * private$scope$res[1]
                
                model$hillslope$sz_band[rw]  <-  min( private$layers$band[idx] )
                
                for(jj in names(private$class$partial)){
                    model$hillslope[[ paste0(jj,"_class") ]] <- unique(
                        private$class$partial[[jj]][idx])
                }
            }
            
            ## do the saem for the channel table
            uid <- private$channel$id
            model$channel <- data.frame(
                id = uid,
                area = unname(tapply(private$layers$channel_area,
                                     private$layers$channel_id,
                                     sum)
                              [ paste(uid) ]),
                length = as.numeric(private$channel$length),
                precip="unknown",
                pet="unknown",
                v_ch="v_ch_default",
                stringsAsFactors=FALSE
            )
            
            ## loop for hillslope flow directions           
            receiving_id <- pmax(private$class$total,private$layers$channel_id,
                                 na.rm=TRUE)
            model$hillslope$sz_dir = rep(list(NULL),length(model$hillslope$id))
            
            bnd <- rep(NA,max(model$hillslope$id)) 
            bnd[model$hillslope$id] <- model$hillslope$sz_band
            bnd[model$channel$id] <- Inf #model$channel$sz_band
            
            for(id in model$hillslope$id){
                rw <- which(model$hillslope$id == id)
                
                ## index of hillslope points
                idx <- which(private$class$total==id)
                
                ## work out flow directions
                rfrc <- rep(0,max(model$hillslope$id))
                for(ii in idx){
                    jj <- private$fFD(ii)
                    jj$cls <- receiving_id[jj$idx]
                    rfrc[jj$cls] <- rfrc[jj$cls] + private$layers$land_area[ii]*jj$frc
                }
                
                kdx <- which( (bnd > bnd[id]) & rfrc>0 )
                model$hillslope$sz_dir[[rw]] <- list(
                    band = model$hillslope$sz_band[[rw]],
                    idx = kdx,
                    frc = rfrc[kdx]/sum(rfrc[kdx])
                )
            }

            ## handle HSU with no downslope by assigning to same class in later bands
            keep <- rep(TRUE,nrow(model$hillslope))
            for(ii in which(sapply(model$hillslop$sz_dir$idx,length)==0)){
                jdx <- which(model$hillslop$sz_band > model$hillslop$sz_band[ii] &
                             model$hillslope$class == model$hillslope$class[ii])
                jdx <- jdx[ which.min(model$hillslop$sz_band[jdx] )]
                if(length(jdx)>0){
                    ## merge to later class
                    model$hillslope$area
                    atb_bar
                    s_bar
                    
                }else{
                    warning("Cannot reassign HSU",model$hillslope$id[ii],"dropping")
                    keep[ii] <- FALSE
                }
            }
            
            ## set surface to saturated for lack of better knwoledge
            model$hillslope$sf_dir <- model$hillslope$sz_dir

            ## remove band from mode
            model$hillslope$sz_band <- NULL
                
            ## work out the channel flow directions
            model$channel$flow_dir = rep(list(NULL),length(model$channel$id))
            for(id in model$channel$id){
                rw <- which( model$channel$id ==id )

                jj <- which( private$channel$id == id )
                
                jdx <- which(private$channel[['startNode']]==
                             private$channel[['endNode']][jj])
                if( length(jdx) > 0 ){
                    model$channel$flow_dir[[rw]] <- list(idx=jdx,
                                                         frc= 1/length(jdx))
                }
            }
            
            ## parameter values
            model$param <- c(q_sfmax_default=Inf,
                             s_rzmax_default=0.05,
                             s_rz0_default=0.99,
                             ln_t0_default=19,
                             m_default=0.004,
                             t_d_default=20,
                             t_sf_default=100,
                             v_ch_default=100)
            
            ## ############################################
            ## Add gauges at all outlets from river network
            ## ############################################
            
            idx <- which(sapply(model$channel$flow_dir,length)==0)
            model$gauge <- data.frame(
                name = paste("channel",model$channel$id[idx],sep="_"),
                id = model$channel$id[idx],
                fraction = rep(1,length(idx)),
                stringsAsFactors=FALSE
            )
            
            ## ##################################
            ## Add point inflow table
            ## ##################################
            ## blank point inflow table
            model$point_inflow <- data.frame(
                name = character(0),
                id = integer(0),
                fraction = numeric(0),
                stringsAsFactors=FALSE
            )
            
            return(model)
            
        },
        
        
        fN = function(k,pflg=FALSE){
            nr <- private$scope$dim[1]
            nc <- private$scope$dim[2]
            if( k > nc*nr | k < 1){ return(rep(NA,8)) }
            
            ## find i,j location
            j <- ((k-1)%/%nc) +1 # row
            i <- k - (j-1)*nc #column

            ## distance
            dx <- private$scope$res[1]
            dxd <- dx*sqrt(2)
            dst <- c(dxd,dx,dxd,dx,dx,dxd,dx,dxd)
            ## contour length
            cl <- dx*c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)
            ## make a matrix
            ii <- i + c(-1,0,1,-1,1,-1,0,1)
            jj <- j + c(-1,-1,-1,0,0,1,1,1)
            idx <- ii>0 & ii<=nc & jj>0 & jj<=nr
            out <- list(idx = (jj[idx]-1)*nc + ii[idx],
                        dx = dst[idx],
                        cl=cl[idx])
            
            ## m <- matrix(c(i-1,i,i+1,i-1,i+1,i-1,i,i+1,
            ##               j-1,j-1,j-1,j,j,j+1,j+1,j+1,
            ##               dst,
            ##               cl),8,4)
            
            ## ## trim to within bounds
            ## m <- m[ (m[,1]>0 & m[,1]<=nc & m[,2]>0 & m[,2]<=nr) ,]
            
            ## out <- list(idx = (m[,2]-1)*nc + m[,1],
            ##             dx = m[,3],
            ##             cl = m[,4])
            return(out)
        },
        fGCL = function(k){
            ## index of neighbours and distances
            out <- private$fN(k)
            ## compute gradient time contour length
            out$grd <- (private$layers$filled_dem[k] - private$layers$filled_dem[out$idx])/out$dx
            out$is_higher <- which( (out$grd < 0) &
                                    !is.finite(private$layers$channel_id[out$idx]) )
            out$is_lower <- numeric(0)
            if( !is.finite(private$layers$channel_id[k]) ){
                out$is_lower <- which( out$grd > 0 )
            }
            
            return(out)
        },
        fFD =function(k){
            tmp <- private$fGCL(k)
            out <- list(idx=numeric(0),
                        frc=numeric(0))
            if( length(tmp$is_lower)>0 ){
                gcl <- tmp$grd[tmp$is_lower]*tmp$cl[tmp$is_lower]
                out$idx <- tmp$idx[tmp$is_lower]
                out$frc <- gcl/sum(gcl)
            }
            return(out)
        }

        
    )
)
