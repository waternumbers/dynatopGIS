#' R6 Class for processing a catchment to make a Dynamic TOPMODEL
#' @examples
#' ## The vignettes contains more examples of the method calls.
#' 
#' ## create temport directory for output
#' demo_dir <- tempfile("dygis")
#' dir.create(demo_dir)
#'
#' ## initialise processing
#' ctch <- dynatopGIS$new(file.path(demo_dir,"meta.json"))
#'
#' ## add digital elevation and channel data
#' dem_file <- system.file("extdata", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
#' dem <- raster::raster(dem_file)
#' ctch$add_dem(dem)
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' sp_lines <- rgdal::readOGR(channel_file)
#' property_names <- c(channel_id="identifier",endNode="endNode",startNode="startNode",length="length")
#' ctch$add_channel(sp_lines,property_names)
#'
#' ## compute properties 
#' ctch$compute_areas()
#' ctch$sink_fill() ## fill sinks in the catchment
#' \donttest{
#' ctch$compute_properties() # like topograpihc index and contour length
#' ctch$compute_flow_lengths()
#' }
#' ## classify and create a model
#' \donttest{
#' ctch$classify("atb_20","atb",cuts=20) # classify using the topographic index
#' ctch$get_class_method("atb_20") ## see the details of the classification
#' ctch$combine_classes("atb_20_band",c("atb_20","band")) ## combine classes
#' ctch$create_model("new_model","atb_20_band","band") ## create a model
#' list.files(demo_dir,pattern="new_model*") ## look at the output files for the model
#' }
#' ## tidy up
#' unlink(demo_dir)
#' @export
dynatopGIS <- R6::R6Class(
    "dynatopGIS",
    public = list(
        #' @description Initialise a project, or reopen an existing project
        #'
        #' @param meta_file filename and path of the meta data file
        #' @param check logical, should checks be performed [TRUE]
        #' @param verbose printing of checking output [TRUE]
        #'
        #' @details This loads the meta data file found at \code{meta_path}, or creates it with a warning if no file is present. It \code{check} is \code{TRUE} then the meta data file contents are checked with the level of returned information being controlled by \code{verbose}.
        #'
        #' @return A new `dynatopGIS` object
        initialize = function(projFile){
            private$apply_initialize( tools::file_path_sans_ext(projFile) )
            invisible(self)
        },
        #' @description Get project meta data
        get_meta = function(){
            private$apply_get_meta()
        },
        #' @description Import a dem to the `dynatopGIS` object
        #'
        #' @param dem a \code{raster} layer object or the path to file containing one which is the DEM
        #' @param fill_na  should NA values in dem be filled. See details
        #' @param verbose Should additional progress information be printed
        #'
        #' @details If not a \code{raster} the DEM is read in using the raster package. If \code{fill_na} is \code{TRUE} all NA values other then those that link to the edge of the dem are filled so they can be identified as sinks.
        #'
        #' @return suitable for chaining              
        add_dem = function(dem,fill_na=TRUE){
            if(!("RasterLayer" %in% class(dem))){ dem <- raster::raster(as.character(dem)) }
            if(!("RasterLayer" %in% class(dem))){ stop("dem is not a RasterLayer") }
            private$apply_add_dem(dem,fill_na)
            invisible(self)
        },
        #' @description Import channel data to the `dynatopGIS` object
        #'
        #' @param channel a SpatialLinesDataFrame, SpatialPolygonsDataFrame or file path containing the channel information
        #' @param property_names named vector of columns of the spatial data frame to use for channel properties - see details
        #' @param default_width default width of a channel if not specified in property_names. Defaults to 2 metres.
        #'
        #' @details Takes the input channel converts it a SpatialPolygonDataFrame with properties length, startNode and endNode. The variable names in the sp_object data frame which corresponding to these properties can be specified in the \code{property_names} vector. In the channel is a SpatialLinesDataFrame (or read in as one) an additional property width is used to buffer the lines and create channel polygons. If required the width property is created using the default value. Note that any columns called length, startNode, endNode  and width are overwritten. Any column called id is copied to a column original_id then overwritten.
        #'
        #' @return suitable for chaining
        add_channel = function(channel){
            if(!is(channel,"SpatialPolygonsDataFrame")){ channel <- raster::shapefile( as.character(channel) ) }
            if(!is(channel,"SpatialPolygonsDataFrame")){ stop("channel is not a SpatialPolygonsDataFrame") }

            private$apply_add_channel(channel)
            invisible(self)
        },
        #' @description Add a layer of geographical information
        #'
        #' @param layer_name name to give to the layer
        #' @param file_path the location of the file containing the new layer
        #'
        #' @details The file given is read by the \code{raster} package and checked against the project meta data. Only layer names not already in use (or reserved) are allowed. If successful the meta data for the project are altered to reflect the new layer name and file location.
        #' @return suitable for chaining
        add_layer = function(layer,layer_name=names(layer)){
            layer_name <- as.character(layer_name)
            if(!("RasterLayer" %in% class(layer))){ layer <- raster::raster(as.character(layer)) }
            if(!("RasterLayer" %in% class(layer))){ stop("layer is not a RasterLayer") }
            private$apply_add_layer(layer,layer_name)
            invisible(self)
        },
        #' @description Get a layer of geographical information or a list of layer names
        #' @param layer_name name of the layer give to the layer
        #' @return a `raster` layer of the requested information if layer_name is given else a vector of layer names
        get_layer = function(layer_name=character(0)){

            ## handle case where a list of layers is requested
            if( length(layer_name) == 0 ){
                return( names(private$brk) )
            }
            ## check layer name exists
            layer_name <- match.arg(layer_name,names(private$brk))
            
            ## make raster and return
            if(layer_name=="channel"){
                return( private$shp )
            }else{
                return( private$brk[[layer_name]] )
            }
        },
        #' @description Plot a layer
        #' @param layer_name the name of layer to plot
        #' @param add_channel should the channel be added to the plot
        #' @return a plot
        plot_layer = function(layer_name,add_channel=TRUE){
            lyr <- self$get_layer(layer_name)
            raster::plot( lyr, main = layer_name)
            if( add_channel ){
                raster::plot(private$shp, add=TRUE )
            }
        },
        #' @description The sink filling algorithm of Planchona and Darboux (2001)
        #'
        #' @param min_grad Minimum gradient between cell centres
        #' @param max_it maximum number of replacement cycles
        #' @param verbose print out additional diagnostic information
        #' @param hot_start start from filled_dem if it exists
        #' @details The algorithm implemented is that described in Planchona and Darboux, "A fast, simple and versatile algorithm to fill the depressions in digital elevation models" Catena 46 (2001). A pdf can be found at (<https://horizon.documentation.ird.fr/exl-doc/pleins_textes/pleins_textes_7/sous_copyright/010031925.pdf>).
        #'
        sink_fill = function(min_grad = 1e-4,max_it=1e6,verbose=FALSE, hot_start=FALSE){
            private$apply_sink_fill(min_grad,max_it,verbose,hot_start)
            invisible(self)
        },
        #' @description Computes statistics e.g. gradient, log(upslope area / gradient) for raster cells
        #'
        #' @param min_grad gradient that can be assigned to a pixel if it can't be computed
        #' @param verbose print out additional diagnostic information
        #'
        #' @details The algorithm passed through the cells in decreasing height. Min grad is applied to all cells. It is also used for missing gradients in pixels which are partially channel but have no upslope neighbours.
        compute_properties = function(min_grad = 1e-4,verbose=FALSE){
            private$apply_compute_properties(min_grad,verbose)
            invisible(self)
        },
        #' @description Computes flow length for each pixel to the channel
        #'
        #' @param verbose print out additional diagnostic information
        #'
        #' @details The algorithm passed through the cells in increasing height. For measures of flow length to the channel are computed. The shortest length (minimum length to channel through any flow path), the dominant length (the length taking the flow direction with the highest fraction for each pixel on the path) and expected flow length (flow length based on sum of downslope flow lengths based on fraction of flow to each cell) and band (strict sequence to ensure that all contributing cell have a higher band value). By definition cells in the channel that have no land area have a length (or band) of NA.
        compute_flow_lengths = function(verbose=FALSE){
            private$apply_flow_lengths(verbose)
            invisible(self)
        }, 
        #' @description Create a catchment classification based cutting an existing layer into classes
        #' @param layer_name name of the new layer to create
        #' @param base_layer name of the layer to be cut into classes
        #' @param cuts values on which to cut into classes. These should be numeric and define either the number of bands (single value) or breaks between band (multiple values).
        #'
        #' @details This applies the given cuts to the supplied landscape layer to produce areal groupings of the catchment. Cuts are implement using \code{raster::cut} with \code{include.lowest = TRUE}. Note that is specifying a vector of cuts values outside the limits will be set to NA.
        classify = function(layer_name,base_layer,cuts){
            private$apply_classify(as.character(layer_name), as.character(base_layer), cuts)
            invisible(self)
        },
        #' @description Combine any number of classifications based on unique combinations and burns
        #' @param layer_name name of the new layer to create
        #' @param pairs a vector of layer names to combine into new classes through unique combinations. Names should correspond to raster layers in the project directory.
        #' @param burns a vector of layer names which are to be burnt on
        #'
        #' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. Burns are added directly in the order they are given. Cuts are implement using \code{raster::cut} with \code{include.lowest = TRUE}. Note that is specifying a vector of cuts values outside the limits will be set to NA.
        combine_classes = function(layer_name,pairs,burns=NULL){
            private$apply_combine_classes(as.character(layer_name),
                                          as.character(pairs),
                                          as.character(burns))
            invisible(self)
        },
        #' @description Compute a Dynamic TOPMODEL
        #'
        #' @param layer_name name for the new model and layers
        #' @param class_layer the layer defining the topographic classes
        #' @param dist_layer the layer defining the distances to the channel
        #' @param transmissivity transmissivity profile to use
        #' @param channel_solver channel solver to use
        #' @param dist_delta used in computing flow linkages see details        
        #' @param rain_layer the layer defining the rainfall inputs
        #' @param rain_label Prepended to rain_layer values to give rainfall series name        
        #' @param pet_layer the layer defining the pet inputs
        #' @param pet_label Prepended to pet_layer values to give pet series name
        #' @param verbose print more details of progress
        #'
        #' @details The \code{class_layer} is used to define the HRUs. Flow between HRUs is based on the distance to a channel. For each HRU the shortest distance to a channel is computed. Flow from a HRU can only go to a HRU with a lower shortest distance to the channel. Flow from a HRU can occur from any raster cell within the HRU whose distance to the channel is within dist_delta of the shortest distance within the HRU.
        #' Setting the transmissivity and channel_solver options ensure the model is set up with the correct parameters present.
        #' The \code{rain_layer} (\code{pet_layer}) can contain the numeric id values of different rainfall (pet) series. If the value of \code{rain_layer} (\code{pet_layer}) is not \code{NULL} the weights used to compute an averaged input value for each HRU are computed, otherwise an input table for the models generated with the value "missing" used in place of the series name.
        create_model = function(layer_name,class_layer,dist_layer,
                                transmissivity = c("exp","bexp","cnst","dexp"),
                                dist_delta=0,
                                rain_layer=NULL,rain_label=character(0),
                                pet_layer=NULL,pet_label=character(0),
                                verbose=FALSE){
            ## check valid transmissivity and channel_solver
            transmissivity <- match.arg(transmissivity)
            
            private$apply_create_model(class_layer, dist_layer,dist_delta,
                                       rain_layer, rain_label,
                                       pet_layer, pet_label,
                                       layer_name,verbose,
                                       transmissivity)

            invisible(self)
        },
        #' @description get the version number
        #' @return a numeric version number
        #' @details the version number indicates the version of the algorithms within the object
        get_version = function(){
            private$version
        },
        #' @description get the cuts and burns used to classify
        #' @param layer_name the name of layer whose classification method is returned
        #' @return a list with two elements, cuts and burns
        get_method = function(layer_name=character(0)){
            if(length(layer_name)==0){
                return( private$layers )
            }
            layer_name <- match.arg(layer_name,names(private$layers))
            return( private$layers[[layer_name]] )
        }               
    ),
    private = list(
        version = "0.3.0",
        projFiles = character(0),
        brk = character(0),
        shp = character(0),
        layers = list(),
        reserved_layers = c("dem","channel","filled_dem",
                            "gradient","upslope_area","atb",
                            "band","shortest_flow_length",
                            "dominant_flow_length","expected_flow_length"),
        writeTIF = function(){
            brk <- terra::rast(private$brk)
            terra::writeRaster(brk,private$projFiles["tif"],overwrite=TRUE)
            NULL
        },
        writeJSON = function(){
            writeLines( jsonlite::toJSON(private$layers), private$projFiles["json"] )
        },
        readJSON = function(fn){
            if(!file.exists(fn)){
                warning("No method json file found")
            }
            jsonlite::fromJSON(fn)
        },
        ## check and read the project files
        apply_initialize = function(projFile){
            ## compute and label file names
            fn <- setNames( paste0(projFile,c(".tif",".shp",".json")), c("tif","shp","json"))
            ## see if files exist
            fExists <- setNames( file.exists( fn ), c("tif","shp","json"))
            ## if no files then start a new project
            if( !any(fExists) ){
                print( paste("Starting new project at", projFile) )
                private$projFiles <- fn
                return(NULL)
            }
            ## if there is no raster data file fail
            if(!fExists["tif"]){ stop("No tif file of raster data") }
            ## read in netcdf and check
            brk <- raster::brick( fn["tif"] )
            if( !all.equal(diff(raster::res(brk)),0) | raster::isLonLat(brk) ){
                stop("tif file is not valid: Processing currently only works on projected data with a square grid")
            }
            if( !("dem" %in% names(brk)) ){ stop("No dem layer in the NetCDF file") }
            ## read in shape file and check
            if(fExists["shp"]){
                shp <- raster::shapefile( fn["shp"] )
                if(!compareCRS(shp,brk)){ stop("Shape and NetCDF projections do not match") }
                if( !("channel" %in% names(brk)) ){
                    stop("No channel layer in the NetCDF file but channel shapefile specified")
                }
            }
            if(fExists["json"]){
                mth <- private$readJSON(fn["json"])
            }
            

            private$brk <- brk
            if(fExists["shp"]){ private$shp <- shp }
            if(fExists["json"]){ private$layers <- mth }
            private$projFiles <- fn
        },
        ## adding dem
        apply_add_dem = function(dem,fill_na){
            
            if("dem" %in% names(private$brk) ){
                stop("The DEM exists, start a new project")
            }
            
            ## add to ensure NA on each edge
            dem <- extend(dem,c(1,1),NA)
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
                dem[!is.na(na_clumps)] <- -1e6+1 # set to low value to indicate missing
            }
            ## covnert to brick and save
            dem <- raster::brick(dem); names(dem) <- "dem"
            private$brk <- dem #raster::brick(private$projFiles["tif"])
            private$writeTIF()
        },
        ## add the channel
        apply_add_channel = function(chn){
            if( length(private$brk) == 0 ){
                stop("Project must be initialised with the DEM")
            }
            
            if( length(private$shp) > 0 ){
                stop("Channel already added - start a new project")
            }
            
            ## check the required field names are present
            if( !all( c("name","length","area","startNode","endNode") %in% names(chn)) ){
                stop("A required property name is not specified")
            }
            
            ## check projection of the chn
            if( !compareCRS(chn,private$brk) ){
                stop("Projection of channel object does not match that of project")
            }
            
            ## check if there is an id feild which will be overwritten
            if( ("id" %in% names(chn)) ){
                warning("The name id is reserved and will be overwritten")
            }

            ## ensure required properties are of correct type
            chn$name <- as.character(chn$name)
            chn$length <- as.numeric(chn$length)
            chn$area <- as.numeric(chn$area)
            chn$startNode <- as.character(chn$startNode)
            chn$endNode <- as.character(chn$endNode)
            if( !all(is.finite(chn$length)) ){
                stop("Some non-finite values of length found!")
            }
            
            ## arrange in id in order of flow direction - so lowest values at outlets of the network
            unds <- unique(c(chn$startNode,chn$endNode)) # unique nodes
            if(any(is.na(unds))){
                stop("The nodes in startNode and endNode must have unique, non-missing codes")
            }
            ## work out number of links from a node
            nFrom <- setNames(rep(0,length(unds)),unds)
            tmp <- table(chn$startNode)
            nFrom[names(tmp)] <- tmp
            max_id <- 0
            chn[["id"]] <- NA ## set id to NA
            while( any(nFrom==0) ){
                idx <- names(nFrom[nFrom==0])
                ## fill if of reaches which are at bottom
                tmp <- chn$endNode %in% idx
                chn$id[tmp] <- max_id + (1:sum(tmp))
                max_id <- max_id + sum(tmp)
                nFrom[idx] <- -1
                ## locate next nodes that are at bootom
                tmp <- table(chn$startNode[ tmp ])
                jj <- intersect(names(tmp),names(nFrom))
                nFrom[jj] <-  nFrom[jj] - tmp[jj]
            }

            ## create a raster of channel id numbers
            ## extract cells index, and fraction of river area in cell
            chn_rst <- private$brk[["dem"]]
            chn_rst[] <- NA
            ch_cell <- raster::extract(chn_rst,chn,weights=TRUE,
                                       cellnumbers=TRUE,na.rm=TRUE)
            ch_cell <- lapply(1:length(chn),function(ii){
                out <- as.data.frame(ch_cell[[ii]][,c("cell","weight")])
                out$weight <- out$weight * chn$area[ii]
                out$id <- chn$id[ii]
                out
            })
            ch_cell <- do.call(rbind,ch_cell)
            ch_cell <- split(ch_cell,ch_cell$cell)
            ch_cell <- lapply(ch_cell,function(x){x[which.max(x$weight),,drop=FALSE]})
            ch_cell <- do.call(rbind,ch_cell)
            chn_rst[ch_cell$cell] <- ch_cell$id

            private$brk[["channel"]] <- chn_rst
            private$shp <- chn
            private$writeTIF()
            raster::shapefile(chn,private$projFiles["shp"])            
        },
        ## Add a layer
        apply_add_layer=function(layer,layer_name){
            if( any(layer_name %in% private$reserved_layers) ){
                stop("Name is reserved")
            }
            if(!compareCRS(layer,private$brk)){ stop("Projection does not match project") }
            if(!all(res(layer)==res(private$brk))){ stop("Resolution does not match project") }
            if(!(extent(layer)==extent(private$brk))){
                ## try buffering it as for dem when read in
                layer <- extend(layer,c(1,1),NA)
            }
            if(!(extent(layer)==extent(private$brk))){ stop("Extent does not match project") }
            private$brk[[layer_name]] <- layer
            private$writeTIF()
            NULL
        },
        ## Sink fill
        apply_sink_fill = function(min_grad,max_it,verbose,hot_start){
            
            rq <- ifelse(hot_start,
                         c("filled_dem","channel"),
                         c("dem","channel"))
            if(!all(rq %in% names(private$brk))){
                stop("Not all required layers are available")
            }

            d <- ifelse(hot_start,"filled_dem","dem")
            d <- values( private$brk[[d]], format="matrix" )
            ch <- values( private$brk[["channel"]], format="matrix" )
           
            ## values that should be valid
            to_be_valid <- !is.na(d) & is.na(ch)  # all values not NA should have a valid height
            is_valid <- is.finite(ch) & !is.na(d) # TRUE if a channel cell for initialisation
            changed <- is_valid # cells changed at last iteration
            fd <- to_be_valid*Inf; fd[is_valid] <- d[is_valid]
            to_eval <- d; to_eval[] <- FALSE

            ## distance between cell centres
            rs <- raster::res( private$brk )
            dxy <- matrix(sqrt(sum(rs^2)),3,3)
            dxy[1,2] <- dxy[3,2] <- rs[2]
            dxy[2,1] <- dxy[2,3] <- rs[1]
            dxy[2,2] <- 0
            dxy <- min_grad*dxy
            
            it <- 1
            ## start of iteration loop
            while(any(changed[]) & it<=max_it){

                ## work out all the the cells to evaluate
                ## should be next to those that are changed
                to_eval[] <- FALSE
                idx <- which(changed,arr.ind=TRUE) # index of changed cells
                jdx <- idx
                for(ii in c(-1,0,1)){
                    for(jj in c(-1,0,1)){
                        if(ii==0 & jj==0){next}
                        ## adjust jdx
                        jdx[,1] <- idx[,1]+ii
                        jdx[,2] <- idx[,2]+jj
                        to_eval[jdx] <- !is_valid[jdx] #TRUE
                    }
                }
                to_eval <- to_eval & to_be_valid #& !is_valid
                if(verbose){
                    cat("Iteration",it,"\n")
                    cat("\t","Cells to evaluate:",sum(to_eval),"\n")
                    cat("\t","Percentage Complete:",
                        round(100*sum(is_valid)/sum(!is.na(d)),1),"\n") #to_be_valid),1),"\n")
                }
                ## alter min value for the evaluation cells
                idx <- which(to_eval,arr.ind=TRUE) # index of changed cells
                jdx <- idx
                mind <- rep(Inf,nrow(idx))
                for(ii in c(-1,0,1)){
                    for(jj in c(-1,0,1)){
                        if(ii==0 & jj==0){next}
                        ## adjust jdx
                        jdx[,1] <- idx[,1]+ii
                        jdx[,2] <- idx[,2]+jj
                        mind <- pmin(mind,fd[jdx] + dxy[ii+2,jj+2],na.rm=TRUE)
                    }
                }
                changed[] <- FALSE
                is_valid[idx] <- d[idx]>mind ## cells where the dem value is valid
                mind <- pmax(d[idx],mind) ## mind is now the replacemnt value
                changed[idx] <- mind < fd[idx]
                fd[idx] <- mind
                
                ## mind <- pmax(d[idx],mind)
                ## cdx <- mind < fd[idx]
                ## fd[idx] <- mind
                ## changed[] <- FALSE
                ## changed[idx[cdx,,drop=FALSE]] <- TRUE
                ## is_valid[idx] <- TRUE
                ## end of loop
                it <- it+1
                
            }

            rfd <- private$brk[["dem"]]; values(rfd) <- fd
            private$brk[["filled_dem"]] <- rfd
            
            private$writeTIF()
            
            if(it>max_it){ stop("Maximum number of iterations reached, sink filling not complete") }
            
        },
        ## Function to compute the properties
        apply_compute_properties = function(min_grad,verbose){

            rq <- c("filled_dem","channel")
            if(!all( rq %in% names( private$brk) )){
                stop("Not all required input layers have been generated \n",
                     "Try running sink_fill first")
            }
            
            ## load raster layer
            ## it is quickest to compute using blocks of dem as a raster
            ## however for small rasters we will just treat as a single block
            ## assumes the raster is padded with NA
            d <- values( private$brk[["filled_dem"]], format="matrix" )
            ch <- values( private$brk[["channel"]], format="matrix" )

            ## work out order to pass through the cells
            idx <- order(d,decreasing=TRUE,na.last=NA)
            n_to_eval <- length(idx)
            
            ## distances and contour lengths
            ## distance between cell centres
            rs <- raster::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
            ## initialise output
            gr <- upa <- atb <- d*NA
            upa[is.finite(d)] <- prod(rs) ## initialise upslope area from resolution

            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }else{
                next_print <- Inf
            }
            
            
            for(ii in idx){
                ngh <- ii + delta ## neighbouring cells
                ## compute gradient
                grd <- (d[ii]-d[ngh])/dxy
                
                to_use <- is.finite(grd) & (grd > 0)
                if(any(to_use)){
                    gcl <- grd[to_use]*dcl[to_use]
                    ## gradient
                    gr[ii] <- max(sum(gcl) / sum(dcl[to_use]),min_grad)
                    ## topographic index
                    atb[ii] <- log(upa[ii]/gr[ii]) #log( upa[ii] / sum(gcl) )
                    ## fraction of flow in each direction
                    frc <- gcl/sum(gcl)
                    ## propogate area downslope
                    upa[ ngh[to_use] ]  <- upa[ ngh[to_use] ] + frc*upa[ii]
                }else{
                    if( !is.finite(ch[ii]) ){
                        ## a hillslope cell that drains nowhere - this is an error
                        stop(paste("Cell",k,"is a hillslope cell with no lower neighbours"))
                    }
                }    
                
                ## verbose output here
                if(it >= next_print){
                    cat(round(100*it / n_to_eval,1),
                        "% complete","\n")
                    next_print <- next_print+print_step
                }
                
                it <- it+1
            }

            ## save raster maps
            out <- private$brk[["dem"]];
            values(out) <- gr; private$brk[["gradient"]] <- out
            values(out) <- upa; private$brk[["upslope_area"]] <- out
            values(out) <- atb; private$brk[["atb"]] <- out
            private$writeTIF()          
        },
        ## work out flow lengths to channel
        apply_flow_lengths = function(verbose){

            rq <- c("filled_dem","channel")
            if(!all( rq %in% names( private$brk) )){
                stop("Not all required input layers have been generated \n",
                     "Try running sink_fill first")
            }

            ## load raster layer
            ## it is quickest to compute using blocks of dem as a raster
            ## however for small rasters we will just treat as a single block
            ## assumes the raster is padded with NA
            d <- values( private$brk[["filled_dem"]], format="matrix" )
            ch <- values( private$brk[["channel"]], format="matrix" )
            
            ## create some distance matrices
            sfl <- d; sfl[] <- NA
            dfl <- d; dfl[] <- NA
            efl <- d; efl[] <- NA
            bnd <- d; bnd[] <- NA

            ## distances and contour lengths
            ## distance between cell centres
            rs <- raster::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
            ## if we go up in height order then we must have looked at all lower
            ## cells first
            idx <- order(d,na.last=NA)
            
            n_to_eval <- length(idx)

            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }else{
                next_print <- Inf
            }

            for(ii in idx){
                if(is.finite(ch[ii])){
                    ## just channel
                    bnd[ii] <- 0
                    sfl[ii]  <- dfl[ii] <- efl[ii] <- 0
                }else{
                    ## it is not a channel
                    jdx <- ii+delta
                    gcl <- (d[jdx]-d[ii])*dcl/dxy
                    is_lower <- is.finite(gcl) & gcl<0
                    kdx <- jdx[is_lower]
                    bnd[ii] <- max(bnd[kdx])+1
                    sfl[ii] <- min(sfl[kdx]+dxy[is_lower])
                    efl[ii] <- sum((efl[kdx]+dxy[is_lower])*gcl[is_lower])/sum(gcl[is_lower])
                    kdx <- which.min(gcl)
                    dfl[ii] <- dfl[jdx[kdx]]+dxy[kdx]
                }

                ## verbose output here
                if(it >= next_print){
                    cat(round(100*it / n_to_eval,1),
                        "% complete","\n")
                    next_print <- next_print+print_step
                }
                
                it <- it+1
            }


            ## save raster maps
            out <- private$brk[["dem"]];
            values(out) <- bnd; private$brk[["band"]] <- out
            values(out) <- sfl; private$brk[["shortest_flow_length"]] <- out
            values(out) <- dfl; private$brk[["dominant_flow_length"]] <- out
            values(out) <- efl; private$brk[["expected_flow_length"]] <- out
            private$writeTIF()
        },
        
        ## split_to_class
        apply_classify = function(layer_name,base_layer,cuts){

            ## check base layer exists
            if(!(base_layer %in% names(private$brk))){
                stop(paste(c("Missing layers:",base_layer,sep="\n")))
            }

            ## check layer_name isn't already used
            if(layer_name %in% names(private$brk)){
                stop("layer_name is already used")
            }

            ## load base layer and mask out channel
            x <-  mask( private$brk[[base_layer]], private$brk[["channel"]], inverse=TRUE)

            ## work out breaks
            brk <- as.numeric(cuts)
            if( length(brk)==1 | is.na(brk) ){
                ## this defines brks in the same way as cut would otherwise
                rng <- cellStats(x,range)
                brk <- seq(rng[1],rng[2],length=brk+1)
            }
            ## cut the raster
            x <- cut(x,breaks=brk,include.lowest=TRUE)

            private$brk[[layer_name]] <- x
            private$writeTIF()
            
            private$layers[[layer_name]] <- list(
                type="classification",
                layer=base_layer,
                cuts=brk)
            private$writeJSON()
        },
        ## split_to_class
        apply_combine_classes = function(layer_name,pairs,burns){
            
            ## check all cuts and burns are in possible layers
            rq <- c(pairs,burns)
            has_rq <- rq %in% names(private$brk)
            if(!all(has_rq)){
                stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
            }

            ## check layer_name isn't already used
            if(layer_name %in% names(private$brk)){
                stop("layer_name is already used")
            }
                
            ## work out new pairings by cantor method then renumber
            init <- TRUE
            for(ii in pairs){
                ## x <- private$brk[[ii]]
                x <-  mask( private$brk[[ii]], private$brk[["channel"]], inverse=TRUE)
                if(init){
                    cp <- x
                    init <- FALSE
                }else{
                    cp <- 0.5*(cp+x)*(cp+x+1)+x
                    uq <- sort(unique(cp))
                    cp <- cut(cp, c(-Inf,(uq[-1]+uq[-length(uq)])/2,Inf))
                }
            }

            ## put all the burns into a single raster
            brn <- cp; brn[] <- NA
            for(ii in burns){
                x <-  mask( private$brk[[ii]], private$brk[["channel"]], inverse=TRUE)
                if(is.null(brn)){
                    brn <- x
                }else{
                    brn <- cover(x,brn)
                }
            }
            ## add burns to pairs
            cp <- cover(brn,cp)
            uq <- sort(unique(cp))
            cts <- c(-Inf,(uq[-1]+uq[-length(uq)])/2,Inf)
            cp <- cut(cp,cts)
            if(length(burns)>0){ brn <- cut(brn,cts) }
            
            
            ## make table of layer values - should be able to combine with above??
            cpv <- raster::getValues(cp) ## quicker when a vector
            uq <- sort(unique(cp)) ## unique values
            uqb <- unique(brn) ## unique burn values
            
            cuq <- rep(NA,length(uq)) ##index of unique values
            for(ii in which(is.finite(cpv))){
                jj <- cpv[ii]
                if(is.na(cuq[jj])){ cuq[jj] <- ii }
            }
            
            if(!all(is.finite(cuq))){
                stop("Error in computing combinations")
            }
            
            ## create data frame
            df <- data.frame(uq); names(df) <- layer_name
            for(ii in pairs){
                df[[ii]] <- private$brk[[ii]][cuq] ## read in raster
            }
            df$burns <- df[[layer_name]] %in% uqb
            

            ## ## add in burns
            ## for(ii in burns){
            ##     x <-  mask( private$brk[[ii]], private$brk[["channel"]], inverse=TRUE)
            ##     ## x <- private$brk[[ii]] ## read in raster
            ##     idx <- Which(is.finite(x))
            ##     cp[idx] <- x[idx]

            ##     ux <- sort(unique(x))
            ##     ## ux that are alreasy layer numbers
            ##     idx <- df[,layer_name] %in% ux
            ##     df[idx,ii] <- df[idx,layer_name]
            ##     ## for ux not in df[,layer_name]
            ##     ux <- ux[!(ux %in% df[,layer_name])]
            ##     y <- matrix(NA,length(ux),ncol(df))
            ##     colnames(y) <- colnames(df)
            ##     y[,layer_name] <- y[,ii] <- ux
            ##     df <- rbind(df,y)
                
            ## }
            
            private$brk[[layer_name]] <- cp
            private$writeTIF()
            private$layers[[layer_name]] <- list(type="combination",
                                                 groups=df)
            private$writeJSON()
        },
        
        ## create a model  
        apply_create_model = function(class_lyr,dist_lyr,delta_dist,
                                      rain_lyr,rainfall_label,
                                      pet_lyr,pet_label,layer_name,verbose,
                                      transmissivity){

            rq <- c("atb","gradient","filled_dem","channel",
                    class_lyr,dist_lyr,
                    rain_lyr,pet_lyr)
            has_rq <- rq %in% names(private$brk)            
            if(!all(has_rq)){
                stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
            }
            if(!(class_lyr %in% names(private$layers) )){
                stop(class_lyr, " is not a classification")
            }

            model <- list()
            
            ## read in classification and distance layers
            cls <-  private$brk[[class_lyr]]## channel values are NA
            dst <- private$brk[[dist_lyr]]

            if(verbose){ cat("Initialising model","\n") }
            ## start model data frame and add minimum distance
            model$hru <- private$layers[[class_lyr]]$groups
            min_dst <- zonal(dst,cls,min) # minimum distance for each classification
            idx <- match(model$hru[[class_lyr]],min_dst[,"zone"])
            names(model$hru) <- paste0("cls_",names(model$hru))
            model$hru$min_dst <- min_dst[idx, "value"]
            ##browser()
            model$hru$zone <- min_dst[idx, "zone"]
            if(!all(is.finite(model$hru$min_dst))){
                stop("Unable to compute a finite minimum distance for all HRUs")
            }

            ## add id and change classification map to match
            model$hru <- model$hru[order(model$hru$min_dst,decreasing=TRUE),] # longest distance in first row
            model$hru$id <- (nrow(model$hru):1) + max(private$shp$id) ## ensure id is greater then number of channels
            map <- raster::subs(cls,model$hru,by="zone",which="id")

            ## add the channel data
            model$hru <- merge(model$hru,private$shp[,c("id","name","length")],by="id",all=TRUE)
            model$hru$is_channel <- model$hru$id <= max(private$shp$id)
            model$hru$min_dst[model$hru$is_channel] <- 0
            map <- cover(private$brk[["channel"]], map)


            
            
            ## add tranmissivity dependent parameters
            par <- switch(transmissivity,
                          "exp" = c(r_sfmax = Inf, c_sf = 0.1, s_rzmax = 0.05, t_d = 7200,
                                    ln_t0 = -2, c_sz = NA, m = 0.04, D= NA, m_2 = NA, omega = NA,
                                    s_rz0 = 0.75, r_uz_sz0 = 1e-6, s_raf=0.0, t_raf=Inf),
                          "bexp" = c(r_sfmax = Inf, c_sf = 0.1, s_rzmax = 0.05, t_d = 7200,
                                     ln_t0 = -2, c_sz = NA, m = 0.04, D=5, m_2 = NA, omega = NA,
                                     s_rz0 = 0.75, r_uz_sz0 = 1e-6, s_raf=0, t_raf=Inf),
                          "dexp" = c(r_sfmax = Inf, c_sf = 0.1, s_rzmax = 0.05, t_d = 7200,
                                     ln_t0 = -2, c_sz = NA, m = 0.04, D= NA, m_2 = 0.0002, omega = 0.5,
                                     s_rz0 = 0.75, r_uz_sz0 = 1e-6, s_raf=0, t_raf=Inf),
                          "cnst" = c(r_sfmax = Inf, c_sf = 0.1, s_rzmax = 0.05, t_d = 7200,
                                     ln_t0 = NA, c_sz = 0.01, m = NA, D= 10, m_2 = NA, omega = NA,
                                     s_rz0 = 0.75, r_uz_sz0 = 1e-6, s_raf=0.0, t_raf=Inf),
                          stop("Unrecognised tranmissivity")
                          )
            model$hru$opt = transmissivity
            model$hru$par <- rep(list(par),nrow(model$hru))

            
            ## work out the properties
            if(verbose){ cat("Computing properties","\n") }
            ## check sorted
            model$hru <- model$hru[order(model$hru$id),]
            if( !all( model$hru$id == 1:nrow(model$hru)) ){
                stop("id should be sequential from 1")
            }
            
            ## it is quickest to compute using blocks of a raster
            ## however for small rasters we will just treat as a single block
            d <- values( private$brk[["filled_dem"]] , format="matrix" )
            mp <- values( map , format="matrix")
            dst <- values( dst , format="matrix")
            gr <- values( private$brk[["gradient"]] , format="matrix")
            atb <- values( private$brk[["atb"]] , format="matrix")

            ## distances and contour lengths
            ## distance between cell centres
            rs <- raster::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(map); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
 
            ## initialise new variables
            model$hru$area <- model$hru$width <- model$hru$atb_bar <- model$hru$s_bar <- as.numeric(0)
            model$hru$width[model$hru$is_channel] <- NA

            flux <- rep( list(numeric(0)), nrow(model$hru))
            
            ## work out order to pass through the cells
            idx <- which(is.finite(mp))
            n_to_eval <- c(length(idx),0,length(idx)/10)

            for(ii in idx){
                #print(ii)
                id <- mp[ii] ## id
                model$hru$area[id] <- model$hru$area[id] + 1
                model$hru$atb_bar[id] <- model$hru$atb_bar[id] + atb[ii]
                model$hru$s_bar[id] <- model$hru$s_bar[id] + gr[ii]

                if( dst[ii] <= (model$hru$min_dst[id] + delta_dist) ){
                    
                    ## look for flow direction      
                    ngh <- ii + delta ## neighbouring cells
                    nghid <- mp[ngh]
                    ## compute gradient
                    grd <- (d[ii]-d[ngh])/dxy
                    ## test if can be used
                    to_use <- is.finite(grd) & (grd > 0) & is.finite(nghid) & (nghid < id)
                    
                    if( any(to_use) ){
                        nghid <- paste(nghid[to_use])
                        tmp <- nghid[!nghid %in% names(flux[[id]])]
                        if(length(tmp)>0){
                            flux[[id]][ tmp ] <- 0
                        }
                        flux[[id]][ nghid ] <- flux[[id]][ nghid ] + grd[to_use]*dcl[to_use]
                        model$hru$width[id] <- model$hru$width[id] + sum(dcl[to_use])
                    }
                }

                n_to_eval[2] <- n_to_eval[2] + 1
                if( verbose & (n_to_eval[2] > n_to_eval[3]) ){
                    cat(round(100*n_to_eval[3] / n_to_eval[1],1),
                        "% complete","\n")
                    n_to_eval[3] <- n_to_eval[3] + n_to_eval[1]/10
                }
                
            }
            model$hru$s_bar <- model$hru$s_bar / model$hru$area
            model$hru$atb_bar <- model$hru$atb_bar / model$hru$area
            model$hru$area <- model$hru$area * prod(rs)

            ## check fluxes are valid and convert to form sz_flux
            nlink <- sapply(flux,length)
            idx <- which(nlink==0)
            if( any(idx %in% model$hru$id[!model$hru$is_channel]) ){
                    stop( "The following HRUs have no valid outflows and are not channels: \n",
                         paste(idx[ idx %in% model$hru$id[!model$hru$is_channel] ], collapse=", "))
            }
            for(ii in 1:length(flux)){
                if(nlink[ii]>0){
                    flux[[ii]] <- data.frame(from = as.integer(ii),
                                             to = as.integer(names(flux[[ii]])),
                                             frc = flux[[ii]]/sum(flux[[ii]]))
                }
            }
            
            model$sz_flow_direction <- do.call(rbind,flux)

            ## alter channel flux for surface
            if(verbose){ cat("Creating channel flow directions","\n") }
            n_to_eval <- c(length(private$shp$id),0,length(private$shp$id)/10)
            for(ii in private$shp$id){
                idx <- model$channel$startNode == model$channel$endNode[ii]
                if(any(idx)){
                    flux[[ii]] <- data.frame(from = as.integer(ii),
                                             to = as.integer( private$shp$id[idx] ),
                                             frc = rep(1/sum(idx),sum(idx)))
                }else{
                    flux[[ii]] <- NULL
                }

                
                n_to_eval[2] <- n_to_eval[2] + 1
                if( verbose & (n_to_eval[2] > n_to_eval[3]) ){
                    cat(round(100*n_to_eval[3] / n_to_eval[1],1),
                        "% complete","\n")
                    n_to_eval[3] <- n_to_eval[3] + n_to_eval[1]/10
                }

            }

            ## create surface flux table
            model$sf_flow_direction <- do.call(rbind,flux)

            ## ############################################
            ## Add outlet at all outlets from river network
            ## ############################################
            idx <- model$hru$id[ !(model$hru$id %in% model$sf_flow_direction$from) ]
            model$output_flux<- data.frame(id = idx, flux = "q_sf", frc = 1)
            model$time_series <- data.frame(
                name = paste("channel",idx,sep="_"),
                id = idx,
                flux = "q_sf"
            )
            
            ## ##################################
            ## Add point inflow table
            ## ##################################
            ## blank point inflow table
            if(verbose){ cat("Adding a blank point_inflow table","\n") }
            model$input_flux<- data.frame(
                name = character(0),
                id = integer(0),
                state = character(0),
                frc = numeric(0)
            )
            
            ## ##################################
            ## Compute precip weights
            ## ##################################
            if(verbose){ cat("Computing rainfall weights","\n") }
            if( is.null(rain_lyr) ){
                model$precip_input <- data.frame(
                    name = "precip",
                    id = model$hru$id,
                    frc = as.numeric(1) )
            }else{
                tmp <- crosstab(map,private$brk[[rain_lyr]],long=TRUE)
                names(tmp[[ii]]) <- c("id","name","cnt")                
                tmp <- tmp[order(tmp$id),]
                tmp$id <- as.integer(tmp$id)
                tmp$name <- paste0(rainfall_label,tmp$name)
                ## tmp is ordered by id so the following returns correct order
                tmp$frc <- unlist(tapply(tmp$cnt,tmp$id,FUN=function(x){x/sum(x)}),use.names=FALSE)
                ## set
                model$precip_input <- tmp[,c("id","name","frc")]
            }
            
            ## #################################
            ## Comnpute PET weights
            ## #################################
            if(verbose){ cat("Computing the pet weights","\n") }
            if( is.null(pet_lyr) ){
                model$pet_input <- data.frame(
                    name = "pet",
                    id =  model$hru$id,
                    frc = as.numeric(1) )
            }else{
                tmp <- crosstab(map,private$brk[[pet_lyr]],long=TRUE)
                names(tmp[[ii]]) <- c("id","name","cnt")                
                tmp <- tmp[order(tmp$id),]
                tmp$id <- as.integer(tmp$id)
                tmp$name <- paste0(pet_label,tmp$name)
                ## tmp is ordered by id so the following returns correct order
                tmp$frc <- unlist(tapply(tmp$cnt,tmp$id,FUN=function(x){x/sum(x)}),use.names=FALSE)
                ## set
                model$pet_input <- tmp[,c("id","name","frc")]
            }
            
            ## ############################
            ## save model
            
            brk <- raster::brick(map); names(brk) <- "hru"
            brk[["class"]] <- private$brk[[class_lyr]]
            brk[["channel"]] <- private$brk[["channel"]]
            if( !is.null(pet_lyr) ){ brk[["pet"]] <- private$brk[[pet_lyr]] }
            if( !is.null(rain_lyr) ){ brk[["precip"]] <- private$brk[[rain_lyr]] }
            brk <- terra::rast(brk)
            terra::writeRaster(brk,paste0(layer_name,".tif"),overwrite=TRUE)
            saveRDS(model,paste0(layer_name,".rds"))
        }        
    )
    )
