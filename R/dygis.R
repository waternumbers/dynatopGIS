#' R6 Class for processing a catchment to make a Dynamic TOPMODEL
#' @examples
#' ## The vignettes contains more examples of the method calls.
#' 
#' ## create temport directory for output
#' demo_dir <- tempfile("dygis")
#' dir.create(demo_dir)
#'
#' ## initialise processing
#' ctch <- dynatopGIS$new(file.path(demo_dir,"test"))
#'
#' ## add digital elevation and channel data
#' dem_file <- system.file("extdata", "SwindaleDTM40m.tif", package="dynatopGIS", mustWork = TRUE)
#' dem <- terra::rast(dem_file)
#' ctch$add_dem(dem)
#' channel_file <- system.file("extdata", "SwindaleRiverNetwork.shp",
#' package="dynatopGIS", mustWork = TRUE)
#' sp_lines <- rgdal::readOGR(channel_file)
#' property_names <- c(name="identifier",endNode="endNode",startNode="startNode",length="length")
#' chn <- convert_channel(sp_lines,property_names)
#' ctch$add_channel(chn)
#'
#' ## compute properties 
#' ctch$sink_fill() ## fill sinks in the catchment
#' \donttest{
#' ctch$compute_properties() # like topograpihc index and contour length
#' ctch$compute_flow_lengths()
#' }
#' ## classify and create a model
#' \donttest{
#' ctch$classify("atb_20","atb",cuts=20) # classify using the topographic index
#' ctch$get_method("atb_20") ## see the details of the classification
#' ctch$combine_classes("atb_20_band",c("atb_20","band")) ## combine classes
#' ctch$create_model(file.path(demo_dir,"new_model"),"atb_20_band","band") ## create a model
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
        #' @param projectFolder folder for data files
        #'
        #' @details This loads the project data files found at \code{projFile} if present. If not the value is stored for later use. The project data files are given by \code{projfile,<tif,shp>}
        #'
        #' @return A new `dynatopGIS` object
        initialize = function(projectFolder){
            private$apply_initialize( projectFolder )
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
            if(!("SpatRaster" %in% class(dem))){ dem <- terra::rast(as.character(dem)) }
            if(!("SpatRaster" %in% class(dem))){ stop("dem is not a SpatRaster") }
            private$apply_add_dem(dem,fill_na)
            invisible(self)
        },
        #' @description Import channel data to the `dynatopGIS` object
        #'
        #' @param channel a SpatialLinesDataFrame, SpatialPolygonsDataFrame or file path containing the channel information
        #'
        #' @details Takes the representation of the channel network as a SpatialPolygonsDataFrame with properties name, length, area, startNode, endNode and overlaying it on the DEM. In doing this a variable called id is created (or overwritten) other variables in the data frame are passed through unaltered.
        #'
        #' @return suitable for chaining
        add_channel = function(channel){
            if(!is(channel,"SpatVector")){ channel <- terra::vect( as.character(channel) ) }
            if(!is(channel,"SpatVector")){ stop("channel is not a SpatialPolygonsDataFrame") }

            private$apply_add_channel(channel)
            invisible(self)
        },
        #' @description Add a layer of geographical information
        #'
        #' @param layer the raster layer to add (see details)
        #' @param layer_name name to give to the layer
        #'
        #' @details The layer should either be a raster layer or a file that can be read by the \code{raster} package. The projection, resolution and extent are checked against the existing project data. Only layer names not already in use (or reserved) are allowed. If successful the layer is added to the project tif file.
        #' @return suitable for chaining
        add_layer = function(layer,layer_name=names(layer)){
            layer_name <- as.character(layer_name)
            if(!("SpatRaster" %in% class(layer))){ layer <- terra::rast(as.character(layer)) }
            if(!("SpatRaster" %in% class(layer))){ stop("layer is not a SpatRaster") }
            if( length(layer_name) != terra::nlyr(layer) ){ stop("Length of names does not match number fo layers") }
            private$apply_add_layer(layer,layer_name)
            invisible(self)
        },
        #' @description Get a layer of geographical information or a list of layer names
        #' @param layer_name name of the layer give to the layer
        #' @return a `raster` layer of the requested information if layer_name is given else a vector of layer names
        get_layer = function(layer_name=character(0)){

            ## handle case where a list of layers is requested
            if( length(layer_name) == 0 ){
                
                return( data.frame(layer = names(private$brk), source=terra::sources(private$brk)) )
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
            stop("TODO - not working")
            if(length(layer_name)==0){
 
                return( private$layers )
            }
            layer_name <- match.arg(layer_name,names(private$layers))
            return( private$layers[[layer_name]] )
        }               
    ),
    private = list(
        version = "0.3.0",
        projectFolder = character(0),
        brk = character(0),
        shp = character(0),
        reserved_layers = c("dem","channel","filled_dem",
                            "gradient","upslope_area","atb",
                            "band","shortest_flow_length",
                            "dominant_flow_length","expected_flow_length"),
        readJSON = function(fn){
            
            if(!file.exists(fn)){
                return(NULL) #warning("No method json file found")
            }
            jsonlite::fromJSON(fn)
        },
        ## check and read the project files
        apply_initialize = function(projectFolder){
            ## see if folder exists
            if( !dir.exists(projectFolder) ){
                print( "Creating new folder" )
                dir.create(projectFolder, showWarnings = TRUE, recursive = TRUE)
            }

            ## read in raster data
            rstNames <- list.files(projectFolder, pattern=".tif$",full.names=TRUE)
            if(length(rstNames) == 0){
                brk <- chn <- NULL
                print( paste("Starting new project at", projectFolder) )
            }else{
                print( paste("Reading existing project at", projectFolder) )
                brk <- terra::rast(list.files(projectFolder, pattern=".tif$",full.names=TRUE))

                if( !all.equal(diff(terra::res(brk)),0) | terra::is.lonlat(brk) ){
                    stop("Raster data not valid: Processing currently only works on projected data with a square grid")
                }
                if( !("dem" %in% names(brk)) ){ stop("No dem layer in Raster data") }
                ## read in shape file and check
                if( file.exists( file.path(projectFolder,"channel.shp")) ){
                    chn <- terra::vect( file.path(projectFolder,"channel.shp") )
                
                    if(  terra::crs(brk, proj=TRUE) != terra::crs(chn, proj=TRUE) ){
                        stop("Shape and Raster projections do not match")
                    }
                    if( !("channel" %in% names(brk)) ){
                        stop("No channel layer in the raster file but channel shapefile specified")
                    }
 
                }else{
                    chn <- NULL
                     if( ("channel" %in% names(brk)) ){
                        stop("Channel layer in the raster file but channel shapefile missing")
                     }
                }
            }

            private$projectFolder <- projectFolder
            private$brk <- brk
            private$shp <- chn
        },
        ## adding dem
        apply_add_dem = function(dem,fill_na){
           
            if("dem" %in% names(private$brk) ){
                stop("The DEM exists, start a new project")
            }
            
            ## add to ensure NA on each edge
            dem <- terra::extend(dem,c(1,1)) # ,NA)
            ## check the projection is OK
            if( !all.equal(diff(terra::res(dem)),0) | terra::is.lonlat(dem) ){
                stop("Processing currently only works on projected dem's with a square grid")
            }
            
            ## fill dem so only NA values are on boundaries
            if(fill_na){
                ## so all catchment dem value a number
                na_clumps <- terra::patches(is.na(dem), directions=8, zeroAsNA=TRUE)
                edge_values <- unique( c(
                    unique( terra::values(na_clumps, mat=TRUE, dataframe=FALSE, row=1,nrows=nrow(na_clumps), col=1, ncols=1, na.rm=TRUE) ),
                    unique( terra::values(na_clumps, mat=TRUE, dataframe=FALSE, row=1,nrows=nrow(na_clumps), col=ncol(na_clumps), ncols=1, na.rm=TRUE) ),
                    unique( terra::values(na_clumps, mat=TRUE, dataframe=FALSE, row=1,nrows=1, col=1, ncols=ncol(na_clumps), na.rm=TRUE) ),
                    unique( terra::values(na_clumps, mat=TRUE, dataframe=FALSE, row=nrow(na_clumps),nrows=1, col=1, ncols=ncol(na_clumps), na.rm=TRUE) )
                ))
                na_clumps <- terra::subst(na_clumps, edge_values, NA) ## clean out edge values since these are not sinks
                sink_values <- terra::unique(na_clumps)$patches
                if(length(sink_values)>0){
                    na_clumps <- terra::subst(na_clumps, unique(na_clumps)$patches, 1 - 1e6) ## clean out edge values since these are not sinks
                }
                dem <- terra::cover(dem,na_clumps) # replace
            }
            ## save
            names(dem) <- "dem"
            demFile <- file.path(private$projectFolder,"dem.tif")
            terra::writeRaster(dem, demFile)
            private$brk <- terra::rast(demFile)
            
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
            if( terra::crs(private$brk, proj=TRUE) != terra::crs(chn, proj=TRUE) ){
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
            ## work out number of links from a node - TODO this is slow.........
            nFrom <- setNames(rep(0,length(unds)),unds)
            tmp <- table(chn$startNode)
            nFrom[names(tmp)] <- tmp
            max_id <- 0
            chn$id <- as.numeric(NA) ## set id to NA
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
            chn_rst <- terra::rasterize(chn,private$brk[["dem"]],field = "id",touches=TRUE)
            names(chn_rst) <- "channel"
            
            shpFile <- file.path(private$projectFolder,"channel.shp")
            rstFile <- file.path(private$projectFolder,"channel.tif")
            terra::writeVector(chn, shpFile)
            terra::writeRaster(chn_rst,rstFile, names = "channel")
            
            private$brk[["channel"]] <- terra::rast(rstFile)
            private$shp <- terra::vect(shpFile)
        },
        ## Add a layer
        apply_add_layer=function(layer,layer_name){
            
            if( any(layer_name %in% private$reserved_layers) ){
                stop("Name is reserved")
            }
            if( any(layer_name %in% names(private$brk)) ){
                stop("Name is already used")
            }
            if( !terra::compareGeom(layer,private$brk,stopOnError=FALSE) ){
                ## try buffering it as for dem when read in
                layer <- terra::extend(layer,c(1,1))
            }
            if( !terra::compareGeom(layer,private$brk,stopOnError=FALSE) ){
                stop("New layer does not match resolution, extent or projection of project")
            }
            names(layer) <- layer_name
            rstFile <- file.path(private$projectFolder,paste0(layer_name[1],".tif"))
            terra::writeRaster(layer, rstFile)
            private$brk[[layer_name]] <- terra::rast( rstFile )
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
            d <- terra::as.matrix( private$brk[[d]] ,wide=TRUE) ##, mat=TRUE)
            ch <-terra::as.matrix( private$brk[["channel"]] , wide=TRUE )
           
            ## values that should be valid
            to_be_valid <- !is.na(d) & is.na(ch)  # all values not NA should have a valid height
            is_valid <- is.finite(ch) & !is.na(d) # TRUE if a channel cell for initialisation
            changed <- is_valid # cells changed at last iteration
            fd <- to_be_valid*Inf; fd[is_valid] <- d[is_valid]
            to_eval <- d; to_eval[] <- FALSE

            ## distance between cell centres
            rs <- terra::res( private$brk )
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

            
            rfd <- terra::rast( private$brk[["dem"]], names="filled_dem", vals=fd )
            rstFile <- file.path(private$projectFolder,"filled_dem.tif")
            terra::writeRaster(rfd, rstFile)
            private$brk[["filled_dem"]] <- terra::rast(rstFile)
            
           
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
            
            d <- terra::as.matrix( private$brk[["filled_dem"]] , wide=TRUE)
            ch <- terra::as.matrix( private$brk[["channel"]] , wide=TRUE)

            ## work out order to pass through the cells
            idx <- order(d,decreasing=TRUE,na.last=NA)
            n_to_eval <- length(idx)
            
            ## distances and contour lengths
            ## distance between cell centres
            rs <- terra::res( private$brk )
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
            out <- terra::rast( private$brk[["dem"]], names="gradient", vals=gr )
            rstFile <- file.path(private$projectFolder,"gradient.tif")
            terra::writeRaster(out, rstFile); private$brk[["gradient"]] <- terra::rast(rstFile)

            out <- terra::rast( private$brk[["dem"]], names="upslope_area", vals=upa )
            rstFile <- file.path(private$projectFolder,"upslope_area.tif")
            terra::writeRaster(out, rstFile); private$brk[["upslope_area"]] <- terra::rast(rstFile)

            out <- terra::rast( private$brk[["dem"]], names="atb", vals=atb )
            rstFile <- file.path(private$projectFolder,"atb.tif")
            terra::writeRaster(out, rstFile); private$brk[["atb"]] <- terra::rast(rstFile)

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
            d <- terra::as.matrix( private$brk[["filled_dem"]], wide=TRUE )
            ch <- terra::as.matrix( private$brk[["channel"]],  wide=TRUE )
            
            ## create some distance matrices
            sfl <- d; sfl[] <- NA
            dfl <- d; dfl[] <- NA
            efl <- d; efl[] <- NA
            bnd <- d; bnd[] <- NA

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
            out <- terra::rast( private$brk[["dem"]], names="band", vals=bnd )
            rstFile <- file.path(private$projectFolder,"band.tif")
            terra::writeRaster(out, rstFile); private$brk[["band"]] <- terra::rast(rstFile)

            out <- terra::rast( private$brk[["dem"]], names="shortest_flow_length", vals=sfl )
            rstFile <- file.path(private$projectFolder,"shortest_flow_length.tif")
            terra::writeRaster(out, rstFile); private$brk[["shortest_flow_length"]] <- terra::rast(rstFile)

            out <- terra::rast( private$brk[["dem"]], names="dominant_flow_length", vals=dfl )
            rstFile <- file.path(private$projectFolder,"dominant_flow_length.tif")
            terra::writeRaster(out, rstFile); private$brk[["dominant_flow_length"]] <- terra::rast(rstFile)

            out <- terra::rast( private$brk[["dem"]], names="expected_flow_length", vals=efl )
            rstFile <- file.path(private$projectFolder,"expected_flow_length.tif")
            terra::writeRaster(out, rstFile); private$brk[["expected_flow_length"]] <- terra::rast(rstFile)

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

            ## TODO - check layer_name isn;t reserved
            
            ## load base layer and mask out channel
            x <-  terra::mask( private$brk[[base_layer]], private$brk[["channel"]], inverse=TRUE)

            ## work out breaks
            brk <- as.numeric(cuts)
            if( length(brk)==1 | is.na(brk) ){
                ## this defines brks in the same way as cut would otherwise
                rng <- as.numeric( terra::global(x, fun="range",na.rm=TRUE) )
                brk <- seq(rng[1],rng[2],length=brk+1)
            }
            ## cut the raster
            outFile <- file.path(private$projectFolder,paste0(layer_name,".tif"))
            private$brk[[ layer_name ]] <- terra::classify(x,rcl=brk,include.lowest=TRUE,
                                 filename= outFile, names=layer_name)

            out <- list(type="classification", layer=base_layer, cuts=brk)
            writeLines( jsonlite::toJSON(out), file.path(private$projectFolder,paste0(layer_name,".json")) )
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
                x <-  terra::mask( private$brk[[ii]], private$brk[["channel"]], inverse=TRUE)
                if(init){
                    cp <- x
                    init <- FALSE
                }else{
                    cp <- 0.5*(cp+x)*(cp+x+1)+x
                    uq <- sort(terra::unique(cp)[[1]])
                    cp <- terra::classify(cp, c(-Inf,(uq[-1]+uq[-length(uq)])/2,Inf),include.lowest=TRUE)
                }
            }

            ## put all the burns into a single raster
            brn <- cp; brn[] <- NA
            for(ii in burns){
                x <-  terra::mask( private$brk[[ii]], private$brk[["channel"]], inverse=TRUE)
                if(is.null(brn)){
                    brn <- x
                }else{
                    brn <- terra::cover(x,brn)
                }
            }
            ## add burns to pairs
            cp <- terra::cover(brn,cp)
            uq <- sort(terra::unique(cp)[[1]])
            cts <- c(-Inf,(uq[-1]+uq[-length(uq)])/2,Inf)
            cp <- terra::classify(cp,cts)
            if(length(burns)>0){ brn <- terra::classify(brn,cts) }
            
            cp <- cp + 1 ## classify returns numeric values starting at 0
            
            ## make table of layer values - should be able to combine with above??
            cpv <- terra::as.matrix(cp, wide=TRUE) ## quicker when a vector
            uq <- sort(terra::unique(cp)[[1]]) ## unique values
            uqb <- terra::unique(brn)[[1]] ## unique burn values
            
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
                df[[ii]] <- terra::as.matrix(private$brk[[ii]],wide=TRUE)[cuq] ## read in raster
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
            outFile <- file.path(private$projectFolder,paste0(layer_name,".tif"))
            private$brk[[ layer_name ]] <- terra::writeRaster( cp, outFile, names=layer_name)
            
            out <- list(type="combination",
                        groups=df)
            writeLines( jsonlite::toJSON(out), file.path(private$projectFolder,paste0(layer_name,".json")) )

 
            ## private$brk[[layer_name]] <- cp
            ## private$writeTIF()
            ## private$layers[[layer_name]] <- list(type="combination",
            ##                                      groups=df)
            ## private$writeJSON()
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

            jsonFile <- paste0(tools::file_path_sans_ext(terra::sources(private$brk[[class_lyr]])),".json")
            if( !file.exists(jsonFile) ){
                stop("No json file giving basis of the classifications")
            }
            
            ## if(!(class_lyr %in% names(private$layers) )){
            ##     stop(class_lyr, " is not a classification")
            ## }

            model <- list()
            
            ## read in classification and distance layers
            cls <-  private$brk[[class_lyr]]## channel values are NA
            dst <- private$brk[[dist_lyr]]

            if(verbose){ cat("Initialising model","\n") }
            ## start model data frame and add minimum distance
            model$hru <- jsonlite::fromJSON(jsonFile)$groups
            
            min_dst <- terra::zonal(dst,cls,min) # minimum distance for each classification
            idx <- match(model$hru[[class_lyr]],min_dst[[class_lyr]])
            names(model$hru) <- paste0("cls_",names(model$hru))
            model$hru$min_dst <- min_dst[idx, dist_lyr]
            ##browser()
            ##model$hru$zone <- min_dst[idx, "zone"]
            if(!all(is.finite(model$hru$min_dst))){
                stop("Unable to compute a finite minimum distance for all HRUs")
            }

            ## add id and change classification map to match
            model$hru <- model$hru[order(model$hru$min_dst,decreasing=TRUE),] # longest distance in first row
            model$hru$id <- (nrow(model$hru):1) + max(private$shp$id) ## ensure id is greater then number of channels
            map <- terra::subst(cls,from=model$hru[[paste0("cls_",class_lyr)]],to=model$hru$id)

            ## add the channel data
            model$hru <- merge(model$hru,private$shp[,c("id","name","length")],by="id",all=TRUE)
            model$hru$is_channel <- model$hru$id <= max(private$shp$id)
            model$hru$min_dst[model$hru$is_channel] <- 0
            map <- terra::cover(private$brk[["channel"]], map)

            ## add the options controlling lateral flow calculations
            model$hru$q_sf <- "constant_velocity"
            model$hru$q_sz <- transmissivity

            
                      
            ## add parameters
            fqsf <- function(x){
                list(opt = x,
                     par = switch(x,
                                  "constant_velocity" = list(r_sfmax = Inf, v_sf = 0.1),
                                  "constant_velocity_raf" = list(r_sfmax = Inf, v_sf = 0.1,s_raf=100,t_raf=10*60*60),
                                  stop("Unrecognised surface option")
                                  )
                     )
            }
            
            fqsz <- function(x){
                list(opt = x,
                     par = switch(x,
                                  "exp" = list( ln_t0 = -2, m = 0.04 ),
                                  "bexp" = list( ln_t0 = -2, m = 0.04, D=5 ),
                                  "dexp" = list( ln_t0 = -2, m = 0.04, m_2 = 0.0002, omega = 0.5 ),
                                  "cnst" = list( v_sz = 0.01, D= 10 ),
                                  stop("Unrecognised tranmissivity")
                                  )
                     )
            }
            
            model$hru$q_sf <- lapply( model$hru$q_sf,fqsf )
            model$hru$q_sz <- lapply( model$hru$q_sz,fqsz )
            model$hru$s_rzmax <- 0.05
            model$hrus_rz0 <- 0.75
            model$hru$t_d <- 7200
            model$hru$r_uz_sz0 <- 1e-6
            
            ##model$hru$par <- Map(c, sf_par, rz_uz_par, sz_par)

            
            ## work out the properties
            if(verbose){ cat("Computing properties","\n") }
            ## check sorted
            model$hru <- model$hru[order(model$hru$id),]
            if( !all( model$hru$id == 1:nrow(model$hru)) ){
                stop("id should be sequential from 1")
            }
            
            ## it is quickest to compute using blocks of a raster
            ## however for small rasters we will just treat as a single block
            d <- terra::as.matrix( private$brk[["filled_dem"]] , wide=TRUE )
            mp <- terra::as.matrix( map  , wide=TRUE )
            dst <- terra::as.matrix( dst  , wide=TRUE )
            gr <- terra::as.matrix( private$brk[["gradient"]]  , wide=TRUE )
            atb <- terra::as.matrix( private$brk[["atb"]]  , wide=TRUE )

            ## distances and contour lengths
            ## distance between cell centres
            rs <- terra::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(map); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
 
            ## initialise new variables
            model$hru$area <- model$hru$width <- model$hru$atb_bar <- model$hru$s_bar <- as.numeric(0)
            
            flux <- rep( list(numeric(0)), nrow(model$hru))
            
            ## work out order to pass through the cells
            idx <- which(is.finite(mp))
            n_to_eval <- c(length(idx),0,length(idx)/10)

            for(ii in idx){
                
                id <- mp[ii] ## id
                model$hru$area[id] <- model$hru$area[id] + 1
                model$hru$atb_bar[id] <- model$hru$atb_bar[id] + atb[ii]
                model$hru$s_bar[id] <- model$hru$s_bar[id] + gr[ii]
                
                if( model$hru$is_channel[id]  |  (dst[ii] <= (model$hru$min_dst[id] + delta_dist)) ){
                    
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
            if(verbose){ cat("Adding a blank inflow table","\n") }
            model$inflow_table <- data.frame(
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
                browser()
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
            
            names(map) <- "hru"
            map[["class"]] <- private$brk[[class_lyr]]
            map[["channel"]] <- private$brk[["channel"]]
            if( !is.null(pet_lyr) ){ brk[["pet"]] <- private$brk[[pet_lyr]] }
            if( !is.null(rain_lyr) ){ brk[["precip"]] <- private$brk[[rain_lyr]] }
            terra::writeRaster(map,paste0(layer_name,".tif"),overwrite=TRUE)
            saveRDS(model,paste0(layer_name,".rds"))
        }        
    )
    )
