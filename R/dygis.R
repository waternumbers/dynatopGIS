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
#' sp_lines <- terra::vect(channel_file)
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
        #' @description Import a dem to the `dynatopGIS` object
        #'
        #' @param dem a \code{raster} layer object or the path to file containing one which is the DEM
        #' @param fill_na  should NA values in dem be filled. See details
        #' @param verbose Should additional progress information be printed
        #'
        #' @details If not a \code{raster} the DEM is read in using the terra package. If \code{fill_na} is \code{TRUE} all NA values other then those that link to the edge of the dem are filled so they can be identified as sinks.
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
            if(!is(channel,"SpatVector")){ stop("channel is not a SpatVector object") }

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
            terra::plot( lyr, main = layer_name)
            if( add_channel & length(private$shp) > 0){
                terra::plot(private$shp, add=TRUE )
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
        #' @details TODO rewrite - The algorithm passed through the cells in increasing height. For measures of flow length to the channel are computed. The shortest length (minimum length to channel through any flow path), the dominant length (the length taking the flow direction with the highest fraction for each pixel on the path) and expected flow length (flow length based on sum of downslope flow lengths based on fraction of flow to each cell) and band (strict sequence to ensure that all contributing cell have a higher band value). By definition cells in the channel that have no land area have a length (or band) of NA.
        compute_flow_lengths = function(flow_routing=c("expected","dominant","d8","shortest"), verbose=FALSE){
            flow_routing = match.arg(flow_routing)
            private$apply_flow_lengths(flow_routing,verbose)
            invisible(self)
        }, 
        #' @description Create a catchment classification based cutting an existing layer into classes
        #' @param layer_name name of the new layer to create
        #' @param base_layer name of the layer to be cut into classes
        #' @param cuts values on which to cut into classes. These should be numeric and define either the number of bands (single value) or breaks between band (multiple values).
        #'
        #' @details This applies the given cuts to the supplied landscape layer to produce areal groupings of the catchment. Cuts are implement using \code{terra::cut} with \code{include.lowest = TRUE}. Note that is specifying a vector of cuts values outside the limits will be set to NA.
        classify = function(layer_name,base_layer,cuts){
            private$apply_classify(as.character(layer_name), as.character(base_layer), cuts)
            invisible(self)
        },
        #' @description Combine any number of classifications based on unique combinations and burns
        #' @param layer_name name of the new layer to create
        #' @param pairs a vector of layer names to combine into new classes through unique combinations. Names should correspond to raster layers in the project directory.
        #' @param burns a vector of layer names which are to be burnt on
        #'
        #' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. Burns are added directly in the order they are given. Cuts are implement using \code{terra::cut} with \code{include.lowest = TRUE}. Note that is specifying a vector of cuts values outside the limits will be set to NA.
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
        #' @param sf_opt Surface solution to use
        #' @param sz_opt transmissivity transmissivity profile to use
        #' @param dist_delta TODO
        #' @param rain_layer the layer defining the rainfall inputs
        #' @param rain_label Prepended to rain_layer values to give rainfall series name        
        #' @param pet_layer the layer defining the pet inputs
        #' @param pet_label Prepended to pet_layer values to give pet series name
        #' @param verbose print more details of progress
        #'
        #' @details The \code{class_layer} is used to define the HRUs. Flow between HRUs is based on the distance to a channel. For each HRU the shortest distance to a channel is computed. Flow from a HRU can only go to a HRU with a lower shortest distance to the channel. Flow from a HRU can occur from any raster cell within the HRU whose distance to the channel is within dist_delta of the shortest distance within the HRU.
        #' Setting the sf_opt and sz_opt options ensures the model is set up with the correct parameters present.
        #' The \code{rain_layer} (\code{pet_layer}) can contain the numeric id values of different rainfall (pet) series. If the value of \code{rain_layer} (\code{pet_layer}) is not \code{NULL} the weights used to compute an averaged input value for each HRU are computed, otherwise an input table for the models generated with the value "missing" used in place of the series name.
        create_model = function(layer_name,class_layer,dist_layer,
                                sf_opt = c("cnst","kin"),
                                sz_opt = c("exp","bexp","cnst","dexp"),
                                dist_delta=0,
                                rain_layer=NULL, rain_label=character(0),
                                pet_layer=NULL, pet_label=character(0),
                                verbose=FALSE){
            
            ## check valid transmissivity and channel_solver
            sf_opt<- match.arg(sf_opt)
            sz_opt <- match.arg(sz_opt)
            
            private$apply_create_model(class_layer, dist_layer,dist_delta,
                                       rain_layer, rain_label,
                                       pet_layer, pet_label,
                                       layer_name,verbose,
                                       sf_opt, sz_opt)

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
        get_method = function(layer_name){
            ## check layer name exists
            layer_name <- match.arg(layer_name,names(private$brk))
            
            jsonFile <- paste0(tools::file_path_sans_ext(terra::sources(private$brk[[layer_name]])),".json")
            if( !file.exists(jsonFile) ){
                stop("No json file giving basis of the classifications")
            }

            return( jsonlite::fromJSON( jsonFile ) )
        }               
    ),
    private = list(
        version = "0.3.0",
        projectFolder = character(0),
        brk = character(0),
        shp = character(0),
        reserved_layers = c("dem","channel","filled_dem",
                            "gradient","upslope_area","atb",
                            "band","flow_length"),
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
                message( "Creating new folder" )
                dir.create(projectFolder, showWarnings = TRUE, recursive = TRUE)
            }

            ## read in raster data
            rstNames <- list.files(projectFolder, pattern=".tif$",full.names=TRUE)
            if(length(rstNames) == 0){
                brk <- chn <- NULL
                message( paste("Starting new project at", projectFolder) )
            }else{
                message( paste("Reading existing project at", projectFolder) )
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
            
            ## TODO
            ## - add check that flow_direction.rds exists if flow direct processed
            ## - add check that channel_direction.rds exists if channel processed
            
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
                    terra::unique( terra::values(na_clumps, mat=TRUE, dataframe=FALSE, row=1,nrows=nrow(na_clumps), col=1, ncols=1, na.rm=TRUE) ),
                    terra::unique( terra::values(na_clumps, mat=TRUE, dataframe=FALSE, row=1,nrows=nrow(na_clumps), col=ncol(na_clumps), ncols=1, na.rm=TRUE) ),
                    terra::unique( terra::values(na_clumps, mat=TRUE, dataframe=FALSE, row=1,nrows=1, col=1, ncols=ncol(na_clumps), na.rm=TRUE) ),
                    terra::unique( terra::values(na_clumps, mat=TRUE, dataframe=FALSE, row=nrow(na_clumps),nrows=1, col=1, ncols=ncol(na_clumps), na.rm=TRUE) )
                ))
                na_clumps <- terra::subst(na_clumps, edge_values, NA) ## clean out edge values since these are not sinks
                sink_values <- terra::unique(na_clumps,na.rm=TRUE)$patches
                if(length(sink_values)>0){
                    na_clumps <- terra::subst(na_clumps, sink_values, 1 - 1e6) ## clean out edge values since these are not sinks
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
            if( !all( c("name","length","area","startNode","endNode","width","slope") %in% names(chn)) ){
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

            ## arrange id in order of flow direction - so lowest values at outlets of the network
            ## This is much quicker as a vector is not constantly accessing via the vect object
            id <- rep(as.integer(0),nrow(chn))
            sN <- chn$startNode
            eN <- chn$endNode
            idx <- !(eN %in%sN) ## outlets are channel lengths whose outlet does not join another channel
            it <- 1
            while(sum(idx)>0){
                id[idx] <- max(id) + 1:sum(idx)
                idx <- eN %in% sN[idx]
                it <- it+1
            }
            chn$id <- id
            if( !(all(chn$id>0))){ stop("error ingesting channel") }
            chn <- chn[ order(chn$id),]
            ##chn$id <- 1:nrow(chn)
                        
            ## create a raster of channel id numbers
            ## TODO - possibly sort so do biggest area first???
            chn_rst <- terra::rasterize(chn,private$brk[["dem"]],field = "id",touches=TRUE)
            chn_rst <- terra::mask(chn_rst,private$brk[["dem"]]) ## make sure value occur only on cells with height
            names(chn_rst) <- "channel"

            
            ## adjust channel to only include those in the raster
            idx <- terra::unique(chn_rst)$channel
            chn$to_keep <- chn$id %in% idx
            sN <- chn$startNode
            eN <- chn$endNode
            for(ii in which(!chn$to_keep)){
                eN[ eN==sN[ii] ] <- eN[ii]
            }            
            chn$endNode <- eN
            chn$startNode <- sN
            chn <- chn[chn$to_keep,]

            ## renumber the channel id
            ## chn_rst <- terra::subst(chn_rst, chn$id, 1:nrow(chn)) - this is very slow
            chn$id <- 1:nrow(chn)
            ## this is quicker but possibly prone to error??
            ##tmp <- chn_rst
            chn_rst <- terra::rasterize(chn,private$brk[["dem"]],field = "id",touches=TRUE)
            chn_rst <- terra::mask(chn_rst,private$brk[["dem"]]) ## make sure value occur only on cells with height
            names(chn_rst) <- "channel"
            
            shpFile <- file.path(private$projectFolder,"channel.shp")
            rstFile <- file.path(private$projectFolder,"channel.tif")
            terra::writeVector(chn, shpFile)
            terra::writeRaster(chn_rst,rstFile, names = "channel")
            
            private$brk <- c(private$brk, terra::rast(rstFile))
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
            rstFile <- file.path(private$projectFolder,paste0(layer_name,".tif"))
            terra::writeRaster(layer, rstFile)
            private$brk <- c(private$brk, terra::rast( rstFile ))
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
            if(hot_start){
                private$brk[["filled_dem"]] <- terra::rast(rstFile)
            }else{
                private$brk <- c( private$brk, terra::rast(rstFile))
            }
            
            
           
            if(it>max_it){ stop("Maximum number of iterations reached, sink filling not complete") }
            
        },
        ## Function to compute the properties
        apply_compute_properties = function(min_grad,verbose){

            rq <- c("filled_dem","channel")
            if(!all( rq %in% names( private$brk) )){
                stop("Not all required input layers have been generated \n",
                     "Try running sink_fill first")
            }

            rq <- c(file.path(private$projectFolder,"flow_direction.rds"),
                    file.path(private$projectFolder,"channel_routing.rds"))
            if( ! all( file.exists(rq) ) ){
                stop("No flow routing method defined\n",
                     "Try running compute_flow_lengths first")
            }else{
                flow_routing <- readRDS(rq[1])
                channel_routing <- readRDS( rq[2] )
            }
            
            ## load raster layer
            ## it is quickest to compute using blocks of dem as a raster
            ## however for small rasters we will just treat as a single block
            ## assumes the raster is padded with NA
            
            d <- terra::as.matrix( private$brk[["filled_dem"]] , wide=TRUE)
            ch <- terra::as.matrix( private$brk[["channel"]] , wide=TRUE)
            
            ## work out order to pass through the cells
##            idx <- order(d,decreasing=TRUE,na.last=NA)
##            n_to_eval <- length(idx)
            n_to_eval <- nrow(flow_routing)
            
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
            
            
            for(rwnum in nrow(flow_routing):1){
                ii <- flow_routing[rwnum,1]
                w <- flow_routing[rwnum,2:9]
                
                ngh <- ii + delta ## neighbouring cells
                ## compute gradient
                grd <- (d[ii]-d[ngh])/dxy

                to_use <- w>0 ##is.finite(grd) & (grd > 0)
                if(any(to_use) & is.na(ch[ii])){
                    gcl <- grd[to_use]*dcl[to_use]
                    ## gradient
                    gr[ii] <- max(sum(gcl) / sum(dcl[to_use]),min_grad)
                    ## topographic index
                    atb[ii] <- log(upa[ii]/gr[ii]) #log( upa[ii] / sum(gcl) )
                    ## if a hillslope propergate upslope area
                    if( !is.finite(ch[ii]) ){ ##then propogate area downslope
                        ## propogate area downslope
                        upa[ ngh ]  <- upa[ ngh ] + w*upa[ii]
                    }
                }else{
                    if( is.na(ch[ii]) ){
                    ## a hillslope cell that drains nowhere - this is an error
                        stop(paste("Cell",ii,"is a hillslope cell with no lower neighbours"))
                    }else{
                        stop("routing area from a channel cell...")
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

            ## merge upslope areas into the channel object
            ch_upa <- tapply(upa,ch,sum)
            ch_upa <- ch_upa[is.finite(ch_upa)]
            if( !all(private$shp$id == as.integer(names(ch_upa))) ){ stop("channels out of order...") }
            private$shp$up_area = as.numeric(ch_upa)

            ## compute catchment area to each reach
            ct_area <- private$shp$up_area
            for(ii in length(channel_routing):1){ ## since in reverse order
                if( length(channel_routing[[ii]] ) > 0 ){
                    ct_area[ channel_routing[[ii]] ] <- ct_area[ channel_routing[[ii]] ] +
                        ct_area[ii]/length( channel_routing[[ii]] )
                }
            }
            private$shp$ct_area <- ct_area

            ## save raster maps
            out <- terra::rast( private$brk[["dem"]], names="gradient", vals=gr )
            rstFile <- file.path(private$projectFolder,"gradient.tif")
            terra::writeRaster(out, rstFile);
            private$brk <- c( private$brk, terra::rast(rstFile))

            out <- terra::rast( private$brk[["dem"]], names="upslope_area", vals=upa )
            rstFile <- file.path(private$projectFolder,"upslope_area.tif")
            terra::writeRaster(out, rstFile); 
            private$brk <- c( private$brk, terra::rast(rstFile))

            out <- terra::rast( private$brk[["dem"]], names="atb", vals=atb )
            rstFile <- file.path(private$projectFolder,"atb.tif")
            terra::writeRaster(out, rstFile); 
            private$brk <- c( private$brk, terra::rast(rstFile))

            shpFile <- file.path(private$projectFolder,"channel.shp")
            terra::writeVector(private$shp, shpFile, overwrite=TRUE)

        },
        ## work out flow lengths to channel
        apply_flow_lengths = function(flow_length,verbose){

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
            fl <- d; fl[] <- NA
            bnd <- d; bnd[] <- NA
            
            

            ## distances and contour lengths
            ## distance between cell centres
            rs <- terra::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)

            ## compute channel bands
            if( verbose ){ print("Computing channel bands") }
            private$shp$band <- as.integer(NA)
            eN <- private$shp$endNode
            sN <- private$shp$startNode
            cbnd <- private$shp$band
            idx <- !(eN %in% sN) ##private$shp$id==1
            cnt <- 1
            while( any(idx) ){
                cbnd[idx] <- cnt
                idx <- eN %in% sN[idx]
                cnt <- cnt+1
            }
            private$shp$band <- cbnd

            ## compute channel routing
            if( verbose ){ print("Computing channel routing") }
            chn_route <- rep(list(integer(0),nrow(private$shp)))
            for(ii in nrow(private$shp):1){
                chn_route[[ii]] <- which( sN==eN[ii] )
            }
            

            
            ## if we go up in height order then we must have looked at all lower
            ## cells first
            idx <- order(d,na.last=NA)
            
            ## create flow direction storage
            fd <- matrix(as.numeric(NA),length(idx),9)
            colnames(fd) <- c("cell","topLeft","left","bottomLeft","top","bottom","topRight","right","bottomRight")
            fd_cnt <- 0
            
            n_to_eval <- length(idx)

            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }else{
                next_print <- Inf
            }

            w <- rep(0,8)
            for(ii in idx){
                if(is.finite(ch[ii])){
                    ## just channel
                    bnd[ii] <- cbnd[ ch[ii] ] ## relies in channels being ordered...
                    fl[ii]  <- 0
                    ## no need to increment count
                }else{
                    
                    ## it is not a channel
                    jdx <- ii+delta
                    grd <- (d[jdx]-d[ii])/dxy
                    gcl <- grd*dcl
                    is_lower <- is.finite(gcl) & gcl<0
                    ## compute weights to sum flow lengths
                    if(flow_length=="shortest"){
                        tmp <- (fl[jdx] + dxy) + (!is_lower * Inf)
                        w[] <- 0
                        w[which.min(tmp)] <- 1
                    }
                    if(flow_length=="d8"){
                        w[] <- 0
                        w[which.min(grd)] <- 1
                    }
                    if(flow_length=="dominant"){
                        w[] <- 0
                        w[which.min(gcl)] <- 1
                    }
                    if(flow_length=="expected"){
                        w[] <- 0
                        w[is_lower] <- gcl[is_lower] / sum( gcl[is_lower] )
                    }

                    ## round the weights
                    w <- round(w,2)
                    w <- w/sum(w)

                    bnd[ii] <- max(bnd[ jdx[w>0] ])+1
                    fl[ii] <- sum( (fl[jdx]+dxy)*w, na.rm=TRUE )
                    
                    fd_cnt <- fd_cnt + 1
                    fd[fd_cnt,] <- c(ii,w)
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
            ## out <- terra::rast( private$brk[["dem"]], names="band", vals=bnd )
            ## rstFile <- file.path(private$projectFolder,"band.tif")
            ## terra::writeRaster(out, rstFile); 
            ## private$brk <- c( private$brk, terra::rast(rstFile))
            
            out <- terra::rast( private$brk[["dem"]], names="band", vals=bnd)
            rstFile <- file.path(private$projectFolder,"band.tif")
            terra::writeRaster(out, rstFile); 
            private$brk <- c( private$brk, terra::rast(rstFile))

            out <- terra::rast( private$brk[["dem"]], names="flow_length", vals=fl )
            rstFile <- file.path(private$projectFolder,"flow_length.tif")
            terra::writeRaster(out, rstFile); 
            private$brk <- c( private$brk, terra::rast(rstFile))

            shpFile <- file.path(private$projectFolder,"channel.shp")
            terra::writeVector(private$shp, shpFile, overwrite=TRUE)

            saveRDS(chn_route,file.path(private$projectFolder,"channel_routing.rds"))
            saveRDS(fd[1:fd_cnt,],file.path(private$projectFolder,"flow_direction.rds"))
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

            
            private$brk <- c( private$brk,
                             terra::classify(x,rcl=brk,include.lowest=TRUE,
                                             filename= outFile, names=layer_name))

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
                brn <- terra::cover(x,brn)
            }
            ## add burns to pairs
            cp <- terra::cover(brn,cp)
            uq <- sort(terra::unique(cp)[[1]])
            cts <- c(-Inf,(uq[-1]+uq[-length(uq)])/2,Inf)
            cp <- terra::classify(cp,cts) + 1 ## classify returns numeric values starting at 0
            if(length(burns)>0){ brn <- terra::classify(brn,cts) +1 }
            
            
            ## TODO - replace with zone taking modal value
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

            outFile <- file.path(private$projectFolder,paste0(layer_name,".tif"))
            private$brk <- c( private$brk, terra::writeRaster( cp, outFile, names=layer_name))
            
            out <- list(type="combination",
                        groups=df)
            writeLines( jsonlite::toJSON(out), file.path(private$projectFolder,paste0(layer_name,".json")) )
            
        },
        
        ## create a model  
        apply_create_model = function(class_lyr,dist_lyr,delta_dist,
                                      rain_lyr,rainfall_label,
                                      pet_lyr,pet_label,layer_name,verbose,
                                      sf_opt,
                                      sz_opt){

            
##            print(system.time({
                
            rq <- c("gradient","filled_dem","channel",
                    class_lyr,dist_lyr,
                    rain_lyr,pet_lyr)
            has_rq <- rq %in% names(private$brk)            
            if(!all(has_rq)){
                stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
            }
            
            if(verbose){ cat("Setting up HRUs","\n") }
            
            ## make basic template based on sf_opt and sz_opt
            tmp_sf <- switch(sf_opt,
                             "cnst" = list(type = "cnst",
                                             parameters = c("c_sf" = 0.3, "d_sf" = 0.0,
                                                            "s_raf" = 0.0, "t_raf" = 999.9)),
                             "kin" = list(type = "kin",
                                          parameters = c("n" = 0.03,
                                                         "s_raf" = 0.0, "t_raf" = 999.9)),
                             stop("Unrecognised surface option")
                             )
            tmp_sz <- switch(sz_opt,
                             "exp" = list(type = "exp",
                                          parameters = c( "t_0" = 0.135, "m" = 0.04 )),
                             "bexp" = list(type = "bexp",
                                           parameters = c( "t_0" = 0.135, "m" = 0.04 , "h_sz_max" = 5)),
                             "dexp" = list(type = "dexp",
                                           parameters = c( "t_0" = 0.135, "m" = 0.04, "m2" = 0.1, "omega"=0.5)),
                             "cnst" = list(type = "cnst",
                                          parameters = c( "v_sz" = 0.1, "h_sz_max" = 5 )),
                             stop("Unrecognised saturated zone option")
                             )
            if(is.null(rain_lyr)){ tmp_prcp <- c("precip"=1) }else{ tmp_prcp <- numeric(0) }
            if(is.null(pet_lyr)){ tmp_pet <- c("pet"=1) }else{ tmp_pet <- numeric(0) }
            tmplate <- list(id = integer(0),
                            states = setNames(as.numeric(rep(NA,4)), c("s_sf","s_rz","s_uz","s_sz")),
                            properties = setNames(rep(0,4), c("area","width","Dx","gradient")),
                            sf = tmp_sf,
                            rz = list(type="orig", parameters = c("s_rzmax" = 0.1)),
                            uz = list(type="orig", parameters = c("t_d" = 8*60*60)),
                            sz = tmp_sz,
                            sf_flow_direction = numeric(0), #list(id = integer(0), fraction = numeric(0)),
                            sz_flow_direction = numeric(0), #list(id = integer(0), fraction = numeric(0)),
                            initialisation = c("s_rz_0" = 0.75, "r_uz_sz_0" = 1e-7),
                            precip = tmp_prcp,
                            pet = tmp_pet
                            )
            
            ## read in classification and distance layers
            cls <-  private$brk[[class_lyr]] ## channel values are NA
            dst <- private$brk[[dist_lyr]] ## !!!!We assume that all channel pixels have a distance that puts them in the correct order!!!!!
                     
            cls <- cls + max(private$shp$id) ## alter class so greater then river channel id
            hmap <-  terra::cover( private$brk[["channel"]], cls) ## make map of HRUs - but numbering not yet correct
            names(hmap) <- "hru"

            ## work out the order of the hrus and add class information
            tbl <- terra::zonal(dst,hmap,min) # minimum distance for each hillslope classification
            tbl$min_dst <- tbl[[dist_lyr]]
            tbl[[dist_lyr]] <- NULL
   
            ## add any missing channels
            if( !all(private$shp$id %in% tbl$hru ) ){ stop("where is the HRU!!!") }

            ## add class information
            jsonFile <- paste0(tools::file_path_sans_ext(terra::sources(private$brk[[class_lyr]])),".json")
            if( !file.exists(jsonFile) ){
                warning("No json file giving basis of the classifications")
            }else{
                tmp <- jsonlite::fromJSON(jsonFile)$groups
                tmp$hru <- tmp[[class_lyr]] + max(private$shp$id)
                tbl <- merge(tbl, tmp, by="hru", all.x=TRUE)
            }
            
            ##tbl$is_channel <- tbl$hru %in% private$shp$id
            tbl <- merge(tbl,as.data.frame(private$shp),by.x="hru",by.y="id",all.x=TRUE) ## add channel information

            ## order and renumber
            tbl <- tbl[ order(tbl$min_dst), ]
            hmap <- terra::subst(hmap, tbl$hru, 0:(nrow(tbl)-1)) ## this line is slow
            tbl$hru <- 0:(nrow(tbl)-1)

            ## make hrus
            hru <- lapply(1:nrow(tbl), function(ii){ ## slow but not as slow as subst..
                tmp <- as.list(tbl[ii,])
                out <- tmplate
                out$id <- tmp$hru
                tmp$hru <- NULL
                out$class <- tmp
                return(out)
            })
            
##            })) ## end of first system.time

##            print(system.time({
                
            ## work out the properties
            if(verbose){ cat("Computing properties","\n") }
     
            ## it is quickest to compute using blocks of a raster
            ## however for small rasters we will just treat as a single block
            d <- terra::as.matrix( private$brk[["filled_dem"]] , wide=TRUE )
            mp <- terra::as.matrix( hmap, wide=TRUE )
            dst <- terra::as.matrix( dst, wide=TRUE )
            gr <- terra::as.matrix( private$brk[["gradient"]], wide=TRUE )
            if( !is.null(rain_lyr) ){
                rain <- terra::as.matrix( private$brk[[rain_lyr]], wide=TRUE )
            }
            if( !is.null(pet_lyr) ){
                pet <- terra::as.matrix( private$brk[[pet_lyr]], wide=TRUE )
            }

##            })) ## end of second system.time
            ##atb <- terra::as.matrix( private$brk[["atb"]]  , wide=TRUE )

##            print(system.time({
                
            ## distances and contour lengths
            ## distance between cell centres
            rs <- terra::res( private$brk )
            dxy <- rep(sqrt(sum(rs^2)),8)
            dxy[c(2,7)] <- rs[1]; dxy[c(4,5)] <- rs[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(rs)
            nr <- nrow(mp); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            cellArea <- prod(rs)
                
            ## work out which cells to pass through
            idx <- which(is.finite(mp))
            n_to_eval <- c(length(idx),0,length(idx)/20)

            
            for(ii in idx){
                id <- mp[ii]
                jj <- id+1 ## index in the hru list (id +1)
                
                hru[[jj]]$properties["area"] <- hru[[jj]]$properties["area"] + 1
                ##hru[[jj]]$properties["atb_bar"] <- hru[[jj]]$properties["atb_bar"] + atb[ii]
                hru[[jj]]$properties["gradient"] <- hru[[jj]]$properties["gradient"] + gr[ii]
                
                ## work out subsurface flow for non-channel
                if( is.na(hru[[jj]]$class$startNode) & (dst[ii] <= (hru[[jj]]$class$min_dst + delta_dist)) ){
                    
                    ## look for flow direction      
                    ngh <- ii + delta ## neighbouring cells
                    nghid <- mp[ngh]
                    ## compute gradient
                    grd <- (d[ii]-d[ngh])/dxy
                    ## test if can be used
                    to_use <- is.finite(grd) & (grd > 0) & is.finite(nghid) & (nghid < id)
                    
                    if( any(to_use) ){
                        nghid <- paste(nghid[to_use])
                        tmp <- nghid[!(nghid %in% names(hru[[jj]]$sz_flow_direction))] 
                        if(length(tmp)>0){
                            hru[[jj]]$sz_flow_direction[ tmp ] <- 0
                        }
                        tmp <- tapply(grd[to_use]*dcl[to_use],nghid,sum)
                        hru[[jj]]$sz_flow_direction[ names(tmp) ] <-
                            hru[[jj]]$sz_flow_direction[ names(tmp) ] + tmp
                        hru[[jj]]$properties["width"] <- hru[[jj]]$properties["width"] + sum(dcl[to_use])
                    }
                    
                }

                ## work out precipitation
                if(!is.null(rain_lyr) && is.finite(rain[ii])){
                    nm <- paste0(rainfall_label,rain[ii])
                    if(!(nm %in% names(hru[[jj]]$precip))){
                        hru[[jj]]$precip[nm] <- 0
                    }
                    hru[[jj]]$precip[nm] <- hru[[jj]]$precip[nm] + 1
                }
                
                ## work out pet
                if(!is.null(pet_lyr) && is.finite(pet[ii])){
                    nm <- paste0(pet_label,pet[ii])
                    if(!(nm %in% names(hru[[jj]]$pet))){ 
                        hru[[jj]]$pet[nm] <- 0
                    }
                    hru[[jj]]$pet[nm] <- hru[[jj]]$pet[nm] + 1
                }   
                    
                n_to_eval[2] <- n_to_eval[2] + 1
                if( verbose & (n_to_eval[2] > n_to_eval[3]) ){
                    cat(round(100*n_to_eval[3] / n_to_eval[1],1),
                        "% complete","\n")
                    n_to_eval[3] <- n_to_eval[3] + n_to_eval[1]/20
                }
                
            }

##            })) ## end of third system.time

##            print(system.time({
                
            ## second pass to correct sumations and compute surface
            sN <- sapply(hru,function(h){h$class$startNode})

            outlet_id <- NULL
            no_outflow <- NULL
            for(ii in 1:length(hru)){

                ## check precip
                if(!is.null(rain_lyr)){
                    if( sum(hru[[ii]]$precip) != hru[[ii]]$properties["area"] ){
                        warning(paste("HRU",hru[[ii]]$id," - Precip cell count is",sum(hru[[ii]]$precip),
                                      "full cell count is",hru[[ii]]$properties["area"]))
                    }
                }
                hru[[ii]]$precip <- list(name = as.character(names(hru[[ii]]$precip)),
                                         fraction = as.numeric(hru[[ii]]$precip/sum(hru[[ii]]$precip)) )
                
                ## check pet
                if(!is.null(pet_lyr)){
                    if( sum(hru[[ii]]$pet) != hru[[ii]]$properties["area"] ){
                        warning(paste("HRU",hru[[ii]]$id," - PET cell count is",sum(hru[[ii]]$precip),
                                      "full cell count is",hru[[ii]]$properties["area"]))
                    }
                }
                
                hru[[ii]]$pet <- list(name = as.character(names(hru[[ii]]$pet)),
                                      fraction = as.numeric(hru[[ii]]$pet/sum(hru[[ii]]$pet)) )
                
                ##hru[[ii]]$properties["atb_bar"] <- hru[[ii]]$properties["atb_bar"] / hru[[ii]]$properties["area"]
                hru[[ii]]$properties["gradient"] <- hru[[ii]]$properties["gradient"] / hru[[ii]]$properties["area"]
                hru[[ii]]$properties["area"] <- hru[[ii]]$properties["area"] * cellArea
                
                if( !is.na(hru[[ii]]$class$startNode) ){
                    ## then a channel
                    hru[[ii]]$properties["Dx"] <- hru[[ii]]$class$length
                    hru[[ii]]$properties["width"] <- hru[[ii]]$class$width
                    hru[[ii]]$properties["gradient"] <- hru[[ii]]$class$slope ## PJS TODO document this and requirement to have it as a variable
                    jdx <- which( sN == hru[[ii]]$class$endNode)
                    if(length(jdx) >0){
                        ## has downstream connenction
                        hru[[ii]]$sf_flow_direction <- list(
                            id = as.integer( jdx - 1 ),
                            fraction = rep( as.numeric( 1 / length(jdx) ),length(jdx) ))
                    }else{
                        ## goes to an outlet
                        hru[[ii]]$sf_flow_direction <- list(id = integer(0),fraction=numeric(0))
                        outlet_id <- c(outlet_id,hru[[ii]]$id)
                    }
                    ## set subsruface to match surface
                    hru[[ii]]$sz_flow_direction <- hru[[ii]]$sf_flow_direction
                    
                }else{
                    if( length(hru[[ii]]$sz_flow_direction)==0 ){
                        no_outflow <- c(no_outflow, hru[[ii]]$id)
                    }
                    
                    hru[[ii]]$sz_flow_direction <- list(
                        id = as.integer(names( hru[[ii]]$sz_flow_direction )),
                        fraction = as.numeric( hru[[ii]]$sz_flow_direction/sum(hru[[ii]]$sz_flow_direction) ))

                    ## set surface to match subsurface
                    hru[[ii]]$sf_flow_direction <- hru[[ii]]$sz_flow_direction
                    hru[[ii]]$properties["Dx"] <- hru[[ii]]$properties["area"]/hru[[ii]]$properties["width"]
                }
            }

##            })) ## end of system.time
            
            ## stop if any points with no outflow that aren't channels
            if( length(no_outflow) >0 ){
                stop( "The following HRUs have no valid outflows and are not channels: \n",
                     paste(no_outflow, collapse=", "))
            }
            
            ## ############################################
            ## Add outlet at all outlets from river network
            ## ############################################
            ##idx <- sapply(hru,function(x){ifelse(length(x$sf_flow_direction$id)==0,x$id,NULL)})
            output_flux<- data.frame(name = paste0("q_sf_",outlet_id),
                                           id = as.integer( outlet_id ), flux = "q_sf")
                        
            ## ############################
            ## save model
            model <- list(hru=hru, output_flux = output_flux)
            model$map <- paste0(layer_name,".tif")
            terra::writeRaster(hmap,model$map,overwrite=TRUE)
            saveRDS(model,paste0(layer_name,".rds"))
        }        
    )
    )
