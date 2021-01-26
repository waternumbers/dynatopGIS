#' R6 Class for processing a catchment to make a Dynamic TOPMODEL
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
        #' @details This loads the meta data file found at \code{meta_path}, or creates it with a warning if no file is present. It \code{check} is \code{TRUE} then the meta data file contents are checked witht he level of returned information being controlled by \code{verbose}.
        #'
        #' @return A new `dynatopGIS` object
        initialize = function(meta_file, check=TRUE,verbose=TRUE){
            ## create directory if it doesn't exist
            
            private$meta_path <- meta_file
            private$wdir <- normalizePath(dirname(meta_file))

            ## if new it is an empty folder so write in meta data
            is_new <- !file.exists(meta_file)
            if(is_new){
                warning("Creating meta file at",meta_file)
                private$write_meta()
            }                
            
            ## load meta data
            if(check){ private$check_meta(verbose) }
            
            invisible(self)
        },
        #' @description Get project emta data
        get_meta = function(){ private$meta },
        #' @description Get current working directory
        #' @details Newly generated layers are added to the working directory. By deafult this is the directory containing the meta date file.
        get_working_directory = function(){private$wdir},
        #' @description Set current working directory
        #'
        #' @param file_path the path to the new directory to create
        #' @param create should the driectory be created if it doesn't exist
        #'
        #' @details Newly generated layers are added to the working directory. By deafult this is the directory containing the meta date file.
        set_working_directory = function(file_path,create=TRUE){
            if(!dir.exists(file_path) & create){ dir.create(file_path,recursive=TRUE) }
            if(!dir.exists(file_path)){ stop("Directory does not exist") }
            private$wdir <- file_path
            invisible(self)
        },
        #' @description Import channel data from an OGR file to the `dynatopGIS` object
        #'
        #' @param sp_object a spatialLinesDataFrame or SpatialPoltgonsDataFrame containing the channel information
        #' @param property_names named vector of columns of the spatial data frame to use for channel properties
        #' @param default_width default width of a channel if not specified in property_names. Defaults to 2 metres.
        #'
        #' @details Takes the input file and covnerts it to polygons with properties length,width,startNode,endNode. The variable names in the sp_object data frame which corresponding to these properties can be specified in the \code{property_names} vector.
        #' The function opulates the channel_id & channel_area layers and alters the land_area layer to correspond.
        #'
        #' @return suitable for chaining              
        add_dem = function(dem,fill_na=TRUE,verbose=FALSE){
            private$apply_add_dem(dem,fill_na,verbose)
            invisible(self)
        },
        #' @description Import channel data from an OGR file to the `dynatopGIS` object
        #'
        #' @param sp_object a spatialLinesDataFrame or SpatialPoltgonsDataFrame containing the channel information
        #' @param property_names named vector of columns of the spatial data frame to use for channel properties
        #' @param default_width default width of a channel if not specified in property_names. Defaults to 2 metres.
        #'
        #' @details Takes the input file and covnerts it to polygons with properties length,width,startNode,endNode. The variable names in the sp_object data frame which corresponding to these properties can be specified in the \code{property_names} vector.
        #' The function opulates the channel_id & channel_area layers and alters the land_area layer to correspond.
        #'
        #' @return suitable for chaining
        add_channel = function(channel,property_names=c(length="length",
                                                        startNode="startNode",
                                                        endNode="endNode",
                                                        width="width"),
                               default_width=2){
            
            if("channel" %in% names(private$find_layer(TRUE))){
                stop("The channel exists, start a new project")
            }
            private$apply_add_channel(channel,property_names,default_width)
            invisible(self)
        },
        #' @description Add a layer of geographical information
        #' @param raster_layer a RasterLayer object to add
        #' @param layer_name name to give to the layer
        #' @param check_name logical indicating if name should be checked against reserved names
        #' @details The default is to check the name so that names reserved for computed varaibles are not overwritten. Disabling to overwrite computed layer may have unintended consequences.
        #' @return suitable for chaining
        add_layer = function(file_path,layer_name){
            layer_name <- as.character(layer_name)
            file_path <- normalizePath(file_path)
            
            pos_val <- private$find_layer(TRUE)
            ## check layer names is acceptable
            if( layer_name %in% private$find_layer() ){ stop("Layer name is already in use") }
            if( file_path %in% pos_val ){
                stop("File already in use as layer", names(pos_val[pos_val==file_path]))
            }

            private$check_rst(file_path,layer_name)
            private$meta$layers[[layer_name]] <- list(file=file_path,type="user")
            private$write_meta()
            invisible(self)
        },
        #' @description Get a layer of geographical information or a list of layer names
        #' @param layer_name name of the layer give to the layer
        #' @return a RasterLayer of the requested information if layer_name is given else a vector of layer names
        get_layer = function(layer_name=character(0)){
            pos_val <- private$find_layer(TRUE)
            
            ## handle case where a list of layers is requested
            if( length(layer_name) == 0 ){
                return( names(pos_val) )
            }
            ## check layer name exists
            layer_name <- match.arg(layer_name,names(pos_val))
            
            ## make raster and return
            fn <- pos_val[layer_name]
            if(layer_name=="channel"){
                return( raster::shapefile(fn) )
            }else{
                return( raster::raster(fn) )
            }
            
        },
        #' @description Plot a layer
        #' @param layer_name the name of layer to plot
        #' @param add_channel shouel the channel be added to the plot
        #' @return a plot
        plot_layer = function(layer_name,add_channel=TRUE){
            lyr <- self$get_layer(layer_name)
            raster::plot( lyr, main = layer_name)
            if( add_channel ){
                chn <- self$get_layer("channel")
                raster::plot( chn, add=TRUE )
            }
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
            private$apply_sink_fill(min_grad,max_it,verbose,hot_start)
            invisible(self)
        },
        #' @description Computes area maps and presence of channel in dem pixels
        #'
        #' @details The algorithm calculates the land and channel area for each DEM pixel assigning a channel_id to each pixel with a channel area.
        compute_areas=function(){
            private$apply_compute_areas()
            invisible(self)
        },
        #' @description Computes statistics e.g. gradient, log(a/tanb) for raster cells
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
        #' @details The algorithm passed through the cells in increaing height. For measures of flow length to the channel are computed. The shortest length (minimum length to channel through any flow path), the dominant length (the length taking the flow diection with the highest fraction for each pixel on the path) and expected flow length (flow length based on sum of downslope flow lengths based on fraction of flow to each cell) and band (strict sequence to ensure that all contributing cell have a higher band value). By definition cells in the channel that have no land area have a length (or band) of NA.
        compute_flow_lengths = function(verbose=FALSE){
            private$apply_flow_lengths(verbose)
            invisible(self)
        }, 
        #' @description Create a catchment classification based on any number of landscape layer cuts or burns
        #' @param cuts A named list of cuts of make to form the HRU. Names should correspond to raster layers in the project directory. Values should be numeric and define either the number of bands (single value) or breaks between band (multiple values)
        #' @param burns a vector of layer names which are to be burnt on
        #'
        #' @details This applies the given cuts to the supplied landscape layers to produce areal groupings of the catchment. Burns are added directly in the order they are given. Not that cuts are impliment using \code{raster::cut} with \code{include.lowest = TRUE}. Note that is specifying a vector of cuts values outside the limits will be set to NA.
        classify = function(layer_name,cuts,burns=NULL){
            private$apply_classify(cuts,burns,layer_name)
            invisible(self)
        },

        #' @description Compute a Dynamic TOPMODEL
        #'
        #' @param class_layer the layer defining the topographic classes
        #' @param dist_layer the layer dfining the disatance to the channel
        #' @param brk breaks for cutting up the distance layer - see classify
        #' @param layer_name name for the new model and layers
        #' @param delta_min see details
        #' @param verbose print more details of progress
        #'
        #' @details The distance layer is cut into classes and combined with \code{class_layer} using classify to produce the HSUs. Flow between HSUs is computed based on the raster cells that drain to HSU of a lower distance class (closer to the bottom of the hillslope). If no such cells exist then the flow is determined by the cells whose distance is less then the lowest diatance within the HSU plus delta_min. 
        #' @param verbose print out additional diagnostic information        
        create_model = function(class_layer,dist_layer,brk,layer_name,
                                delta_min=sqrt(sum(private$meta$resolution^2)),
                                verbose=FALSE){
            private$apply_create_model(class_layer,dist_layer,brk,layer_name,delta_min,verbose)            
            invisible(self)
        },
        #' @description Compute HSU fractional inputs
        #'
        #' @param hsu_layer layer of the model
        #' @param input_raster raster of values containing the input class
        #' @param hsu_label string containing the label to add the values of hsu_layer
        #' @param input_label string containing the label to add the values of input_layer
        #'
        #' @details Fractions are computed for all HSUs (including the channel) using raster::crosstab. It does not accout for unequal cells areas or HSUs which only cover a fraction of a cell
        input_frac = function(hsu_layer,input_layer,hsu_label="hsu_",input_label=""){
            private$apply_input_frac(hsu_layer,input_layer,hsu_label,input_label)
        },            
        #' @description get the version number
        #' @return a numeric version number
        #' @details teh version number indicates the version of the algorithms within the object
        get_version = function(){
            private$version
        },
        #' @description get the cuts and burns used to classify
        #' @return a list with two elements, cuts and burns
        get_class_method = function(layer_name){
            layer_name <- match.arg(layer_name,private$find_layer())
            if( private$meta$layers[[layer_name]]$type == "classification" ){
                private$meta$layers[[layer_name]]$method
            }else{
                stop(layer_name," is not a classification")
            }
        },
        
        #' @description DEPRECIATED Get the model is a form suitable for dynatop
        #' @return a Dynamic TOPMODE description suitable for the \code{dyantop} package
        get_model = function(){
            stop("This function is depreciated - model is written out on generation")
        },
        #' @description DEPRECIATED Plot the spatial properties of a model
        #' @param add_channel should the channel be added to the plot
        plot_model = function(add_channel=TRUE){
            stop("This function is depreciated\n",
                 " - use $plot_layer(<layer_name>) where layer_name is name of the model"
                 )
        },
        #' @description DEPRECIATED Return the index of neighbouring cells
        #' @param x index of cells for which to find neighbours
        #' @return a list of indexes of neighbours
        neighbour = function(x){
            stop("Function is depreciated")
        },
        #' @description DEPECIATED Get a classification layer of geographical information or a list of layer names
        #' @param layer_name name of the layer give to the layer
        #' @return a RasterLayer of the requested information if layer_name is given else a vector of layer names
        get_class = function(layer_name){
            stop("This function is no longer valid: \n",
                 " - use $get_layer(<layer_name>) to see output \n",
                 " - use $geta_class_method(<layer_name>) to see the splits and cuts")
        },
        
        #' @description DEPECIATED Plot a classification layer
        #' @param layer_name the name of layer to plot
        #' @param add_channel shouel the channel be added to the plot
        #' @return a plot
        plot_class = function(layer_name="final",add_channel=TRUE){
            warning("This function is no longer valid",
                    "- use $plot_layer(<layer_name>) to see output")
            self$plot_layer(layer_name,add_channel)
        },
        #' @description DEPRECIATED Get the channel description
        #' @return an SpatialPolygonDataFrame containing the channels 
        get_channel = function(){
            warning("$get_channel() Depreciated use $get_layer('channel')")
            self$get_layer("channel")
        }
               
    ),
    private = list(
        version = "0.2.0.9010",
        wdir=character(0),
        meta_path = character(0),
        meta=list(
            crs=character(0),
            extent=numeric(0),
            resolution=numeric(0),
            layers= list(
                dem = list(type="reserved",file=character(0)),
                channel = list(type="reserved",file=character(0)),
                filled_dem = list(type="reserved",file=character(0)),
                land_area = list(type="reserved",file=character(0)),
                channel_area = list(type="reserved",file=character(0)),
                channel_id = list(type="reserved",file=character(0)),
                gradient = list(type="reserved",file=character(0)),
                upslope_area = list(type="reserved",file=character(0)),
                contour_length = list(type="reserved",file=character(0)),
                atb = list(type="reserved",file=character(0)),
                band = list(type="reserved",file=character(0)),
                shortest_flow_length = list(type="reserved",file=character(0)),
                dominant_flow_length = list(type="reserved",file=character(0)),
                expected_flow_length = list(type="reserved",file=character(0))
            )
        ),
        make_filename = function(fn,is_shp=FALSE){
            if(!is_shp){
                return( file.path(private$wdir,paste0(fn,'.tif')) )
            }else{
                return( file.path(private$wdir,paste0(fn,'.shp')) )
            }
        },
        read_meta = function(){
            #browser()
            if(!file.exists(private$meta_path)){ stop("Missing metadata file") }
            meta <- jsonlite::fromJSON(private$meta_path)            
            if(length(meta$crs)>0){meta$crs <- crs(meta$crs)}
            if(length(meta$extent)>0){meta$extent <- extent(meta$extent)}
            private$meta <- meta
        },
        write_meta = function(){
            meta <- private$meta
            meta$crs <- crs(meta$crs,asText=TRUE)
            if("Extent" %in% class(meta$extent)){
                meta$extent <- c(meta$extent@xmin,meta$extent@xmax,meta$extent@ymin,meta$extent@ymax)
            }
            writeLines(jsonlite::toJSON(meta,null="null",pretty=TRUE),
                       file.path(private$meta_path))
        },
        find_layer = function(file_name=FALSE){
            if(file_name){
                unlist(sapply(private$meta$layers,function(x){as.character(x$file)}))
            }else{
                names(private$meta$layers)
            }
        },
        
                
        ## check raster layer
        check_rst = function(fn,nm){
            rst <- fn
            if(!("RasterLayer" %in% class(fn))){ rst <- raster::raster(fn) }
            if(!("RasterLayer" %in% class(rst))){ stop(nm," is not a RasterLayer") }
            
            if(!compareCRS(rst,private$meta$crs)){ stop(nm," projection does not match meta data") }
            if(!(extent(rst)==private$meta$extent)){ stop(nm," extent does not match meta data") }
            if(!all(res(rst)==private$meta$resolution)){ stop(nm," resolution does not match meta data") }
            NULL
        },
        ## function to check meta data
        check_meta = function(verbose){
            ## browser()
            private$read_meta()
            ## TODO
            warning("No checks on the meta are currently performed")
        },

        ## adding dem
        apply_add_dem = function(dem,fill_na,verbose){
            
            if("dem" %in% names(private$find_layer(TRUE))){
                stop("The DEM exists, start a new project")
            }

            if(!("RasterLayer" %in% class(dem))){
                ## then dem should be a character string points to the file
                dem <- as.character(dem)
                if(file.exists(as.character(dem))){
                    dem <- raster::raster(dem)
                }else{
                    stop("dem should be either a RasterLayer or the path to a file which can be read as a RasterLayer")   
                }
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
                dem[!is.na(na_clumps)] <- -Inf # set to low value to indicate missing
            }

            
            if(length(private$meta$crs)==0){ private$meta$crs = crs(dem) }
            if(length(private$meta$extent)==0){ private$meta$extent <- extent(dem) }
            if(length(private$meta$resolution)==0){ private$meta$resolution <- res(dem) }
            private$check_rst(dem,"DEM")

            fn <- private$make_filename("dem")
            writeRaster(dem,fn); private$meta$layers[["dem"]]$file <- fn
            private$write_meta()
        },

        ## add the channel
        apply_add_channel = function(sp_object,property_names,default_width){
            ##browser()
            ## check property_names
            if( !all( c("length","startNode","endNode") %in% names(property_names)) ){
                stop("A required property name is not specified")
            }

            if(!is(sp_object,"SpatialLinesDataFrame") && !is(sp_object,"SpatialPolygonsDataFrame")){
                ## then sp_object should be a character string points to the file
                sp_object <- as.character(sp_object)
                if(file.exists(as.character(sp_object))){
                    sp_object <- raster::shapefile(sp_object)
                }else{
                    stop("sp_object is not a spatial lines or polygon data frame object (even when read in)")
                }
            }

            ## populate meta data from channel if missing
            if(length(private$meta$crs)==0){ private$meta$crs = crs(sp_object) }
            
            ## check project of the sp_object
            if( !compareCRS(crs(sp_object),private$meta$crs) ){
                stop("Projection of channel object does not match")
            }
                        
            ## check if there is an id feild which will be overwritten
            if( ("id" %in% names(sp_object)) ){
                warning("The name id is reserved and will be overwritten",
                        "Original id values moved to original_id")
                sp_object[["original_id"]] <- sp_object[["id"]]
            }

            ## ensure correct columns are in the sp_object
            sp_object[["id"]] <- 1:length(sp_object)
            for(ii in names(property_names)){
                sp_object[[ii]] <- sp_object[[property_names[ii]]]
            }

            ## convert factors to strings
            idx <- sapply(sp_object@data, is.factor)
            sp_object@data[idx] <- lapply(sp_object@data[idx], as.character)

            ## add width if required
            if(!is(sp_object,"SpatialPolygonsDataFrame")){
                if(!("width" %in% names(sp_object))){
                    warning("Modifying to spatial polygons using default width")
                    sp_object[['width']] <- default_width
                }else{
                    warning("Modifying to spatial polygons using provided width")
                }
                
                ## buffer to convert to spatial polygons object
                sp_object <- rgeos::gBuffer(sp_object, byid=TRUE, width=sp_object[['width']])
            }

            file_name <- private$make_filename("channel",TRUE)
            raster::shapefile(sp_object,file_name)
            private$meta$layers[["channel"]]$file <- file_name
            private$write_meta()
        },

        ## calculate area maps
        apply_compute_areas = function(){
            
            rq <- c("dem","channel")
            pos_val <- private$find_layer(TRUE)
            has_rq <- rq %in% names(pos_val)
            if(!all(has_rq)){
                stop("Missing files:\n",paste(rq[!has_rq],sep="\n"))
            }
            
            dem <- raster::raster(pos_val["dem"])
            chn <- raster::shapefile(pos_val["channel"])

            ## compute land area ignoring channels
            land_area <- raster::area(dem, na.rm=TRUE)
            ## compute the areas taken up by the channels
            ch_area <- raster::area(chn) # areas of river channels
            
            ## extract cells index, and fraction of river area in cell
            ch_cell<- raster::extract(land_area,chn,weights=TRUE,
                                      cellnumbers=TRUE,na.rm=TRUE)

            
            ## initialise the final rasters
            channel_area <- dem; channel_area[!is.na(dem)] <- 0
            channel_id <- dem*NA

            ## Loop and add the chanel id and channel area
            for(ii in 1:length(ch_cell)){
                idx <- ch_cell[[ii]][,'cell']
                la <- ch_cell[[ii]][,'value']
                ca <- ch_area[ii]*ch_cell[[ii]][,'weight']
                ca <- pmin(ca,la)
                jdx <- !is.na(channel_area[idx]) & (ca > channel_area[idx])
                channel_area[ idx[jdx] ] <- ca[jdx]
                channel_id[ idx[jdx] ] <- chn[["id"]][ii]
            }

            ## correct land area
            land_area <- land_area - channel_area
            
            ## write out rasters
            fn <- private$make_filename("land_area")
            writeRaster(land_area,fn); private$meta$layers[['land_area']]$file <- fn
            fn <- private$make_filename("channel_area")
            writeRaster(channel_area,fn); private$meta$layers[['channel_area']]$file <- fn
            fn <- private$make_filename("channel_id")
            writeRaster(channel_id,fn); private$meta$layers[['channel_id']]$file <- fn
            private$write_meta()
        },

        ## Sink fill
        ## TODO: this is not handling is_valid as for original looping algorithm
        ## this would make it much more efficent
        apply_sink_fill = function(min_grad,max_it,verbose,hot_start){
            
            rq <- ifelse(hot_start,
                         c("filled_dem","channel_id"),
                         c("dem","channel_id"))
            pos_val <- private$find_layer(TRUE)
            has_rq <- rq %in% names(pos_val)
            if(!all(has_rq)){
                stop("Not all required input files have been generated \n",
                     "Try running compute_areas first")
            }

            d <- raster::raster(pos_val[rq[1]])
            ch <- raster::raster(pos_val["channel_id"])

            rfd <- d ## initialise the output

            ## convert to matrices - much quicker but will fail for large dem's
            d <- as.matrix(d)
            ch <- as.matrix(ch)

            ## values that should be valid
            to_be_valid <- !is.na(d) & is.na(ch)  # all values not NA should have a valid height
            is_valid <- is.finite(ch) & !is.na(d) # TRUE if a channel cell for initialisation
            changed <- is_valid # cells changed at last iteration
            fd <- to_be_valid*Inf; fd[is_valid] <- d[is_valid]
            to_eval <- d; to_eval[] <- FALSE

            ## distance between cell centres
            dxy <- matrix(sqrt(sum(private$meta$resolution^2)),3,3)
            dxy[1,2] <- dxy[3,2] <- private$meta$resolution[2]
            dxy[2,1] <- dxy[2,3] <- private$meta$resolution[1]
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
                cat("Iteration",it,"\n")
                cat("\t","Cells to evaluate:",sum(to_eval),"\n")
                cat("\t","Percentage Complete:",
                    round(100*sum(is_valid)/sum(!is.na(d)),1),"\n") #to_be_valid),1),"\n")
                
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

            rfd[] <- fd
            ## write out filled dem
            fn <- private$make_filename("filled_dem")
            writeRaster(rfd,fn); private$meta$layers[['filled_dem']]$file <- fn
            private$write_meta()
            
            if(it>max_it){ stop("Maximum number of iterations reached, sink filling not complete") }
            
        },

        ## Function to compute the properties
        apply_compute_properties = function(min_grad,verbose){

            rq <- c("filled_dem","channel_id","land_area")
            pos_val <- private$find_layer(TRUE)
            has_rq <- rq %in% names(pos_val)
            if(!all(has_rq)){
                stop("Not all required input files have been generated \n",
                     "Try running compute_areas first and sink_fill first")
            }
            
            ## load rasters
            d <- raster::raster(pos_val["filled_dem"])
            ch <- raster::raster(pos_val["channel_id"])
            la <- raster::raster(pos_val["land_area"])

            n_to_eval <- cellStats(is.finite(d),sum)
            out <- d ## template for output
            
            ## it is quickest to compute using blocks of dem as a raster
            ## however for small rasters we will just treat as a single block
            ## assumes the raster is padded with NA
            d <- as.matrix(d)
            ch <- as.matrix(ch)
            la <- as.matrix(la)

            idx <- order(d,decreasing=TRUE,na.last=NA)
            idx <- idx[ la[idx]>0 ]
            
            n_to_eval <- length(idx)
            
            ## distances and contour lengths
            ## distance between cell centres
            dxy <- rep(sqrt(sum(private$meta$resolution^2)),8)
            dxy[c(2,7)] <- private$meta$resolution[1]; dxy[c(4,5)] <- private$meta$resolution[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(private$meta$resolution)
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
            ## initialise output
            gr <- cl <- atb <- d*NA
            upa <- la

            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }
            
            for(ii in idx){
                is_channel <- is.finite(ch[ii]) ## if a channel
                ngh <- ii + delta ## neighbouring cells
                
                ## compute gradient
                grd <- (d[ii]-d[ngh])/dxy
                
                ## process depending upon whether it is a channel cell...
                if( is_channel ){
                    ## TODO these cells should have ath and cl as well!!!
                    ## only need to evaluate gradient and atanb
                    ## cells that drain *into* the one being evaluated
                    to_use <- is.finite(grd) & (grd < 0) &
                        !is.finite(ch[ngh])
                    if(any(to_use)){
                        ## if upstream cells draining in
                        gr[ii] <- max(min_grad,-sum(grd[to_use]*dcl[to_use]) / sum(dcl[to_use]))
                    }else{
                        ## if no upstream cells set default gradient
                        gr[ii] <- min_grad
                    }
                    cl[ii] <- dxy[1] # diagonal distance
                    atb[ii] <- log( upa[ii] / gr[ii] )
                }else{
                    ## is not a channel - need gradient and to pass area along
                    to_use <- is.finite(grd) & (grd > 0)
                    if(any(to_use)){
                        gcl <- grd[to_use]*dcl[to_use]
                        ## gradient
                        gr[ii] <- max(sum(gcl) / sum(dcl[to_use]),min_grad)
                        ## contour length
                        cl[ii] <- sum(dcl[to_use])
                        ## topographic index
                        atb[ii] <- log(upa[ii]/gr[ii]) #log( upa[ii] / sum(gcl) )
                        ## fraction of flow in each direction
                        frc <- gcl/sum(gcl)
                        ## propogate area downslope
                        
                        upa[ ngh[to_use] ]  <- upa[ ngh[to_use] ] + frc*upa[ii]
                    }else{
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
            
            

            ## write out raster maps
            fn <- fn <- private$make_filename("gradient")
            out[] <- gr; writeRaster(out,fn); private$meta$layers[['gradient']]$file <- fn
            
            fn <- private$make_filename("upslope_area")
            out[] <- upa; writeRaster(out,fn); private$meta$layers[['upslope_area']]$file <- fn
            
            fn <- private$make_filename("contour_length")
            out[] <- cl; writeRaster(out,fn); private$meta$layers[['contour_length']]$file <- fn

            
            fn <- private$make_filename("atb")
            out[] <- atb; writeRaster(out,fn); private$meta$layers[['atb']]$file <- fn

            private$write_meta()
        },

        ## work out flow lengths
        apply_flow_lengths = function(verbose){

            rq <- c("filled_dem","channel_id","land_area")
            pos_val <- private$find_layer(TRUE)
            has_rq <- rq %in% names(pos_val)
            if(!all(has_rq)){
                stop("Not all required input files have been generated \n",
                     "Try running compute_areas first and sink_fill first")
            }

            ## load rasters
            d <- raster::raster(pos_val["filled_dem"])
            ch <- raster::raster(pos_val["channel_id"])
            la <- raster::raster(pos_val["land_area"])

            out <- d

            ## convert to matrix for speed
            d <- as.matrix(d)
            ch <- as.matrix(ch)
            la <- as.matrix(la)
            
            ## create some distance matrices
            sfl <- d; sfl[] <- NA
            dfl <- d; dfl[] <- NA
            efl <- d; efl[] <- NA
            bnd <- d; bnd[] <- NA

            ## distances and contour lengths
            ## distance between cell centres
            dxy <- rep(sqrt(sum(private$meta$resolution^2)),8)
            dxy[c(2,7)] <- private$meta$resolution[1]; dxy[c(4,5)] <- private$meta$resolution[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(private$meta$resolution)
            nr <- nrow(d); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
            ## if we go up in height order then we must have looked at all lower
            ## cells first
            idx <- order(d,na.last=NA)
            
            n_to_eval <- length(idx)
            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }
            for(ii in idx){
                if(is.finite(ch[ii])){
                    ## it is a channel
                    if(la[ii]>0){
                        ## mixed cell
                        bnd[ii] <- 1
                        sfl[ii]  <- dfl[ii] <- efl[ii] <- sqrt(sum(private$meta$resolution^2))/2
                    }else{
                        ## just channel
                        bnd[ii] <- 0
                        sfl[ii]  <- dfl[ii] <- efl[ii] <- 0
                    }
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

            ## set channel cells to NA
            idx <- is.finite(ch) & !(la>0)
            bnd[idx] <- sfl[idx] <- dfl[idx] <- efl[idx] <- NA
            
                        
            ## write out raster files
            fn <- private$make_filename("band")
            out[] <- bnd; writeRaster(out,fn); private$meta$layers[['band']]$file <- fn
            fn <- private$make_filename("shortest_flow_length")
            out[] <- sfl; writeRaster(out,fn);
            private$meta$layers[["shortest_flow_length"]]$file <- fn

            fn <- private$make_filename("dominant_flow_length")
            out[] <- dfl; writeRaster(out,fn);
            private$meta$layers[['dominant_flow_length']]$file <- fn
            
            fn <- private$make_filename("expected_flow_length")
            out[] <- efl; writeRaster(out,fn);
            private$meta$layers[['expected_flow_length']]$file <- fn

            private$write_meta()
        },
        
        ## split_to_class
        apply_classify = function(cuts,burns,layer_name){

            ## check all cuts and burns are in possible layers
            rq <- c(names(cuts),burns)
            layer_name <- as.character(layer_name)
            pos_val <- private$find_layer(TRUE)

            has_rq <- rq %in% names(pos_val)
            if(!all(has_rq)){
                stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
            }

            ## check layer_name isn;t already used
            if(layer_name %in% names(pos_val)){
                stop("layer_name is already used")
            }
                
            ## work out new cuts by cantor_pairing
            init <- TRUE
            for(ii in names(cuts)){
                brk <- as.numeric(cuts[[ii]])
                fn <- pos_val[ii] ## work out filename
                x <- raster::raster(fn)
                if( !(length(brk)==1 & is.na(brk[1])) ){
                    
                    if(length(brk)==1){
                        ## this defines brks in the same way as cut would otherwise
                        rng <- cellStats(x,range)
                        brk <- seq(rng[1],rng[2],length=brk+1)
                        ## the following does quantile based cuts
                        #brk <- brk+1
                        #brk <- quantile(x,probs=(0:brk)/brk)
                        cuts[[ii]] <- brk
                    }
                    ## else{
                    ##     brk <- unqiue(brk,rng)
                    ##     brk <- brk[brk<=rng[2] & brk>=rng[1]]
                    ## }
                    x <- cut(x,breaks=brk,include.lowest=TRUE)
                }
                if(init){
                    cp <- x
                    init <- FALSE
                }else{
                    cp <- 0.5*(cp+x)*(cp+x+1)+x
                    uq <- sort(unique(cp))
                    cp <- cut(cp, c(-Inf,(uq[-1]+uq[-length(uq)])/2,Inf))
                }
            }

            ## add in burns
            for(ii in burns){
                fn <- fn <- pos_val[ii,"file"] ## work out filename
                x <- raster:raster(fn)
                idx <- Which(is.finite(x))
                cp[idx] <- x[idx]
            }
            
            fn <- private$make_filename(layer_name)
            writeRaster(cp,fn)
            private$meta$layers[[layer_name]] <- list(file=fn,type="classification",
                                                    method=list(cuts=cuts,burns=burns))
            private$write_meta()
        },
        
        ## create a model  
        apply_create_model = function(class_lyr,dist_lyr,brk,layer_name,delta_min,verbose){

            rq <- c("atb","gradient","land_area","filled_dem",
                    "channel","channel_area","channel_id",
                    class_lyr,dist_lyr)
            pos_val <- private$find_layer(TRUE)

            has_rq <- rq %in% names(pos_val)
            if(!all(has_rq)){
                stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
            }
            if(!private$meta$layers[[class_lyr]]$type == "classification"){
                stop(class_lyr, " is not a classification")
            }

            
            ## create the new layer of HSU
            if(verbose){ cat("Creating new HSU layer","\n") }
            if(tolower(dist_lyr)=="band"){
                warning("Cell order being used for distance - brk being ingnored")
                brk <- NA
            }
            private$apply_classify(setNames(list(NA,brk),c(class_lyr,dist_lyr)),NULL,layer_name)
            pos_val <- private$find_layer(TRUE)
            brk <- self$get_class_method(layer_name)$cuts[[dist_lyr]]

            ## initialise the model
            model <- list(map = list(hillslope = pos_val[layer_name],
                                     channel = pos_val["channel"],
                                     channel_id = pos_val["channel_id"]))

            ## Create the channel table
            if(verbose){ cat("Creating Channel table","\n") }
            
            channel_id <- raster::raster(pos_val["channel_id"])
            channel_area <- raster::raster(pos_val["channel_id"])
            chn <- raster::shapefile(pos_val["channel"])
            
            model$channel <- chn@data
            model$channel[["precip"]] <- "unknown"
            model$channel[["pet"]] <- "unknown"
            model$channel[["v_ch"]] <- "v_ch_default"            
            tmp <- zonal(channel_area,channel_id,sum)
            model$channel[["area"]] <- 0
            for(ii in 1:nrow(tmp)){
                model$channel[["area"]][model$channel[["id"]]==tmp[ii,"zone"]] <- tmp[ii,"value"]
            }

            rm(channel_area,chn) # remove to clear memory since not used again

            ## Creaate channel flow directions
            model$channel$flow_dir = rep(list(NULL),length(model$channel$id))
            for(rw in 1:nrow(model$channel)){
                jdx <- which( model$channel[['startNode']]==
                              model$channel[['endNode']][rw])
                if( length(jdx) > 0 ){
                    model$channel$flow_dir[[rw]] <- list(idx= model$channel$id[jdx],
                                                         frc= rep(1/length(jdx),length(jdx)))
                }
            }

            ## adapt the hillslope map to have the id as its number
            if(verbose){ cat("Adapting hillslope map to match id","\n") }
            hsu <- raster::raster(pos_val[layer_name])
            hsu <- hsu + max(model$channel$id)
            raster::writeRaster(hsu,pos_val[layer_name],overwrite=TRUE)
            
            ## Create the hillslope table
            if(verbose){ cat("Loading data to create hillslope table","\n") }
            
            ## load the hsu map and work out index of each hsu in it
            hsu <- as.matrix(hsu)
            la <- as.matrix(raster::raster(pos_val["land_area"]))
            gr <- as.matrix(raster::raster(pos_val["gradient"]))
            atb <- as.matrix(raster::raster(pos_val["atb"]))
            cls <- as.matrix(raster::raster(pos_val[class_lyr]))
            dst <- raster::raster(pos_val[dist_lyr])
            if(all(is.na(brk))){
                cdst <- as.matrix(dst)
            }else{
                cdst <- as.matrix(cut(dst,brk,include.lowest=TRUE))
            }
            dst <- as.matrix(dst)
            dem <- as.matrix(raster::raster(pos_val["filled_dem"]))
            channel_id <- as.matrix(channel_id)
            
            ## set all channel cells with no land area to band 0            
            cdst[is.na(cdst) & is.finite(channel_id) & !(la>0)] <- 0
            
            ## distance betweek breaks
            if(tolower(dist_lyr)=="band"){
                dbrk <- rep(mean(private$meta$resolution),max(cdst,na.rm=TRUE))
            }else{
                dbrk <- diff(brk)
            }

            
            if(verbose){ cat("Computing properties of Hillslope HSUs","\n") }
            
            ## unique hsu numbers - go 1,2,3,...
            uhsu <- sort(unique(as.numeric(hsu))) # sort drop sNA lesft by unique

            ## writing into data frames directly is really slow so create the variables outside first
            area <- rep(0,length(uhsu))
            atb_bar <- rep(0,length(uhsu))
            s_bar <- rep(0,length(uhsu))
            class <- rep(NA,length(uhsu))
            delta_x <- rep(NA,length(uhsu))
            band <- rep(NA,length(uhsu))
            width <- rep(0,length(uhsu))
            flow_dir = rep(list(list(idx=integer(0),frc=integer(0))),length(uhsu))

            ## distances and contour lengths
            ## distance between cell centres
            dxy <- rep(sqrt(sum(private$meta$resolution^2)),8)
            dxy[c(2,7)] <- private$meta$resolution[1]; dxy[c(4,5)] <- private$meta$resolution[2]
            dcl <- c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)*mean(private$meta$resolution)
            nr <- nrow(dem); delta <- c(-nr-1,-nr,-nr+1,-1,1,nr-1,nr,nr+1)
            
            ## flags to use in loop
            goes_to_lower_band <- rep(FALSE,length(uhsu))
            min_dst <- tapply(dst,hsu,min)
            min_dst <- min_dst[paste(uhsu)]+delta_min # add wiggle on distance
            set_cnst <- rep(FALSE,length(uhsu))

            idx <- order(dem,na.last=NA)
            n_to_eval <- length(idx)
            it <- 1
            if(verbose){
                print_step <- round(n_to_eval/20)
                next_print <- print_step
            }
            for(ii in idx){
                ##print(ii)
                if(!(la[ii]>0)){ next }
                
                rw <- match(hsu[ii], uhsu) 

                ## These only need to be set once!
                if(!set_cnst[rw]){
                    ## band and class
                    band[rw]  <-  cdst[ii]
                    ## min_dst[rw] <- brk[ cdst[ii] ]
                    class[rw] <- cls[ii]
                    ## length of hillslope - rubbish method...
                    delta_x[rw] <- dbrk[ band[rw] ]
                    set_cnst[rw] <- TRUE
                }

                
                ## area
                area[rw] <- area[rw] + la[ii]
                ## average atanb
                atb_bar[rw] <- atb_bar[rw] + atb[ii]*la[ii]
                ## average gradient
                s_bar[rw] <- s_bar[rw] + gr[ii]*la[ii]

                ## flow direction
                if(is.finite(channel_id[ii])){
                    ## partial channel cell so all flow goes to that channel reach
                    kk <- channel_id[ii]; gcl <- gr[ii]*dxy[1]
                    to_lower_band <- TRUE
                    sw <- dxy[1]
                }else{
                    ## a land cell
                    jj <- ii+delta ## neighbouring cell numbers
                    gcl <- (dem[ii]-dem[jj])*dcl/dxy
                
                    to_lower <- is.finite(gcl) & gcl>0
                    to_use <- to_lower & is.finite(cdst[jj]) & cdst[jj] < cdst[ii]
                    to_lower_band <- TRUE
                    if(!any(to_use) & !goes_to_lower_band[rw]){
                        to_use <- to_lower & is.finite(dst[jj]) & dst[jj] < min_dst[rw]
                        to_lower_band <- FALSE
                    }
                    kk <- hsu[ jj[to_use] ]
                    gcl <- gcl[to_use]
                    sw <- sum(dcl[to_use])
                }
                if(length(kk)>0){
                    if(to_lower_band & !goes_to_lower_band[rw]){
                        ## reset the flow direction
                        flow_dir[[rw]] <- list(idx=integer(0),frc=integer(0))
                        goes_to_lower_band[rw] <- TRUE
                        width[[rw]] <- 0
                    }
                    width[[rw]] <- width[[rw]] + sw

                    ekk <- match(kk,flow_dir[[rw]]$idx)
                    if(any(is.na(ekk))){
                        ## add these to the list
                        flow_dir[[rw]]$idx <-
                            c(flow_dir[[rw]]$idx,kk[is.na(ekk)])
                        flow_dir[[rw]]$frc <-
                            c(flow_dir[[rw]]$frc,rep(0,sum(is.na(ekk))))
                        ekk <- match(kk,flow_dir[[rw]]$idx)
                    }
                    flow_dir[[rw]]$frc[ekk] <-
                        flow_dir[[rw]]$frc[ekk] + gcl
               }
                
                ## verbose output here
                if(it >= next_print){
                    cat(round(100*it / n_to_eval,1),
                        "% complete","\n")
                    next_print <- next_print+print_step
                }
                
                it <- it+1
            }
            ## check every HSU has its constants
            if(!all(set_cnst)){
                stop("HSU with the following map_id cannot be generated: \n",
                     paste(uhsu[!set_cnst],collapse=", "))
            }

            ## rescale by area
            atb_bar <-  atb_bar/ area
            s_bar <-  s_bar/ area

            ## warn if any HSU does not go to lower band
            if(!all(goes_to_lower_band)){
                warning("The following HSUs with the following map_id do not drain to a lower distance class: \n",
                     paste(uhsu[!goes_to_lower_band],collapse=", "))
            }
            
            ## tidy up flow
            for(rw in 1:length(flow_dir)){
                if( length(flow_dir[[rw]]$frc) != length(flow_dir[[rw]]$idx) ){
                    stop("Flow length and fractions are different lengths for hsu: ",uhsu[rw])
                }
                if( length(flow_dir[[rw]]$frc) == 0 ){
                    stop("Flow length and fractions are 0 for hsu: ",uhsu[rw])
                }
                flow_dir[[rw]]$frc <-
                    flow_dir[[rw]]$frc / sum(flow_dir[[rw]]$frc)
            }

            ## some checks
            if(!all(set_cnst)){
                stop("HSU with the following map_id cannot be generated: \n",
                     paste(uhsu[!set_cnst],collapse=", "))
            }
            if(!all(set_cnst)){
                stop("HSU with the following map_id cannot be generated: \n",
                     paste(uhsu[!set_cnst],collapse=", "))
            }
            ## initiailise the hillslope table
            model$hillslope <- data.frame(
                id = as.integer(uhsu),
                area = area,
                atb_bar = atb_bar,
                s_bar = s_bar,
                band = band,
                class = class,
                delta_x = delta_x,
                width = width,
                precip="unknown",
                pet="unknown",
                r_sfmax="r_sfmax_default",
                s_rzmax="s_rzmax_default",
                s_rz0="s_rz0_default",
                ln_t0="ln_t0_default",
                m="m_default",
                t_d="t_d_default",
                c_sf="c_sf_default",
                stringsAsFactors=FALSE
            )
            model$hillslope$flow_dir = flow_dir
            
            ## parameter values
            if(verbose){ cat("Setting default parameter values","\n") }
            model$param <- c(r_sfmax_default=Inf,
                             s_rzmax_default=0.05,
                             s_rz0_default=0.75,
                             ln_t0_default=-2,
                             m_default=0.04,
                             t_d_default=2*60*60,
                             c_sf_default=0.1,
                             v_ch_default=1)
            
            ## ############################################
            ## Add gauges at all outlets from river network
            ## ############################################
            if(verbose){ cat("Adding gauges at the outlets","\n") }
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
            if(verbose){ cat("Adding a blank point_inflow table","\n") }
            model$point_inflow <- data.frame(
                name = character(0),
                id = integer(0),
                fraction = numeric(0),
                stringsAsFactors=FALSE
            )
            
            ## ##################################
            ## Add diffuse inflow table
            if(verbose){ cat("Adding a blank diffuse_inflow table","\n") }
            model$diffuse_inflow <- data.frame(
                name = character(0),
                id = integer(0),
                stringsAsFactors=FALSE
            )

            ## add model to record
            
            fn <- file.path(private$wdir,paste0(layer_name,".rds"))
            saveRDS(model,fn)
            private$meta$layers[[layer_name]]$model_file <- fn
            private$meta$layers[[layer_name]]$type <- "model"
            private$write_meta()
            
        },
        apply_input_frac = function(hsu_layer,input_layer,hsu_label,input_label){
            
            rq <- c(hsu_layer,input_layer,"channel_id")
            pos_val <- private$find_layer(TRUE)

            has_rq <- rq %in% names(pos_val)
            if(!all(has_rq)){
                stop(paste(c("Missing layers:",rq[!has_rq],sep="\n")))
            }
            if(!(private$meta$layers[[hsu_layer]]$type %in% "model")){
                stop(hsu_layer, " is not a model layer")
            }
            if(!(private$meta$layers[[input_layer]]$type %in% c("user","classification"))){
                stop(input_layer, " is not a user defined or classification layer")
            }
            ## load data
            hsu_grid <- raster::raster(pos_val[hsu_layer])
            chn_grid <- raster::raster(pos_val["channel_id"])
            input_grid <- raster::raster(pos_val[input_layer])
            ## make the matrix
            out <- rbind(raster::crosstab(chn_grid,input_grid),
                         raster::crosstab(hsu_grid,input_grid))
            ## relabel rows and columns
            ##browser()
            rownames(out) <- paste0(hsu_label,rownames(out))
            colnames(out) <- paste0(input_label,colnames(out))
            ## check for NA
            if(!all(is.finite(out))){
                stop("Non finite values in the cross tabulation \n",
                     "Check the input is defined for all hsu cells")
            }
            ## sweep out the rowSums to get scaled to 1
            out <- sweep(out,1,rowSums(out),"/")
            return(out)
        }
        
    )
    )
