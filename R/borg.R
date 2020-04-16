#' update a dynatopGIS class object
#' 
#' @description Update a dynatopGIS object while keeping layers and channel
#'
#' @param old_object the current dynatopGIS object
#' @param force force conversion even to an older version
#' 
#' @return a dynatopGIS updated to the current version of the methods
#' @details Setting is_class to TRUE and not providing a layer_name returns the overall classification
#' @export
borg <- function(old_object,force=FALSE){
    if(!("dynatopGIS" %in% class(old_object))){
        stop("Not a dynatopGIS class object")
    }

    ## dem and create object
    new_object <- dynatopGIS$new(old_object$get_layer("dem"))

    if( !force ){
        if( new_object$get_version() < old_object$get_version() ){
            stop(paste("This is converting to an older version.",
                       "Set the force flag in the call if you really want to do this",
                       sep="\n"))
        }
    }
    
    ## add channel
    new_object$add_channel(old_object$get_channel())

    ## copy other layers
    for(nm in old_object$get_layer()){
        new_object$add_layer(old_object$get_layer(nm),nm,FALSE)
    }

    ## get classification
    cls <- old_object$get_class_method()
    if(length(cls)>0){
        new_object$classify(cuts=cls$cuts,burns=cls$burns)
    }

    return(new_object)
}
