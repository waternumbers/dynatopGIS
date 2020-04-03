#' update a dynatopGIS class object
#' 
#' @description Update a dynatopGIS object while keeping layers and channel
#'
#' @param old_object the current dynatopGIS object
#' @param is_class logical if the class based on this layer should be returned
#' @return a RasterLayer of the requested information
#' @details Setting is_class to TRUE and not providing a layer_name returns the overall classification
#' @export
borg <- function(old_object){
    if(!("dynatopGIS" %in% class(old_object))){
        stop("Not a dynatopGIS class object")
    }

    ## dem and create object
    new_object <- dynatopGIS$new(old_object$get_layer("dem"))

    ## add channel
    new_object$add_channel(old_object$get_channel())

    ## copy other layers
    for(nm in old_object$list_layers()){
        new_object$add_layer(old_object$get_layer(nm),nm,FALSE)
    }

    ## get classification
    cls <- old_object$get_class()
    new_object$classify(cuts=cls$cuts,burns=cls$burns)

    return(new_object)
}
