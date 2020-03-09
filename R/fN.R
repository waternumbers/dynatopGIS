#' Function to give neighbours of cells
#' @description Returns index of neighbouring cells and possibly some properties
#'
#' @param k cell index (only single values accepted)
#' @param nr number fo rows in original raster
#' @param nc number fo columnc in original raster
#' @param dx resolution of raster grid (presumed square)
#'
#' @return The index of neighbouring cells. If dx is a single finite value then values or contour length, and distance to neighbours is returned.
#'
#' @details This is in effect a cut down version of raster adjacent for use with the vectors exctraced by raster::getValues. See example below.
#'
#' @examples
#' ## create a test raster
#' rst <- raster::raster(matrix(1:200,20,10))
#' nc <- ncol(rst)
#' nr <- nrow(rst)
#' ## specify some locations
#' kvec <- c(1, # top left
#'           floor(nc/2), # top row middle column
#'           nc, # top right
#'           (nc*floor(nr/2))+1, # first column, middle row
#'           (nc*floor(nr/2))+floor(nc/2), # middle column middle row
#'           (nc*floor(nr/2))+nc, # last column, middle row
#'           (nc*(nr-1))+1, # bottom left
#'           (nc*nr)- floor(nc/2), # bottom row middle
#'           (nc*nr)) # bottom right
#'
#' ## test matching of locations
#' for(k in kvec){
#'     ra <- sort( raster::adjacent(rst,k,direction=8,pairs=FALSE) )
#'     fn <- sort( fN(k,nrow(rst),ncol(rst)) )
#'     if( length(ra)!=length(fn) || any( ra!=fn ) ){
#'         stop(paste("fN does not match raster adjacent for value",k))
#'     }
#' }
#' @export
fN <- function(k,nr,nc,dx=NA){
    if( k > nc*nr | k < 1){ return(numeric(0)) }
    pflg <- length(dx)==1 && is.finite(dx)

    ## find i,j location
    j <- ((k-1)%/%nc) +1 # row
    i <- k - (j-1)*nc #column
    m <- matrix(c(i-1,i,i+1,i-1,i+1,i-1,i,i+1,
                  j-1,j-1,j-1,j,j,j+1,j+1,j+1),8,2)
    ## compute properties if required
    if(pflg){
        dxd <- dx*sqrt(2)
        p <- matrix( c(dxd,dx,dxd,dx,dx,dxd,dx,dxd, # distance
                       dx*c(0.35,0.5,0.35,0.5,0.5,0.35,0.5,0.35)), # from Quinn
                       8,2)
        ## Ocatagon based countout length is
        ##rep( dx /(1+sqrt(2)),8))
        m <- cbind(m,p)
    }
    ## trim to within bounds
    m <- m[ m[,1]>0 & m[,1]<=nc & m[,2]>0 & m[,2]<=nr ,]

    if( pflg ){
        return( list(idx = (m[,2]-1)*nc + m[,1],
                     dx = m[,3],
                     cl = m[,4]) )
    }else{
        return((m[,2]-1)*nc + m[,1])
    }
}
