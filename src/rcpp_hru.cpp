#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for computation of redistribution matrices
//' 
//' @param dem Digital elevation model
//' @param grad average downslope gradient
//' @param land_area Area of land surface in pixel
//' @param channel_area Area of channel surface in pixel
//' @param channel_id UID of channel present in the pixel
//' @param hillslope_id hillslope class of each pixel
//' @param offset - difference between cell index of adjacent cells and current cell index - clockwise from top left
//' @param dx distance between cell centres - from top left in clockwise direction
//' @param cl contour length - from top left in a clockwise direction. The 9th value is used for cells split beteen land and channel
//' @param max_index maximum value of the hillslope and channel id's
//' @return list of hillslope and channel properties
//'
// [[Rcpp::export]]
List rcpp_hru(NumericVector dem,
	      NumericVector grad,
	      NumericVector land_area,
	      NumericVector channel_area,
	      IntegerVector channel_id,
	      IntegerVector hillslope_id,
	      IntegerVector offset,
	      NumericVector dx,
	      NumericVector cl,
	      int max_index){
  
  // initialise the output - default filled with zeros
  NumericVector area(max_index);
  NumericVector av_grad(max_index);
  NumericMatrix W(max_index,max_index);

  for( int i=0; i<dem.length(); i++){
    if( !(IntegerVector::is_na(hillslope_id(i))) ){
      // then a hillslope element
      area( hillslope_id(i) ) += land_area(i); // add area
      av_grad( hillslope_id(i) ) += land_area(i)*grad(i); //weighted sum of average gradient
      // redistribute downstream
      if( !(IntegerVector::is_na(channel_id(i))) ){
	// then part cell with channel - all flow to channel
	W( channel_id(i), hillslope_id(i) ) += land_area(i);
      }else{
	// work out downslope gradients
	IntegerVector ngh = offset + i;
	LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);
	NumericVector gl(8,0.0); // gradient
	for(int j=0;j<ngh.length();j++){	
	  if( in_range(j) ){
	    if( !(NumericVector::is_na(dem(ngh(j)))) &&
		(land_area(i) > 0) &&
		dem(ngh(j)) < dem(i) ){
	      gl(j) = cl(j)*( (dem(i) - dem(ngh(j))) / dx(j) ); //tan(beta)*L
	    }
	  }
	}
	double sum_gl = sum(gl);
	// distribute downstream
	for(int j=0;j<ngh.length();j++){	
	  if( in_range(j) ){
	    // drain to hillslope if next downstream has hillslope_id
	    if( !(IntegerVector::is_na(hillslope_id(ngh(j)))) &&
		dem(ngh(j)) < dem(i) ){
	      W( hillslope_id(ngh(j)), hillslope_id(i) ) += gl(j)*land_area(i)/sum_gl;
	    }else{
	      // if not got a hillslope_id then drain to channel if available
	      if( !(IntegerVector::is_na(channel_id(ngh(j)))) &&
		  dem(ngh(j)) < dem(i) ){
		W( channel_id(ngh(j)), hillslope_id(i) ) += gl(j)*land_area(i)/sum_gl;
	      }
	    }
	  }
	}
      }
    }
    // handle if has a channel_id
    if( !(IntegerVector::is_na(channel_id(i))) ) {
      // then a hillslope element
      area( channel_id(i) ) += channel_area(i); // add area
    }
  }
  
  return List::create(Named("area") = area,
		      Named("av_grad") = av_grad,
		      Named("W") = W);
}
