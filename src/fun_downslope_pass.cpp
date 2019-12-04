#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for passing up the catchments from river nodes
//' 
//' @param dem Digital elevation model as a vector
//' @param order initial order values as a vector, internally starts at one and move upstream
//' @param nc number of columns in the matrix
//' 
//' @return a list with the filled dem and order
//'
// [[Rcpp::export]]
List fun_downslope_pass(NumericVector dem, IntegerVector order,
			IntegerVector seq,
			IntegerVector offset,
			NumericVector area, NumericVector dx,
			NumericVector cl){

  NumericVector upslope_area = area;
  NumericVector contour_length(dem.length(),NA_REAL);
  NumericVector gradient(dem.length(),NA_REAL);
  NumericVector atanb(dem.length(),NA_REAL);
  
  IntegerVector ngh(offset.length());
  
  
  // distance, area and contour length
  
  // loop down sequence
  for(int s=0;s < seq.length(); s++){
    
    NumericVector w(8,0.0); // initialise weight
    double sum_w = 0; // number of weights
    int num_w = 0;

    int i = seq(s);
    
    contour_length(i) = 0;
	
    ngh = offset + i;
    LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);

    if( order(i) > 1 ){ // everything except last order
      for(int j=0;j<ngh.length();j++){
	if( in_range(j) ){
	  if( !(NumericVector::is_na(dem(ngh(j)))) &&
	      dem(ngh(j)) < dem(i) ){
	    w(j) = (dem(i) - dem(ngh(j))) / dx(j);
	    sum_w += w(j);
	    num_w += 1;
	    contour_length(i) += cl(j);
	  }
	}
      }
      gradient(i) = (sum_w/num_w);
      atanb(i) = log( (upslope_area(i)/contour_length(i)) / gradient(i) );
      for(int j=0;j<8;j++){
	if( w(j) > 0 ){
	  upslope_area(ngh(j)) = upslope_area(ngh(j)) +
	    upslope_area(i)*(w(j)/sum_w);	
	}
      }
    }else{
      // for bottom order take upslope garient to work out atanb
      for(int j=0;j<8;j++){
	if( in_range(j) ){
	  if( !(NumericVector::is_na(dem(ngh(j)))) &&
	      dem(ngh(j)) > dem(i) ){
	    w(j) = (dem(ngh(j)) - dem(i)) / dx(j);
	    sum_w += w(j);
	    num_w += 1;
	  }
	}
      }
      gradient(i) = sum_w/num_w;
      contour_length(i) = cl(9);
      atanb(i) = log( (upslope_area(i)/contour_length(i)) / gradient(i) );
    }
  }
  
  // create output
  List out=List::create(Named("upslope_area")=upslope_area,
			Named("contour_length")=contour_length,
			Named("gradient")=gradient,
			Named("atanb")=atanb);
  
  return out;
}
