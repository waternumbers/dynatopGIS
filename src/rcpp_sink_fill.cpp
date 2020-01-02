#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for passing up the catchments from river nodes
//' 
//' @param dem Digital elevation model as a vector
//' @param channel_id UID of channel in the pixel (id any) as a vector
//' @param offset difference between index of neighbours and current cell - clockwise from top left
//' 
//' @return a list with the filled dem
//'
// [[Rcpp::export]]
NumericVector rcpp_sink_fill(NumericVector dem,
			     IntegerVector channel_id,
			     IntegerVector offset){

  // store locations of sinks and lowest finite neighbour
  NumericVector filled_dem(dem.length(),NA_REAL);
  NumericVector lowest_ngh(dem.length(),NA_REAL);
  

  // look at all cells
  for(int i=0;i < dem.length(); i++){
    if( !(NumericVector::is_na(dem(i))) ){
      // then should be finite
      filled_dem(i) = dem(i);
      
      // work out lowest valid neighbour
      double min_ngh=R_PosInf;
      IntegerVector ngh = offset + i; // all possible neighbours
      LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);
      int n_finite = 0;
      for(int j=0;j<ngh.length();j++){
	if( in_range(j) ){
	  int jdx = ngh(j);
	  if( !(NumericVector::is_na(dem(jdx))) ){
	    n_finite += 1;
	    if( !(dem(jdx)==R_NegInf) & (dem(jdx) < min_ngh) ){
	      min_ngh = dem(jdx);
	    }
	  }
	}
      }
      Rcout << i<< " " << n_finite << " " << min_ngh << "\n";
      if( n_finite==8 ){
	lowest_ngh(i) = min_ngh;
      }else{
	// on edge so set lowest nieghbour very low
	lowest_ngh(i) = R_NegInf;
      }
    }
    Rcpp::checkUserInterrupt();
  }
  
  // loop possible sinks
  int iter=0;
  while(iter<100){
    // find possible sink with lowest neighbour
    int idx=-1;
    double idx_min=R_PosInf;   
    for(int j=0;j<filled_dem.length();j++){
      // Rcout << j << " " << idx_min << " " << lowest_ngh(j) << " " << filled_dem(j) << "\n";
      if( (lowest_ngh(j) < idx_min) & (lowest_ngh(j) > filled_dem(j)) ){
	idx = j;
	idx_min = lowest_ngh(j);
      }
    }
    Rcout << "idx " << idx << " " << idx_min << "\n";

    if(idx < 0){
      // if no possible value then set to large numebr to exit
      iter = INT_MAX;
      Rcout << "stopping" << "\n";
    }else{
      // process the sink
      IntegerVector ngh = offset + idx; // all possible neighbours
      LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);
      
      // populate
      double sum_dem=0;
      double n_dem =0;
      for(int j=0;j<ngh.length();j++){
	int jdx = ngh(j);
	if( !(filled_dem(jdx)==R_NegInf) ){
	  sum_dem += filled_dem(jdx);
	  n_dem += 1;
	}
      }
      //Rcout << "refill " << sum_dem << " " << n_dem << "\n";
      filled_dem(idx) = sum_dem / n_dem;
      
      // change lowest neighbour for neighbours
      for(int j=0;j<ngh.length();j++){
	int jdx = ngh(j);
	if( lowest_ngh(jdx) > filled_dem(idx) ){
	  lowest_ngh(jdx) = filled_dem(idx);
	}
      }
      
      // increment iter
      iter += 1;
      //Rcout << iter << " " << n_sink << "\n";
    }
    Rcpp::checkUserInterrupt();
  }
  
  return filled_dem; 
}
