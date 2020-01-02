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
NumericVector fun_sink_fill(NumericVector dem,
		   IntegerVector channel_id,
		   IntegerVector offset,
		   NumericVector delta){

  // store locations of sinks and lowest finite neighbour
  LogicalVector is_sink(dem.length(),false);
  NumericVector min_ngh(dem.length(),NA_REAL);
  int n_sink = 0;
  
  // work out number of lower cells
  for(int i=0;i < dem.length(); i++){
    if( !(NumericVector::is_na(dem(i))) ){
      // has a finite dem value - check if a sink
      
      IntegerVector ngh = offset + i; // all possible neighbours
      LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);

      int n_finite=0;
      double min_dem=R_PosInf;
      
      for(int j=0;j<ngh.length();j++){
	if( in_range(j) ){
	  int jdx = ngh(j);

	  if( !(NumericVector::is_na(dem(jdx))) ){
	    n_finite += 1;
	    if( !(dem(jdx)==R_NegInf) & (dem(jdx) < min_dem) ){
	      min_dem = dem(jdx);
	    }
	  }
	}
      }
      min_ngh(i) = min_dem;
      if( (n_finite == 8) & (min_ngh(i)>dem(i)) ){
	  //then a sink
	  is_sink(i) = true;
	  n_sink += n_sink;
      }
    }
  }

  // fix sinks
  for(int i=0;i < n_sink; i++){
    // find sink with lowest value
    int idx=-1;
    double idx_min=R_PosInf;
    
    for(int j=0;j<is_sink.length();j++){
      if( !NumericVector::is_na(min_ngh(j)) & (min_ngh(j) < idx_min) ){
	idx = j;
	idx_min = min_ngh(j);
      }
    }

    if( idx>-1 ){
      // then populate sink
      IntegerVector ngh = offset + idx; // all possible neighbours

      double sum_dem=0;
      double n_dem =0;
      for(int j=0;j<ngh.length();j++){
	int jdx = ngh(j);
	if( !(dem(jdx)==R_NegInf) ){
	  sum_dem += dem(jdx);
	  n_dem += 1;
	}
      }
      dem(idx) = sum_dem / n_dem;

      // change mins for neighbours
      for(int j=0;j<ngh.length();j++){
	int jdx = ngh(j);
	if( is_sink(jdx) and min_ngh(jdx) > dem(idx) ){
	  min_ngh(jdx) = dem(idx);
	}
      }
    }
  }
  
  return dem; 
}
