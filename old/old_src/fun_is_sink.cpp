#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for filling of sinks
//' 
//' @param dem Digital elevation model
//' @param is_channel TRUE is a channel pixel
//' @param delta 3x3 matrix of minimum evelation drop to each adjacent pixel
//' @param max_iter maximum number of iterations
//'
//' @return matrix containing filled dem
//'
// [[Rcpp::export]]
LogicalVector fun_is_sink(NumericVector dem, int nc ){

  LogicalVector is_sink(dem.length());
  IntegerVector ngh(8);
    
  double na_test_val = -10000; // test the for NAN against this - if NAN will return false

  // determine if a cell if it is a sink
  for(int i=0;i < dem.length(); i++){
    // presume not a sink
    is_sink(i) = false;
    if( dem(i) > na_test_val ){
      // finite value so initialy presume it is a sink
      is_sink(i) = true;

      // work out index of neighbours
      ngh(0) = i-nc-1;
      ngh(1) = i-nc;
      ngh(2) = i-nc+1;
      ngh(3) = i-1;
      ngh(4) = i+1;
      ngh(5) = i+nc-1;
      ngh(6) = i+nc;
      ngh(7) = i+nc+1;

      // loop neighbours
      for(int j=0; j < 8; j++){
	if( ngh(j)< dem.length() &&
	    ngh(j)>-1 ){
	  int k = ngh(j);
	  if( dem(k) > na_test_val ){
	    // then a finite value so test height
	    if( dem(k) < dem(i) ){
	      is_sink(i) = false;
	    }
	  }else{
	    // na value so do know if a sink..
	    is_sink(i)=false;
	  }
	}
      }
    }
  }
    
  return is_sink; 
}
