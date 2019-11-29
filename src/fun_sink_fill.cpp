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
NumericVector fun_sink_fill(NumericVector dem, LogicalVector is_channel,
			    NumericVector delta, int nc, int max_iter){

  
  // int n_sink = 0;
  // int n_finite = 0;
  // int n_iter = 0;
  // int ngh(8);
  // int n_lower=0;
  // double min_ngh = 1e32;
  
  // LogicalVector is_sink(dem.length(),true);
  
  // double na_test_val = -10000; // test the for NAN against this - if NAN will return false
  
  // // loop
  // n_iter = 1;
  // while( (n_iter < max_iter) ){ //any(is_sink) ){ //& (n_iter < max_iter) ){
  //   n_sink=0;
  //   for(int i=0;i < dem.length(); i++){
  //     if( is_sink(i) ){
  // 	if( (dem(i) > na_test_val) & not is_channel(i) ){
  // 	  // work out index of neighbours
  // 	  ngh(0) = i-nc-1;
  // 	  ngh(1) = i-nc;
  // 	  ngh(2) = i-nc+1;
  // 	  ngh(3) = i-1;
  // 	  ngh(4) = i+1;
  // 	  ngh(5) = i+nc-1;
  // 	  ngh(6) = i+nc;
  // 	  ngh(7) = i+nc+1;
	  
  // 	  // set the lowest neighbour value and number of finite neighbours
  // 	  min_ngh = 1e32;
  // 	  n_finite = 0;
  // 	  n_lower = 0;
  // 	  // loop neighbours
  // 	  for(int j=0; j < 8; j++){
  // 	    if( ngh(j) < dem.length() &&
  // 		ngh(j) > -1 &&
  // 		dem(ngh(j)) > na_test_val &&
  // 		dem(ngh(j)) < dem(i) ){
  // 	      // increase number of finite values
  // 	      n_finite = n_finite + 1;
  // 	      // compute new lowest neighbour
  // 	      n_lower = n_lower + 1;
  // 	    }
  // 	  }
	  
  // 	  // work out if changed
	  
  // 	  // set new value and neighbours to is_sink=true
	  
  // 	  if( (n_finite < 8) && (n_lower==0) ){
  // 	    // potential edge drain set order to 1
  // 	    //order(i) = 1;
  // 	  }
  // 	}else{
  // 	  //can_eval(i)=false;
  // 	}
  // 	Rcpp::checkUserInterrupt();
  // 	Rcout << "The number of sinks handled in iteration " << n_iter << " is " << n_sink << std::endl;
	
  //     }
  //   }
  // }
  
  return dem; 
}
