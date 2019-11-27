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
IntegerVector fun_increase_order(NumericVector dem, IntegerVector order, int nc ){

  IntegerVector ngh(8), j_ngh(8);
  IntegerVector n_lower(dem.length(),10);
  
  double na_test_val = -10000; // test the for NAN against this - if NAN will return false
  int cnt = 0, dp=0, n_finite=0;

  // work out number of lower cells
  for(int i=0;i < dem.length(); i++){
    if( dem(i) > na_test_val ){
      // work out index of neighbours
      ngh(0) = i-nc-1;
      ngh(1) = i-nc;
      ngh(2) = i-nc+1;
      ngh(3) = i-1;
      ngh(4) = i+1;
      ngh(5) = i+nc-1;
      ngh(6) = i+nc;
      ngh(7) = i+nc+1;

      // set number of lower neighbours to 0
      n_lower(i) = 0;
      n_finite = 0;
      // loop neighbours
      for(int j=0; j < 8; j++){
	if( ngh(j) < dem.length() &&
	    ngh(j) > -1 &&
	    dem(ngh(j)) > na_test_val &&
	    dem(ngh(j)) < dem(i) ){
	  n_lower(i) = n_lower(i) + 1;
	}
      }
      if( (n_finite < 8) && (n_lower(i)==0) ){
	// potential edge drain set order to 1
	order(i) = 1;
      }
      
    }
  }
  
  // loop depths
  dp = 1;
  cnt = 1;
  while( (cnt > 0) && (dp<10000) ){
    //    for(int dp=1;dp<10000;dp++){
    cnt = 0;
    // loop cells
    for(int i=0;i < dem.length(); i++){
      
      if( order(i) == dp ){
	
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
	  if( ngh(j) < dem.length() &&
	      ngh(j) > -1 &&
	      dem(ngh(j)) > na_test_val &&
	      dem(ngh(j)) > dem(i) ){
	    // cell has finite dem and is higher then cell i
	    // take one of the n_lower value
	    n_lower(ngh(j)) = n_lower(ngh(j)) - 1;

	    // if n_lower is 0 then set to next depth
	    if(n_lower(ngh(j)) == 0){
	      order(ngh(j)) = dp + 1;
	      cnt = cnt + 1;
	    }
	  }
	}
      }
    }
    dp = dp +1;
    Rcout << "Evaluated depth " << dp+1 << " and set " << cnt << " cells " <<"\n";
  }
  
  return order; 
}
