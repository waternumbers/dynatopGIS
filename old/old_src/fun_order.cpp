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
List fun_increase_order(NumericVector dem, IntegerVector order, int nc ){

  List out=List::create(Named(dem)=dem,Named(order)=order);
    
  IntegerVector ngh_offset = {-nc-1,-nc,-nc+1,-1,1,nc-1,nc,nc+1};
  IntegerVector ngh(8);
  IntegerVector n_lower(dem.length(),10);

  
  //  double na_test_val = -10000; // test the for NAN against this - if NAN will return false
  //int cnt = 0, dp=0, n_finite=0;

  //double delta_min = 0.00;
  
  // work out number of lower cells
  for(int i=0;i < dem.length(); i++){
    if( !(NumericVector::is_na(dem(i))) ){
      // neighbours
      ngh = ngh_offset + i;
      
      IntegerVector good_ngh = ngh[!is_na(dem(ngh)) & (ngh<dem.length()) & (ngh>-1)];
      //NumericVector good_dem = dem(good_ngh);

      n_lower(i) = sum( dem(good_ngh) < dem(i) );
      if( (good_ngh.length() < 8) && (n_lower(i)==0) ){
	// potential edge drain set order to 1
	order(i) = 1;
      }
      if( (good_ngh.length() == 8) && (n_lower(i)==0) ){
	// potential sink drain set order to 1
	Rcout << "Sink at cell number " << i << "\n";
      }
      
    }
  }
  
  // loop depths
  dp = 1;
  cnt = 1;
  while( (cnt > 0) && (dp<max_depth) ){
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
	      (dem(ngh(j))+delta_min) > dem(i) ){
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
    Rcout << "Evaluated depth " << dp << " and set " << cnt << " cells " <<"\n";
  }
  
  return order; 
}
