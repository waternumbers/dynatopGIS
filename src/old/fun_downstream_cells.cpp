#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for working out downstream cells
//' 
//' @param dem Digital elevation model as vector
//' @param is_channel logical vector of channel pixels
//'
//' @return matrix containing fil
//'
// [[Rcpp::export]]
IntegerMatrix fun_edges(NumericVector dem, LogicalVector channel, int nc ){

  IntegerVector ngh(8);
  IntegerMatrix edges(dem.length()*8,2);
  Rcout << "initial declarations\n";
  double na_test_val = -10000; // test the for NAN against this - if NAN will return false
  long cnt = 0;

  //double delta_min = 0.00;
  
  for(long i=0;i < dem.length(); i++){
    Rcout << i << " " << cnt <<"\n";
    if( (dem(i) > na_test_val) and not channel(i) ){
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
	    dem(ngh(j)) < dem(i) ){
	  edges(cnt,1)=i;
	  edges(cnt,2)=ngh(j);
	  cnt=cnt+1;
	}
      }
    }
  }
  return edges; 
}
