// Sink filling algorithm applied using Rcpp //

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fun_tanb(NumericMatrix dem, LogicalMatrix is_channel,
		       NumericMatrix dist){
  
  NumericMatrix tanb(dem.nrow(),dem.ncol());
  numeric n_above, s_above, n_below, s_below;
  numeric g;
  
  numeric na_test_val = -10000; // test for NAN against this - if NAN will return false
    

  // determine if we can route a cell area downstream
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      if( dem(i,j) > na_test_val ){
	// initialise sums
	n_above=0;
	s_above=0;
	n_below=0;
	s_below=0;

	// loop neighbours
	for(int ii=-1;ii<2;ii++){
	  for(int jj=-1;jj<2;jj++){
	    if( i+ii > -1 && i+ii < dem.nrow() && 
		j+jj > -1 && j+jj < dem.ncol() &&
		jj !=0 && ii !=0 ){
	      if( dem(i+ii,j+jj) > na_test_val ){
		// tanb
		tb = (dem(i,j) - dem(i+ii,j+jj))/dist(ii+1,jj+1);
		if( tb >= 0 ){
		  n_below = n_below+1;
		  s_below = s_below+g;
		}else{
		  n_above = n_above+1;
		  s_above = s_above-g;
		}
	      }
	    }
	  }
	}
	// see if we shoould use above or below gradient
	if( is_channel(i,j) | n_above+n_below < 8 ){
	  tanb(i,j) = s_above / n_above;
	}
      }
      Rcpp::checkUserInterrupt();
    }
  }
  
  return tanb;
}
