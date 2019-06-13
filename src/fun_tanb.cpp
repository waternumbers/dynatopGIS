// Sink filling algorithm applied using Rcpp //

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fun_tanb(NumericMatrix dem, LogicalMatrix is_channel,
		       NumericMatrix dist){
  
  NumericMatrix tanb(dem.nrow(),dem.ncol());
  double n_above, s_above, n_below, s_below;
  double tb;
  int n_finite;
  
  double na_test_val = -10000; // test for NAN against this - if NAN will return false
    

  // determine if we can route a cell area downstream
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      if( dem(i,j) > na_test_val ){
	// initialise sums
	n_above=0;
	s_above=0;
	n_below=0;
	s_below=0;
	n_finite = 0;
	// loop neighbours
	for(int ii=-1;ii<2;ii++){
	  for(int jj=-1;jj<2;jj++){
	    if( i+ii > -1 && i+ii < dem.nrow() && 
		j+jj > -1 && j+jj < dem.ncol() &&
		!(jj==0 && ii==0)){
	      n_finite = n_finite + 1;
	      if( dem(i+ii,j+jj) > na_test_val ){
		//n_finite = n_finite + 1;
		// tanb
		tb = (dem(i,j) - dem(i+ii,j+jj))/dist(ii+1,jj+1);
		if( tb >= 0 ){
		  n_below = n_below+1;
		  s_below = s_below+tb;
		}else{
		  n_above = n_above + 1;
		  s_above = s_above - tb;
		}
	      }
	    }
	  }
	}
	//Rcout << n_above + n_below  << " " << n_finite << std::endl;
	// see if we should use above or below gradient
	if( is_channel(i,j)==true || (n_above+n_below) < 8 ){
	  //Rcout << "above"  << std::endl;
	  tanb(i,j) = s_above / n_above;
	}else{
	  //Rcout << "below"  << std::endl;
	  tanb(i,j) = s_below / n_below;
	}
	
      }
      Rcpp::checkUserInterrupt();
    }
  }
  
  return tanb;
}
