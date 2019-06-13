// Sink filling algorithm applied using Rcpp //

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fun_sink_fill(NumericMatrix dem, LogicalMatrix is_channel,
			    NumericMatrix delta, int max_iter){

  int n_sink = 0;
  int n_finite = 0;
  int niter = 0;
  LogicalMatrix can_eval(dem.nrow(),dem.ncol());

  double na_test_val = -10000; // test the for NAN against this - if NAN will return false
    
  // determine if we can fill a cell if it is a sink
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      // set to can't be filled
      can_eval(i,j) = false;
      // to fill it there must be 9 finite values in the block centred on i,j
      n_finite = 0;
      for(int ii=-1;ii<2;ii++){
	for(int jj=-1;jj<2;jj++){
	  if( i+ii > -1 && i+ii < dem.nrow() && 
	      j+jj > -1 && j+jj < dem.ncol() ){
	    if( dem(i+ii,j+jj) > na_test_val ){
	      n_finite = n_finite + 1;
	    }
	  }
	}
      }
      // see if there are 9 finite values
      if(n_finite==9){
	can_eval(i,j) = true;
      }
    }
  }
  
  // loop to fill
  n_sink = 1;
  while( n_sink > 0 && niter < max_iter){
    n_sink = 0;
    // remember cpp is 0 base
    for(int i=0;i < dem.nrow(); i++){
      for(int j=0; j < dem.ncol(); j++){
	// check it can be filled and isn't a channel
	if( can_eval(i,j)==true && is_channel(i,j)==false){
	  // presume it is a sink
	  bool is_sink = true;
	  double lowest_neighbour = 1e32;
	  
	  // check its status
	  for(int ii=-1;ii<2;ii++){
	    for(int jj=-1;jj<2;jj++){
	      if( i+ii > -1 && i+ii < dem.nrow() && 
		  j+jj > -1 && j+jj < dem.ncol() &&
		  !(ii == 0 && jj == 0) ){
		// value that dem(i,j) but excedd if the gradient in greater then min required
		double neighbour = dem(i+ii,j+jj) + delta(ii+1,jj+1);
		// see if this neighbour is low enough
		if( neighbour <= dem(i,j) ){
		  is_sink = false;
		}
		// if it is a sink then update the possible value
		if( is_sink==true && neighbour < lowest_neighbour ){
		  lowest_neighbour = neighbour;
		}
	      }
	    }
	  }
	  
	  // if it is still a sink
	  if(is_sink==true){
	    dem(i,j) = lowest_neighbour;
	    n_sink = n_sink + 1 ;
	  }
	}
	Rcpp::checkUserInterrupt();
      }
    }
    niter = niter + 1;
    Rcpp::checkUserInterrupt();
    Rcout << "The number of sinks handled in iteration " << niter << " is " << n_sink << std::endl;
    
  }
  
  return dem; 
}
