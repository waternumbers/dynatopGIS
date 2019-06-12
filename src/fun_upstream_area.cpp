// Sink filling algorithm applied using Rcpp //

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fun_upslope_area(NumericMatrix dem, NumericMatrix area, LogicalMatrix is_channel,
				NumericMatrix dist, int max_iter){
  
  //NumericMatrix upstream_area(dem.nrow(),dem.ncol());
  LogicalMatrix can_eval(dem.nrow(),dem.ncol());
  IntegerMatrix n_above(dem.nrow(),dem.ncol());
  int niter = 0;
  int n_propergate = 0;

  NumericMatrix W(3,3);
  numeric sW;
    
  
  numeric na_test_val = -10000; // test for NAN against this - if NAN will return false
    

  // determine if we can route a cell area downstream
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      // set to can't be evaluated
      can_eval(i,j) = false;
      // to route downstream a pixel must:
      // - have 9 finite values in the block centred on i,j
      // - not be a channel
      // - have positive area

      // set so zero upstream for later
      n_above(i,j) = 0;

      
      // work out finite neighbours
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
      // see if tests passed
      if(n_finite==9 && !is_channel(i,j) && area(i,j)>0){
	can_eval(i,j) = true;
	n_propergate = n_propergate + 1;
      }
    }
  }
  
  // loop to determin how many upstream cells each routing cell has
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      if(can_eval(i,j)==true){
	for(int ii=-1;ii<2;ii++){
	  for(int jj=-1;jj<2;jj++){
	    if( i+ii > -1 && i+ii < dem.nrow() && 
		j+jj > -1 && j+jj < dem.ncol() ){
	      if( dem(i+ii,j+jj) < dem(i,j) ){
		n_above(i+ii,j+jj) = n_above(i+ii,j+jj)+1;
	      }
	    }
	  }
	}
      }
    }
  }

  // loop to propogate area
  while( n_propergate > 0 && niter < max_iter){
    // remember cpp is 0 base
    for(int i=0;i < dem.nrow(); i++){
      for(int j=0; j < dem.ncol(); j++){
	// check it can be propogaed
	if( can_eval(i,j)==true && n_above(i,j)==0 ){

	  // work out weights
	  sW = 0; 
	  for(int ii=-1;ii<2;ii++){
	    for(int jj=-1;jj<2;jj++){
	      W(ii,jj) = max(0,(dem(i,j)-dem(i+ii,j+jj))/dist(ii+1,jj+1));
	      sW=sW+W(ii,jj);
	    }
	  }

	  // propogate
	  for(int ii=-1;ii<2;ii++){
	    for(int jj=-1;jj<2;jj++){
	      if( W(i,j) > 0 ){
		area(i+ii,j+jj) = area(i+ii,j+jj) + area(i,j)*W(ii,jj)/sW;
		n_above(ii+ii,j+jj) = n_above(i+ii,j+jj)-1;
	      }
	    }
	  }
	  // update count
	  n_propergate = n_propergate-1;
	}
	Rcpp::checkUserInterrupt();
      }
    }

    // record iteration information
    niter = niter + 1;
    Rcpp::checkUserInterrupt();
    Rcout << "The number of cells left to propogate after iteration " << niter << " is " << n_propergate << std::endl;
    
  }

  return area;
}
