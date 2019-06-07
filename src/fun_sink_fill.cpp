// Sink filling algorithm applied using Rcpp //

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fun_sink_fill(NumericMatrix dem, LogicalMatrix is_channel,
			    NumericMatrix delta, int max_iter){

  int n_sink = 10;
  int niter = 0;
  LogicalMatrix can_fill(dem.nrow(),dem.ncol());
  
  // determine if we can fill it if it is a sink
  for(int i=0;i < dem.nrow(); i++){
      for(int j=0; j < dem.ncol(); j++){
	can_fill(i,j) = false;
	if( dem(i,j) > -1000000 ){ // will be FALSE for any NAN dem values
	  can_fill(i,j)=true;
	  for(int ii=-1;ii<2;ii++){
	      for(int jj=-1;jj<2;jj++){
		if( i+ii > -1 && i+ii < dem.nrow() && 
		    j+jj > -1 && j+jj < dem.ncol() &&
		    ii != 0 && jj != 0){
		  if( dem(i+ii,j+jj) > -10000 ){
		    can_fill(i,j)=true;
		  }else{
		    can_fill(i,j)=false; // can't fill its on the edge
		  }
		}
	      }
	  }
	}
      }
  }

  // loop to fill
  while( n_sink > 0 && niter < max_iter){
    n_sink = 0;
    // rember cpp is 0 base
    for(int i=0;i < dem.nrow(); i++){
      for(int j=0; j < dem.ncol(); j++){
	//Rcout << "dem(i,j) " << dem(i,j) << std::endl;
	// check if finite
	if( can_fill(i,j)==true ){ // we can fill it if it is a sink
	  //check if a channel
	  if(is_channel(i,j)==false){
	    // presume it is a sink
	    bool is_sink = true;
	    double lowest_neighbour = 1e32;
	    double new_value = 1e32;
	    
	    // check its status
	    for(int ii=-1;ii<2;ii++){
	      for(int jj=-1;jj<2;jj++){
		if( i+ii > -1 && i+ii < dem.nrow() && 
		    j+jj > -1 && j+jj < dem.ncol() &&
		    ii != 0 && jj != 0){
		  double neighbour = dem(i+ii,j+jj);
		  //Rcout << "sink " << is_sink << std::endl;
		  //Rcout << "neighbour " << neighbour << std::endl;
		  if( neighbour < lowest_neighbour ){
		    lowest_neighbour = neighbour;
		    new_value = neighbour + delta(ii+1,jj+1);
		  }
		  if( neighbour < dem(i,j) ){
		    is_sink=false;
		  }
		}
	      }
	    }
	
	    // if it is still a sink
	    if(is_sink==true){
	      dem(i,j) = new_value;
	      n_sink = n_sink + 1 ;
	    }
	  }
	}
	Rcpp::checkUserInterrupt();
      }
    }
    niter = niter + 1;
    Rcpp::checkUserInterrupt();
    Rcout << "The number of sinks after iteration is " << niter << " is " << n_sink << std::endl;

  }
  
  return dem; 
}
