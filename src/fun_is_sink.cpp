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
LogicalMatrix fun_is_sink(NumericVector dem, int nc ){

  LogicalVector is_sink(dem.length(),dem.ncol());
  int ngh(8);
  
  double na_test_val = -10000; // test the for NAN against this - if NAN will return false
    
  // determine if a cell if it is a sink
  for(int i=0;i < dem.length(); i++){
    if( dem[i] > na_test_val ){
      // if finite see if lower cell
      is_

      // work out index of neighbours
      ngh[0] = i-nc-1;
      ngh[1] = i-nc;
      ngh[2] = i-nc+1;
      ngh[3] = i-1;
      ngh[4] = i+1;
      ngh[5] = i+nc-1;
      ngh[6] = i+nc;
      ngh[7] = i+nc+1;
    
      for(int j=0; j < 8; j++){
	if( ngh(j)< dem.length() &&
	    ngh(j)>-1 ){
	  if( dem(ngh(j)) > na_test_val ){
	    if( dem(ngh(j)) < dem(i) ){
	      is_sink = false;
	    }
	  }else{
	    is_sink=false;
	  }
	}
      }
    }else{
      is_sink=false;
    }
  }
    return(is_sink)
      }
	// then cell is in catchment
	// presume it is a sink
	is_sink(i,j) = true;
	// look at neighboring cells
	for(int ii=-1;ii<2;ii++){
	  for(int jj=-1;jj<2;jj++){
	    if( i+ii > -1 && i+ii < dem.nrow() && 
		j+jj > -1 && j+jj < dem.ncol() ){
	      // only those in range
	      if( dem(i+ii,j+jj) > na_test_val ){
		// finite value
		if( dem(i+ii,j+jj) < dem(i,j) ){
		  is_sink(i,j)=false;
		}
	      }else{
		// non finite value so boundary and not sink
		is_sink(i,j)=false;
	      }
	    }
	  }
	}
      }else{
	// na cell so can;t be a sink
	is_sink(i,j)=false;
      }
      Rcpp::checkUserInterrupt();
    }
  }
    
  return is_sink; 
}
