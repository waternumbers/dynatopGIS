#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for computation of redistribution matrices
//' 
//' @param dem Digital elevation model
//' @param land_area Area of land surface in pixel
//' @param hillslope hillslope class of each pixel
//' @param channel channel class of each pixel
//' @param number_hillslope_class number of hill slope classes
//' @param number_channel_class number of channel classes
//' @param dist 3x3 matrix of distances to ajoining cells
//'
//' @return list of hillslope and channel properties
//'
// [[Rcpp::export]]
List fun_redistribution(NumericMatrix dem, NumericMatrix land_area,
			IntegerMatrix hillslope, IntegerMatrix channel,
			int number_hillslope_class,
			int number_channel_class,
			NumericMatrix dist){
  
  // initialise the output
  NumericMatrix slope_redist(number_hillslope_class,number_hillslope_class);
  NumericMatrix channel_redist(number_channel_class,number_hillslope_class);
  // logical test matrix
  LogicalMatrix can_eval(dem.nrow(),dem.ncol());
  int uc, dc, n_finite;
  NumericMatrix W(3,3);
  double sW;
  
  double na_test_val = -10000; // test for NAN against this - if NAN will return false
  
  //Rcout << "Start of function" << std::endl;
  
  // loop to determin cells can be drained from
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      // to be drained froma pixle must have:
      // - a  hillslope class
      // - have a positive land area
      // - either:
      //     - a channel class
      //     - or 9 finite values in the block around i,j
            
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
      if( (n_finite==9 || channel(i,j)>na_test_val) &&
	  land_area(i,j)>0 &&
	  hillslope(i,j) > na_test_val ){
	can_eval(i,j) = true;
      }
    }
  }

  //Rcout << "Computed valid sources and downsteams" << std::endl;
  
  // the see how cells redistribute
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      // check it can be propogaed
      if( can_eval(i,j)==true ){
	
	if( channel(i,j) > na_test_val ){
	  // drains straight to the channel
	  uc = hillslope(i,j);
	  dc = channel(i,j);
	  channel_redist( dc,uc ) = channel_redist( dc,uc ) + land_area(i,j) ;
	}else{
	  // drains to various areas

	  // work out weights in each direction
	  sW = 0; 
	  for(int ii=-1;ii<2;ii++){
	    for(int jj=-1;jj<2;jj++){
	      W(ii+1,jj+1) = (dem(i,j)-dem(i+ii,j+jj))/dist(ii+1,jj+1);
	      if( W(ii+1,jj+1) > 0 ){
		sW=sW+W(ii+1,jj+1);
	      }
	    }
	  }

	  // propogate
	  for(int ii=-1;ii<2;ii++){
	    for(int jj=-1;jj<2;jj++){
	      if( W(ii+1,jj+1) > 0 ){
		if( hillslope(i+ii,j+jj) > na_test_val ){
		  uc = hillslope(i,j);
		  dc = hillslope(i+ii,j+jj);
		  // pass to hillslope class
		  slope_redist(dc,uc) = slope_redist(dc,uc) + 
		    land_area(i,j)*W(ii+1,jj+1)/sW;
		}else{
		  if( channel(i+ii,j+jj) > na_test_val ){
		    uc = hillslope(i,j);
		    dc = channel(i+ii,j+jj);
		    channel_redist(dc,uc) = channel_redist(dc,uc) + 
		      land_area(i,j)*W(ii+1,jj+1)/sW;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return List::create(Named("hillslope") = slope_redist , Named("channel") = channel_redist);
}
