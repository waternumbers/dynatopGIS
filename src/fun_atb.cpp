#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for computation of atb topographic wetness index
//' 
//' @param dem Digital elevation model
//' @param area Area of land surface in pixel
//' @param is_channel TRUE is a channel pixel
//' @param dist 3x3 matrix of distances to each adjacent pixel
//' @param max_iter maximum number of iterations
//'
//' @return list of hillslope and channel properties
//'
//' @details Computes topographic index area divided by tan(beta)
// [[Rcpp::export]]
NumericMatrix fun_atb(NumericMatrix dem, NumericMatrix area, LogicalMatrix is_channel,
		      NumericMatrix dist, int max_iter){

  NumericMatrix atb(dem.nrow(),dem.ncol());
  NumericMatrix av_tanb(dem.nrow(),dem.ncol());
  LogicalMatrix is_valid(dem.nrow(),dem.ncol());
  IntegerMatrix n_above(dem.nrow(),dem.ncol());
  int n_to_eval = 0;

  // loop to determin cells with valid dem and area values
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      atb(i,j) = R_NaN;
      av_tanb(i,j) = R_NaN;
      is_valid(i,j)=false;
      if(dem(i,j) > -10000){
	if(area(i,j) > 0){
	  is_valid(i,j) = true;
	}
      }
    }
  }

  // find the number of cells upslope that are valid and drain to the cell
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      if(is_valid(i,j)==false){
	n_above(i,j) = R_NaN; // dummy value
      }
      if(is_valid(i,j)==true){
	n_above(i,j) = 0;
	n_to_eval = n_to_eval + 1 ;
	// then look for all ustream cells
	for(int ii=-1;ii<2;ii++){
	  for(int jj=-1;jj<2;jj++){
	    if( i+ii > -1 && i+ii < dem.nrow() && 
		j+jj > -1 && j+jj < dem.ncol() &&
		ii != 0 && jj != 0 &&
		is_valid(i+ii,j+jj) ==true){
	      if( dem(i,j) < dem(i+ii,j+jj) ){
		n_above(i,j) = n_above(i,j) + 1;
	      }
	    }
	  }
	}
      }
    }
  }
  Rcout << "The number of cells to evaluate is " << n_to_eval << std::endl;
  
  // keep iterating data until all is done
  int niter=0;
  while( (n_to_eval > 0) & (niter < max_iter) ){
    
    // loop to cells
    for(int i=0;i < dem.nrow(); i++){
      for(int j=0; j < dem.ncol(); j++){
	if(is_valid(i,j)==true && n_above(i,j)==0){
	  // compute av_tanb, area and atb
	  double n_downslope = 0.0;
	  double sum_tanb = 0.0;
	  
	  if( is_channel(i,j) == false ){
	    // then evaluate the cell using downslope gradients
	    for(int ii=-1;ii<2;ii++){
	      for(int jj=-1;jj<2;jj++){
		if( i+ii > -1 && i+ii < dem.nrow() && 
		    j+jj > -1 && j+jj < dem.ncol() &&
		    ii != 0 && jj != 0 &&
		    is_valid(i+ii,j+jj) ==true){
		  if( dem(i,j) > dem(i+ii,j+jj) ){
		    sum_tanb = sum_tanb + ( (dem(i,j) - dem(i+ii,j+jj))/dist(ii+1,jj+1) );
		    n_downslope = n_downslope + 1;
		  }
		}
	      }
	    }
	    if( n_downslope > 0 ){
	      // compute av_tanb and atb
	      av_tanb(i,j) = sum_tanb / n_downslope;
	      atb(i,j) = log( area(i,j) / av_tanb(i,j) );
	      // assign area downstream
	      for(int ii=-1;ii<2;ii++){
		for(int jj=-1;jj<2;jj++){
		  if( i+ii > -1 && i+ii < dem.nrow() && 
		      j+jj > -1 && j+jj < dem.ncol() &&
		      ii != 0 && jj != 0 &&
		      is_valid(i+ii,j+jj) ==true){
		    if( dem(i,j) > dem(i+ii,j+jj) ){
		      double w = ( (dem(i,j) - dem(i+ii,j+jj))/dist(ii+1,jj+1) ) / sum_tanb;
		      area(i+ii,j+jj) = area(i+ii,j+jj) + w*area(i,j);
		      n_above(i+ii,j+jj) = n_above(i+ii,j+jj) - 1;
		    }
		  }
		}
	      }
	    }
	  }else{
	    // then evaluate the cell using the upslope gradients
	    for(int ii=-1;ii<2;ii++){
	      for(int jj=-1;jj<2;jj++){
		if( i+ii > -1 && i+ii < dem.nrow() && 
		    j+jj > -1 && j+jj < dem.ncol() &&
		    ii != 0 && jj != 0 &&
		    is_valid(i+ii,j+jj) ==true){
		  if( dem(i,j) < dem(i+ii,j+jj) ){
		    sum_tanb = sum_tanb + ( (dem(i+ii,j+jj) - dem(i,j))/dist(ii+1,jj+1) );
		    n_downslope = n_downslope + 1;
		  }
		}
	      }
	    }
	    if( n_downslope > 0 ){
	      av_tanb(i,j) = sum_tanb / n_downslope;
	      atb(i,j) = log( area(i,j) / av_tanb(i,j) );
	    }
	  }
	  n_to_eval = n_to_eval - 1 ;
	  n_above(i,j) = -999 ;
	}
      }
    }

    niter = niter + 1;
    Rcout << "The number of cells left to evaluate on iteration " << niter << " is " << n_to_eval << std::endl;
  }

  return atb;
}
