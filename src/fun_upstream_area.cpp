// Sink filling algorithm applied using Rcpp //

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fun_av_tanb(NumericMatrix dem, NumericMatrix area, LogicalMatrix is_channel,
		      NumericMatrix dist, int max_iter){

  NumericMatrix atb(dem.nrow(),dem.ncol());
  NumericMatrix av_tanb(dem.nrow(),dem.ncol());
  LogicalMatrix is_valid(dem.nrow(),dem.ncol());
  //IntegerMatrix n_above_not_eval(dem.nrow(),dem.ncol());
  // int n_to_eval = 0;

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

  // loop to compute average tan beta
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      if(is_valid(i,j)==true){
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
		  //Rcout << i << " " <<j<< " " <<ii<< " " <<jj<< std::endl;
		  //Rcout << dem(i,j) << " " << dem(i+ii,j+jj) << " " << dist(ii+1,jj+1)  << std::endl;
		  //Rcout << ( (dem(i,j) - dem(i+ii,j+jj))/dist(ii+1,jj+1) ) << std::endl;
		  sum_tanb = sum_tanb + ( (dem(i,j) - dem(i+ii,j+jj))/dist(ii+1,jj+1) );
		  n_downslope = n_downslope + 1;
		}
	      }
	    }
	  }
	  if( n_downslope > 0 ){
	    av_tanb(i,j) = sum_tanb / n_downslope;
	  }
	}else{
	  // then evaluate the cell using teh upslope gradients
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
	  }
	}
      }
    }
  }

  return av_tanb;
}

// 	  // evaluate using the upslope gradients
	  
//       atb(i,j) = R_NaN;
//       is_valid(i,j)=false;
//       if(dem(i,j) > -10000){
// 	if(area(i,j) > 0){
// 	  is_valid(i,j) = true;
// 	}
//       }
//     }
//   }

  
//   // initialise atb and the number of upstream cells n_not_evaluated
//   n_to_eval = 0;
//   for(int i=0;i < dem.nrow(); i++){
//     for(int j=0; j < dem.ncol(); j++){
//       // set to NaN to start with
//       n_above_not_eval(i,j) = R_NaN;
//       atb(i,j) = R_NaN;
      
//       // check if finite
//       if( dem(i,j) > -1000000 ){// will be FALSE for any NAN dem values
// 	if( area(i,j) > 0 ){ // no point evaluating any points without rivers
// 	  // counter number of bordering cells whch are higher and not rivers
// 	  int n_higher = 0;
// 	  for(int ii=-1;ii<2;ii++){
// 	    for(int jj=-1;jj<2;jj++){
// 	      if( i+ii > -1 && i+ii < dem.nrow() && 
// 		  j+jj > -1 && j+jj < dem.ncol() &&
// 		  ii != 0 && jj != 0){
// 		if( dem(i+ii,j+jj) > -10000 ){
// 		  if( area(i+ii,j+jj) > 0 ){
// 		    // then valid
// 		    if( dem(i+ii,j+jj) > dem(i,j) && is_channel(i+ii,j+jj) == false ){
// 		      // then higher
// 		      n_higher = n_higher + 1;
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	n_above_not_eval(i,j) = n_higher;
// 	n_to_eval = n_to_eval + 1;
//       }
//     }
//   }
//   Rcout << "The number of cells to evaluate is " << n_to_eval << std::endl;

//   // Now loop and process those with no upstream cells left to process
//   int niter=0;
//   while( (n_to_eval > 0) & (niter < max_iter) ){
    
//     Rcpp::checkUserInterrupt();

//     for(int i=0;i < dem.nrow(); i++){
//       for(int j=0; j < dem.ncol(); j++){
//   	if( n_above_not_eval(i,j) == 0 ){
// 	  double n_downslope = 0;
// 	  double sum_tanb = 0;
// 	  if( is_channel(i,j) == false ){
// 	    // then evaluate the cell using downstream gradients
// 	    for(int ii=-1;ii<2;ii++){
// 	      for(int jj=-1;jj<2;jj++){
// 		if( i+ii > -1 && i+ii < dem.nrow() && 
// 		    j+jj > -1 && j+jj < dem.ncol() &&
// 		    ii != 0 && jj != 0){
// 		  if( dem(i+ii,j+jj) > -10000 ){ // will be false if NAN
// 		    if( dem(i,j) > dem(i+ii,j+jj) && is_channel(i+ii,j+jj) == false){ // will be false if higher or NAN
// 		      sum_tanb = sum_tanb + ( (dem(i,j) - dem(i+ii,j+jj))/dist(ii,jj) );
// 		      n_downslope = n_downslope + 1;
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	    if( n_downslope > 0 ){
// 	      // propogate the area and evaluation
// 	      for(int ii=-1;ii<2;ii++){
// 		for(int jj=-1;jj<2;jj++){
// 		  if( i+ii > -1 && i+ii < dem.nrow() && 
// 		      j+jj > -1 && j+jj < dem.ncol() &&
// 		      ii != 0 && jj != 0){
// 		    if( dem(i+ii,j+jj) > -10000 ){ // will be false if NAN
// 		      if( dem(i,j) > dem(i+ii,j+jj) ){ // will be false if higher or NAN
// 			double w = ( (dem(i,j) - dem(i+ii,j+jj))/dist(ii,jj) ) / sum_tanb;
// 			area(i+ii,j+jj) = area(i+ii,j+jj) + w*area(i,j);
// 			n_above_not_eval(i+ii,j+jj) = n_above_not_eval(i+ii,j+jj)-1;
// 		      }
// 		    }
// 		  }
// 		}
// 	      }
// 	      // compute atb
// 	      Rcout << area(i,j) << " " << sum_tanb << " " << n_downslope << std::endl;
// 	      atb(i,j) = log(area(i,j)) / (sum_tanb/n_downslope);
// 	    }
// 	  }else{
// 	    // channel cell use upstream gradients
// 	    for(int ii=-1;ii<2;ii++){
// 	      for(int jj=-1;jj<2;jj++){
// 		if( i+ii > -1 && i+ii < dem.nrow() && 
// 		    j+jj > -1 && j+jj < dem.ncol() &&
// 		    ii != 0 && jj != 0){
// 		  if( dem(i+ii,j+jj) > -10000 ){ // will be false if NAN
// 		    if( dem(i,j) < dem(i+ii,j+jj) ){ // will be false if higher or NAN
// 		      sum_tanb = sum_tanb + ( (dem(i+ii,j+jj) - dem(i,j))/dist(ii,jj) );
// 		      n_downslope = n_downslope + 1;
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	    if( n_downslope > 0 ){
// 	      atb(i,j) = log(area(i,j)) / (sum_tanb/n_downslope);
// 	    }
// 	  }
// 	  // alter counts and flags
// 	  n_to_eval = n_to_eval - 1;
// 	  n_above_not_eval(i,j)=-999;
// 	}
//       }
//     }
    
//     niter = niter + 1;
//     Rcout << "The number of cells left to evaluate on iteration " << niter << " is " << n_to_eval << std::endl;
//   }

//   return atb;
// }
