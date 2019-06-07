// Sink filling algorithm applied using Rcpp //

#include <Rcpp.h>
using namespace Rcpp;

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
  LogicalMatrix is_valid_from(dem.nrow(),dem.ncol());
  LogicalMatrix is_valid_to(dem.nrow(),dem.ncol());

  Rcout << "Start of function" << std::endl;
  
  // loop to determin cells with valid dem, area and hillslope class values
  // these are the cells we can drain from 
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      // evaluate if we can drain from this cell
      // valid dem, land_area and hillslope class
      is_valid_from(i,j)=false;
      if(dem(i,j) > -10000){
	if(land_area(i,j) > 0){
	  if(hillslope(i,j) > 0){
	    is_valid_from(i,j) = true;
	  }
	}
      }
      // evaluate if we can drain to this cell
      // valid dem, and either hillslope class & area > 0 or river class
      is_valid_to(i,j)=false;
      if(dem(i,j) > -10000){
	if(land_area(i,j) > 0){
	  if(hillslope(i,j) > 0){
	    is_valid_to(i,j) = true;
	  }
	}
	if(channel(i,j) > 0){
	  is_valid_to(i,j) = true;
	}
      }
    }
  }
  Rcout << "Computed valid sources and downsteams" << std::endl;
  
  // Loop to compute the fluxes
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      if(is_valid_from(i,j)==true){
	// Rcout << "Start to determine initial index" << std::endl;
	int upslope_index = hillslope(i,j);
	int downslope_index = -1;
	// Rcout << "Determined initial indexes" << std::endl;
	if( channel(i,j) > -10000 ){ // then it is a channel pixel
	  // all land area goes to channel
	  downslope_index = channel(i,j) ;
	  //Rcout << downslope_index << " " << upslope_index << std::endl;
	  channel_redist(downslope_index,upslope_index) =
	    channel_redist(downslope_index,upslope_index) + land_area(i,j) ;
	  //Rcout << "Evaluated as downsteam channel" << std::endl;
	}else{
	  // sent to downslope pixels proportional to gradient
	  // compute total gradient
	  double sum_tanb = 0;
	  for(int ii=-1;ii<2;ii++){
	    for(int jj=-1;jj<2;jj++){
	      if( i+ii > -1 && i+ii < dem.nrow() && 
		  j+jj > -1 && j+jj < dem.ncol() &&
		  ii != 0 && jj != 0 &&
		  is_valid_to(i+ii,j+jj) ==true){
		if( dem(i,j) > dem(i+ii,j+jj) ){
		  sum_tanb = sum_tanb + ( (dem(i,j) - dem(i+ii,j+jj))/dist(ii+1,jj+1) );
		}
	      }
	    }
	  }
	  //Rcout << "Computed sum_tanb" << std::endl;
	  // redistribute
	  for(int ii=-1;ii<2;ii++){
	    for(int jj=-1;jj<2;jj++){
	      if( i+ii > -1 && i+ii < dem.nrow() && 
		  j+jj > -1 && j+jj < dem.ncol() &&
		  ii != 0 && jj != 0 &&
		  is_valid_to(i+ii,j+jj) ==true){
		if( dem(i,j) > dem(i+ii,j+jj) ){
		  double w = ( (dem(i,j) - dem(i+ii,j+jj))/dist(ii+1,jj+1) ) / sum_tanb;
		  if( is_valid_from(i+ii,j+jj)==true ){
		    // send to the hillslope class
		    int downslope_index = hillslope(i+ii,j+jj) ;
		    //Rcout << downslope_index << " " << upslope_index << std::endl;
		    slope_redist( downslope_index,upslope_index) =
			slope_redist(downslope_index,upslope_index) + w*land_area(i,j) ;
		    //Rcout << "Evaluated as hillslope class" << std::endl;
		  }else{
		    // send to the channel class if there is one
		    if( channel(i+ii,j+jj) > 0 ){
		      downslope_index = channel(i+ii,j+jj);
		      channel_redist(downslope_index,upslope_index) =
			channel_redist(downslope_index,upslope_index) + w*land_area(i,j) ;
		      //Rcout << "Evaluated as channel" << std::endl;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //Rcout << "Got to output" << std::endl;
  //return slope_redist;
  return List::create(Named("hillslope") = slope_redist , Named("channel") = channel_redist);
}
