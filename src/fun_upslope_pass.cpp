#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for passing up the catchments from river nodes
//' 
//' @param dem Digital elevation model as a vector
//' @param order initial order values as a vector, internally starts at one and move upstream
//' @param nc number of columns in the matrix
//' 
//' @return a list with the filled dem and order
//'
// [[Rcpp::export]]
List fun_upslope_pass(NumericVector dem,
		      IntegerVector channel_id,
		      IntegerVector offset){

  // store of the order
  IntegerVector order(dem.length(),NA_INTEGER);
  // store the computational sequenence
  IntegerVector seq(dem.length(),NA_INTEGER);
  int seq_loc = 0;
    
  IntegerVector ngh(offset.length());
  IntegerVector n_lower(dem.length(),NA_INTEGER);

  int n_finite = 0;
  
  // work out number of lower cells
  for(int i=0;i < dem.length(); i++){
    if( !(NumericVector::is_na(dem(i))) ){
      // has finite dem value
      
      // check if cell is in a channel channel
      if( !(IntegerVector::is_na(channel_id(i))) ){
	// then a channel - set order to 1 ans n_lower so not evaluated again
	order(i) = 1;
	n_lower(i) = -99;
      }else{
	// then not a channel so wokr out how many lower cells
	n_lower(i) = 0;
	n_finite = 0;
      
	// neighbours
	ngh = offset + i;
	LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);

	for(int j=0;j<ngh.length();j++){
	  if( in_range(j) ){
	    if( !(NumericVector::is_na(dem(ngh(j)))) ){
	      n_finite += n_finite;
	      if( dem(ngh(j)) < dem(i) ){
		n_lower(i) = n_lower(i) + 1;
	      }
	    }
	  }
	}
	
	// if there are no lower values then either sinkor  edge drain
	if( n_lower(i) == 0 ){
	  if( n_finite == 8 ){
	    //sink....
	    Rcout << "Sink at pixel" << i << "\n";
	  }else{
	    Rcout << "Edge drain at pixel" << i << "\n";
	    order(i) = 1;
	  }
	}
      }
    }
  }

  Rcout << "got to seq_loc" << "\n";
  
  // loop order to get sequence vector
  seq_loc = 0;
  for(int i=0;i < dem.length(); i++){
    if( order(i) == 1 ){
      seq(seq_loc) = i;
      seq_loc += 1;
    }
  }

  Rcout << "got to populate" << "\n";
  
  // loop to populate all of order
  for(int i=0;i < seq.length(); i++){
    //Rcout << i << "\n";
    if( !(IntegerVector::is_na(seq(i))) ){
      // work out neighbours and if they are in range
      ngh = offset + seq(i);
      LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);
      // loop neighbours
      for(int j=0;j<ngh.length();j++){	
	if( in_range(j) ){
	  if( !(NumericVector::is_na(dem(ngh(j)))) &&
	      dem(ngh(j)) > dem(seq(i)) ){
	    n_lower(ngh(j)) = n_lower(ngh(j)) - 1;
	    if( n_lower(ngh(j))==0 ){
	      order(ngh(j)) = order(seq(i))+1;
	      seq(seq_loc) = ngh(j);
	      seq_loc += 1;
	    }
	  }
	}
      }
    }
  }
  Rcout << "got to output" << "\n";
  // create output list
  List out=List::create(Named("dem_filled")=dem,
			Named("order")=order,
			Named("seq")=seq);
  
  return out; 
}
