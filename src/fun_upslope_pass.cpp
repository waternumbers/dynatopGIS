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
		      IntegerVector order,
		      IntegerVector offset){
  
  // store the computational sequenence
  IntegerVector seq(dem.length(),NA_INTEGER);
  int seq_loc = 0;
    
  IntegerVector ngh(offset.length());
  IntegerVector n_lower(dem.length(),NA_INTEGER);

  int cnt,dp;
  
  // work out number of lower cells
  for(int i=0;i < dem.length(); i++){
    if( !(NumericVector::is_na(dem(i))) ){
      n_lower(i) = 0;
      
      // neighbours
      ngh = offset + i;
      LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);

      for(int j=0;j<ngh.length();j++){
	if( in_range(j) ){
	  if( !(NumericVector::is_na(dem(ngh(j)))) &&
	      dem(ngh(j)) < dem(i) ){
	    n_lower(i) = n_lower(i) + 1;
	  }
	}
      }
    }
  }

  // loop order to get sequence vector
  seq_loc = 0;
  for(int i=0;i < dem.length(); i++){
    if( order(i) == 1 ){
      seq(seq_loc) = i;
      seq_loc += 1;
    }
  }
  
  // loop to populate all of order
  for(int i=0;i < seq.length(); i++){
    if( !(NumericVector::is_na(seq(i))) ){
      // work out nieghbours an if they are in range
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
  
  // create output list
  List out=List::create(Named("dem_filled")=dem,
			Named("order")=order,
			Named("seq")=seq);
  
  return out; 
}
