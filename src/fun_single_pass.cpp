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
List fun_single_pass(NumericVector dem,
		     IntegerVector channel_id,
		     NumericVector land_area,
		     IntegerVector offset,
		     NumericVector dx,
		     NumericVector cl){

  // store of the returned values
  IntegerVector order(dem.length(),NA_INTEGER);
  NumericVector upslope_area = land_area ;
  NumericVector contour_length(dem.length(),NA_REAL);
  NumericVector gradient(dem.length(),NA_REAL);
  NumericVector atanb(dem.length(),NA_REAL);
  
  // store the computational sequenence
  IntegerVector seq(dem.length(),NA_INTEGER);
  int seq_loc = 0;
    
  IntegerVector ngh(offset.length());
  IntegerVector n_higher(dem.length(),NA_INTEGER);

  int n_finite = 0;
  
  // work out number of higher cells
  for(int i=0;i < dem.length(); i++){
    if( !(NumericVector::is_na(dem(i))) &&
	(land_area(i) > 0) ){
      // has finite dem value and finite land area so needs computing
      n_higher(i) = 0;
      n_finite = 0;
      
      // neighbours
      ngh = offset + i;
      LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);

      for(int j=0;j<ngh.length();j++){
	if( in_range(j) ){
	  if( !(NumericVector::is_na(dem(ngh(j)))) ){
	    n_finite += 1;
	    if( dem(ngh(j)) > dem(i) &&
		(IntegerVector::is_na(channel_id(ngh(j))))){
	      //Rcout << "incriment higher" << "\n";
	      n_higher(i) += 1;
	    }
	  }
	}
      }
	
      // if there are no higher values then a peak flag as such
      if( n_higher(i) == 0 ){
	seq(seq_loc) = i;
	seq_loc += 1;
	order(i) = 1 ;
	if( n_finite < 8 ){
	  Rcout << "Boundary peak at cell " << i << "\n";
	}
      }
    }
  }

  Rcout << "got to populate" << "\n";
  Rcout << seq_loc << "\n";
  // loop to populate all of order
  for(int s=0;s < seq.length(); s++){
    
    if( !(IntegerVector::is_na(seq(s))) ){
      int i = seq(s);
      
      // Rcout << s << " " << i << "\n";
    
      // common parts

      // work out neighbours and if they are in range
      ngh = offset + i;
      LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);

      // initialise the weights
      NumericVector gl(8,0.0); // gradient
      contour_length(i) = 0;
      
      // see if the cell is a river cell
      if( !(IntegerVector::is_na(channel_id(i))) ){
	// this cell has partial land area and parital channel
	// use gradient draining to cell and special contour length
	// Rcout << "channel cell" << "\n";
	// loop neighbours to compute summary statisics
	for(int j=0;j<ngh.length();j++){	
	  if( in_range(j) ){
	    if( !(NumericVector::is_na(dem(ngh(j)))) &&
		(land_area(i) > 0) &&
		dem(ngh(j)) > dem(i) ){
	      gl(j) = cl(j)*( (dem(ngh(j))-dem(i)) / dx(j) ); //tan(beta)*L
	      contour_length(i) += cl(j); // Sum{L}
	    }
	  }
	}
	double sum_gl = sum(gl); 
	gradient(i) = sum_gl/contour_length(i);
	atanb(i) = log(upslope_area(i) / (cl(8)*gradient(i)));
	// all goes to river so no need to pass on

      }else{
	// this is a complete land cell
	// use downslope gradient and contour lengths
	// Rcout << "land cell" << "\n";
	
	// loop neighbours to compute summary statisics
	for(int j=0;j<ngh.length();j++){	
	  if( in_range(j) ){
	    if( !(NumericVector::is_na(dem(ngh(j)))) &&
		(land_area(i) > 0) &&
		dem(ngh(j)) < dem(i) ){
	      gl(j) = cl(j)*( (dem(i) - dem(ngh(j))) / dx(j) ); //tan(beta)*L
	      contour_length(i) += cl(j); // Sum{L}
	    }
	  }
	}
	double sum_gl = sum(gl); 
	double C = upslope_area(i)/sum_gl; //label as per Quinn 1991
	if( contour_length(i) > 0.0 ){
	  gradient(i) = sum_gl/contour_length(i);
	  atanb(i) = log(C);
	  
	  // loop neighbours to distribute downstream
	  for(int j=0;j<ngh.length();j++){	
	    if( in_range(j) ){
	      if( !(NumericVector::is_na(dem(ngh(j)))) &&
		  dem(ngh(j)) < dem(i) ){
		// spread area downstream
		upslope_area(ngh(j)) += C*gl(j);
		// update number of higher points
		n_higher(ngh(j)) -= 1;
		// add to list and assign order
		if( n_higher(ngh(j)) == 0 ){
		  order(ngh(j)) = order(i) + 1;
		  seq(seq_loc) = ngh(j);
		  seq_loc += 1;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  
  Rcout << "got to output" << "\n";
  // create output list
  List out=List::create(Named("order")=order,
			Named("upslope_area")=upslope_area,
			Named("contour_length")=contour_length,
			Named("gradient")=gradient,
			Named("atanb")=atanb);
  
  return out; 
}
