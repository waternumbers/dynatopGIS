#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List rcpp_fN(int k, NumericVector rst_prop){
  //Rcout << "in fN ffss" << "\n";
  
  int nc = rst_prop(1);
  int nr = rst_prop(0);
  double dx = rst_prop(2);
  double dxd = dx*sqrt(2);
  double cl = dx*0.5;
  double cld = dx*0.35;
  //double cl = dx /(1+sqrt(2));
  
  int row = (k/nc); // row - should be integer division
  //Rcout << row << "\n";
  int col = k - row*nc; //column
  //Rcout << col << "\n";
  
  IntegerVector new_col = IntegerVector::create(col-1,col,col+1,col-1,col+1,col-1,col,col+1);
  IntegerVector new_row = IntegerVector::create(row-1,row-1,row-1,row,row,row+1,row+1,row+1);
  IntegerVector new_idx = new_row*nc + new_col;
  NumericVector new_dx = NumericVector::create(dxd,dx,dxd,dx,dx,dxd,dx,dxd);
  NumericVector new_cl = NumericVector::create(cld,cl,cld,cl,cl,cld,cl,cld);

  LogicalVector is_valid = (new_col>-1) & (new_col < nc) & (new_row>-1) & (new_row<nr);

  //Rcout << new_col.length() << "\n";
  //Rcout << new_row.length() << "\n";
  //Rcout << new_idx.length() << "\n";
  //Rcout << new_dx.length() << "\n";
  //Rcout << new_cl.length() << "\n";
  //Rcout << is_valid.length() << "\n";
  
  return List::create(Named("idx") = new_idx[is_valid],
		      Named("dx") = new_dx[is_valid],
		      Named("cl") = new_cl[is_valid]);
}


//' cpp wrapper function for computation of hru properties
//' 
//' @param dem Digital elevation model
//' @param hillslope_id hillslope class of each pixel
//' @param channel_id UID of channel present in the pixel
//' @param land_area Area of land surface in pixel
//' @param channel_area Area of channel surface in pixel
//' @param grad gradient of the cell
//' @param atb topographic index ln(upslope area / tan beta) of the cell
//' @param W saturated zone redistribution matrix
//' @param area HSU area
//' @param s_bar HSU average gradient
//' @param atb_bar HSU averable topographic index
//' @param rst_prop dimension ans resolution of raster converted to vectors
// [[Rcpp::export]]
arma::sp_mat rcpp_hru(NumericVector dem,
		      IntegerVector hillslope_id,
		      IntegerVector channel_id,
		      NumericVector land_area,
		      NumericVector channel_area,
		      NumericVector grad,
		      NumericVector atb,
		      arma::sp_mat W,
		      NumericVector area,
		      NumericVector s_bar,
		      NumericVector atb_bar,
		      NumericVector rst_prop){
  
  
  for( int i=0; i<dem.length(); i++){
    
    
    //Rcpp::Rcout << dem(i) << "\n";
    if( !(NumericVector::is_na(dem(i))) ){
      //Rcout << i << "\n";
      // then a hillslope or channel cell
      bool is_channel = !(IntegerVector::is_na( channel_id(i) ));
      bool has_land = land_area(i) > 0 ;
      
      //Rcpp::Rcout << is_channel << "\n";
      //Rcpp::Rcout << channel_id(i) << "\n";
      //Rcpp::Rcout << has_land << "\n";
      //Rcpp::Rcout << land_area(i) << "\n";
      if( is_channel ){
	//Rcout << "in channel" << "\n";
	area( channel_id(i)-1 ) += channel_area(i);
	if( has_land ){
	  area( hillslope_id(i)-1 ) += land_area(i);
	  s_bar( hillslope_id(i)-1 ) += grad(i)*land_area(i);
	  atb_bar( hillslope_id(i)-1 ) += atb(i)*land_area(i);
	  //Rcout << "assigning to W" << "\n";
	  W(channel_id(i)-1,hillslope_id(i)-1) += land_area(i);
	  
	}
      }else{
	//Rcout << "in the not channel part" << "\n";
	
	area( hillslope_id(i)-1 ) += land_area(i);
	s_bar( hillslope_id(i)-1 ) += grad(i)*land_area(i);
	atb_bar( hillslope_id(i)-1 ) += atb(i)*land_area(i);
	
	List j = rcpp_fN(i,rst_prop);
	
	//Rcout << "Got neighbours" << "\n";
	//if( i == 1946138 ){
	//  IntegerVector idx = j["idx"];
	//  Rcout << idx << "\n";
	//}

	NumericVector cl = j["cl"];
	NumericVector dx = j["dx"];
	IntegerVector idx = j["idx"];
	NumericVector ddn = dem[idx];
	NumericVector ddi (idx.length(),dem(i));

	//Rcout << idx << "\n";
	//Rcout << dx << "\n";
	//Rcout << cl << "\n";

	//Rcout << ddn << "\n";
	//Rcout << ddi << "\n";
	
	NumericVector grd = (ddi - ddn) / dx ;
	//Rcout << grd << "\n";
		
	LogicalVector is_lower = (grd>0) ;
	is_lower[is_na(is_lower)] = 0; // set na values to false
	
	//LogicalVector is_shit = is_na(grd);
	//Rcout << is_shit << "\n";				      
	//Rcout << is_lower << "\n";
	if( any(is_lower).is_true() ){

	  NumericVector trim_cl = cl[is_lower];
	  IntegerVector trim_idx = idx[is_lower];
	  NumericVector gcl = s_bar(hillslope_id(i)-1)*trim_cl;
	  double sum_gcl = Rcpp::sum( gcl );
	  NumericVector frc = gcl / sum_gcl;
	  //Rcout << "frc: " << frc << "\n";
	  //Rcout << "is_lower: " << is_lower << "\n";

	  for(int k=0; k< trim_idx.length(); k++){
	    if( land_area(trim_idx(k)) > 0 ){
	      W(hillslope_id(trim_idx(k))-1,hillslope_id(i)-1) += frc(k)*land_area(i);
	    }else{
	      W(channel_id(trim_idx(k))-1,hillslope_id(i)-1) += frc(k)*land_area(i);
	    }
	  }
	  
	}
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return W;
}
