#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector get_Ck_v(NumericVector expr, NumericMatrix spots) {
	
	int Cknrow = spots.nrow();
	
	NumericVector Ck(Cknrow);
	
	for (int i = 0; i < Cknrow; ++i) {
		
		int spoti = spots(i, 0);
		int spotj = spots(i, 1);
		
		double lexpr = expr[spoti-1];
		double rexpr = expr[spotj-1];
		
		Ck[i] = lexpr * rexpr;
		
	}
	
	return Ck;
	
}