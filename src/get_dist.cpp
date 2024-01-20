#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericVector get_dist(NumericVector x, NumericVector y, NumericMatrix spots) {
	
	int nrow = spots.nrow();
	NumericVector out(nrow);
	
	for (int i = 0; i < nrow; ++i) {
		
		int spoti = spots(i, 0);
		int spotj = spots(i, 1);
		
		double dist = sqrt(pow(x[spoti - 1] - x[spotj - 1], 2.0) + pow(y[spoti - 1] - y[spotj - 1], 2.0));
		out[i] = dist;
	}
	
	return out;
	
}