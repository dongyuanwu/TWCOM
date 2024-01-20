#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
LogicalVector if_neighbors_samesample(NumericVector X, NumericVector Y, double dist) {
	
	int nspot = X.length();
	LogicalVector if_neighbor(nspot * nspot);
	int k = 0;
	
	for (int i = 0; i < nspot; ++i) {
		
		for (int j = 0; j < nspot; ++j) {
			
			if ((abs(X[i] - X[j]) <= dist) & (abs(Y[i] - Y[j]) <= dist)) {
			
					if_neighbor[k] = 1;
			
			} else if_neighbor[k] = 0;
				
			
			k += 1;
		
		}
		
	}
	
	return if_neighbor;
	
}