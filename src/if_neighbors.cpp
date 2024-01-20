#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

LogicalVector if_neighbors(NumericVector X, NumericVector Y, NumericVector sample, double dist) {
	
	int nspot = X.length();
	LogicalVector if_neighbor(nspot * nspot);
	int k = 0;
	
	for (int i = 0; i < nspot; ++i) {
		
		for (int j = 0; j < nspot; ++j) {
			
			if (sample[i] == sample[j]) {
				
				if ((abs(X[i] - X[j]) <= dist) & (abs(Y[i] - Y[j]) <= dist)) {
			
					if_neighbor[k] = 1;
			
				} else if_neighbor[k] = 0;
				
			} else if_neighbor[k] = 0;
			
			k += 1;
		
		}
		
	}
	
	return if_neighbor;
	
}