#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix get_Ms(NumericMatrix &M, NumericMatrix &spots) {
	
	int ngrp = M.ncol();
	int npair = spots.nrow();
	
	NumericMatrix Ms(npair, ngrp*ngrp);
	
	for (int k = 0; k < npair; ++k) {
		
		int spoti = spots(k, 0);
		int spotj = spots(k, 1);
		
		NumericVector Mij(ngrp*ngrp);
		int posi = 0;
		
		for (int i = 0; i < ngrp; ++i) {
			
			for (int j = 0; j < ngrp; ++j) {
				
				Mij[posi] = M(spoti - 1, i) * M(spotj - 1, j);
				posi += 1;
				
			}
			
		}
		
		Ms(k, _) = Mij;
		
	}
	
	return Ms;
	
}