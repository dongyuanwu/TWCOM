#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_Ck_strict(NumericMatrix expr, NumericMatrix spots, NumericVector ligand, NumericVector receptor) {
	
	int Cknrow = spots.nrow();
	
	NumericVector Ck(Cknrow);
	
	for (int i = 0; i < Cknrow; ++i) {
		
		int spoti = spots(i, 0);
		int spotj = spots(i, 1);
		
		double lexpr = 1;
		for (int a = 0; a < ligand.length(); ++a) {
			lexpr = lexpr * expr(ligand[a]-1, spoti-1);
		}
		lexpr = pow(lexpr, (1.0/ligand.length()));
		
		double rexpr = 1;
		for (int a = 0; a < receptor.length(); ++a) {
			rexpr = rexpr * expr(receptor[a]-1, spotj-1);
		}
		rexpr = pow(rexpr, (1.0/receptor.length()));
		
		Ck[i] = lexpr * rexpr;
		
	}
	
	return Ck;
	
}


// [[Rcpp::export]]
NumericVector get_Ck(NumericMatrix expr, NumericMatrix spots, NumericVector ligand, NumericVector receptor) {
	
	int Cknrow = spots.nrow();
	
	NumericVector Ck(Cknrow);
	
	for (int i = 0; i < Cknrow; ++i) {
		
		int spoti = spots(i, 0);
		int spotj = spots(i, 1);
		
		double lexpr = 0;
		for (int a = 0; a < ligand.length(); ++a) {
			lexpr = lexpr + expr(ligand[a]-1, spoti-1);
		}
		lexpr = lexpr / ligand.length();
		
		double rexpr = 0;
		for (int a = 0; a < receptor.length(); ++a) {
			rexpr = rexpr + expr(receptor[a]-1, spotj-1);
		}
		rexpr = rexpr / receptor.length();
		
		Ck[i] = lexpr * rexpr;
		
	}
	
	return Ck;
	
}
