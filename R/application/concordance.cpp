#include "Rcpp.h"

// [[Rcpp::export]]
Rcpp::List cindex(Rcpp::NumericVector Y, Rcpp::NumericVector delta, 
	Rcpp::NumericVector eta) 
{
    double num = 0.0;
    double den = 0.0;

    for (int i = 0; i < Y.length(); i++) {
	for (int j=i+1; j < Y.length(); j++) {
	    double a = (Y(i) < Y(j)) && (delta(i) == 1); 
	    double b = (Y(j) < Y(i)) && (delta(j) == 1); 
	    den += a + b;
	    num += (a && (eta(i) > eta(j))) + (b && (eta(j) > eta(i)));
	}
    }

    return Rcpp::List::create(
	Rcpp::Named("cindex") = num / den,
	Rcpp::Named("concordant") = num,
	Rcpp::Named("discordant") = den - num
    );
}
