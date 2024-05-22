#include <RcppArmadillo.h>   
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]




// This gives the inverse wishart distribution only upto a constant

double diwishart(arma::mat X, arma::mat psi, double nu){
  
  double lkhd = pow(det(X), -(nu+X.n_rows+1)/2)*exp(-0.5*trace(psi*X.i()));
  
  return(log(lkhd));
  
}
