#include <RcppArmadillo.h>   
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]


// When cor = 0 we have that the kalman filter of the joint is giving the same values
// as when running two independent kalman filters



double dmvnormC(arma::vec X, arma::vec mean, arma::mat cov){
  double lkhd = 0;
  if(X.has_nan() == TRUE){
    return(0);
  }else{
    arma::vec res = X-mean;
    lkhd = as_scalar(res.t()*cov.i()*res);
  }
  
  
  return(-0.5*lkhd-0.5*log(det(2*datum::pi*cov)));
  
}
