#include <RcppArmadillo.h>   
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]


arma::mat myround(arma::mat A, int digits = 0) {
  arma::mat B = A;
  
  int row_size = A.n_rows;
  int col_size = A.n_cols;
  for (int r = 0; r < row_size; r++) {
    for (int c = 0; c < col_size; c++) {
      if(r!=c){
        B(r,c) = ::Rf_fround(A(r,c), digits);
        }
    }
    
  }
  return B;
}
