#include <RcppArmadillo.h>   
#include <Rcpp.h>

using namespace arma;
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]


arma::mat makeSymmetric(arma::mat A) {
  
  
  int row_size = A.n_rows;
  int col_size = A.n_cols;
  for (int r = 0; r < row_size; r++) {
    for (int c = r+1; c < col_size; c++) {
        A(c,r) = A(r,c);
    }
    
  }
  return A;
}
