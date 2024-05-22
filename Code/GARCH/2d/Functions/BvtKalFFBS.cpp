#include <RcppArmadillo.h>   
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]
List MVTKalffbsC(arma::mat Y, List filter_output){
  
  // Model specification
  arma::mat GG(2,2);
  GG(0,0) = 1;
  GG(1,1) = 1;
  GG(0,1) = 0;
  GG(1,0) = 0;
  
  List theta_t = filter_output["theta_t"];
  List P_t = filter_output["P_t"];
  List theta_t_1 = filter_output["theta_t_1"];
  List P_t_1 = filter_output["P_t_1"];
  
  int n = Y.n_rows;
  
  List theta(n+1) ;
  List h_t(n+1) ;
  List H_t(n+1) ;
  
  // initializing the values
  h_t[n] = theta_t[n];      // m(n)
  H_t[n] = P_t[n];          // C(n)
  
  // dimension of state vector
  int p = as<vec>(h_t[n]).n_rows;
  
  // initializing theta
  arma::mat L = chol(as<mat>(H_t[n]));
  theta[n] = L*as<vec>(Rcpp::rnorm(p, 0, 1))  + as<vec>(h_t[n]);

  
  for(int i = 1; i <= n; i++){    // FFBS Iterations Start
    // Rcout << i << std::endl; // !!! REMOVE JUST FOR TEST
    arma::vec tt = theta_t[n-i];
    arma::vec tt1 = theta_t_1[n-i];
    arma::vec t = theta[n-i+1];
    arma::mat Pt = P_t[n-i];
    arma::mat Pt1 = P_t_1[n-i];
    
    h_t[n-i] = tt+Pt*GG.t()*inv(Pt1)*(t-tt1);
    H_t[n-i] = Pt-Pt*GG.t()*inv(Pt1)*GG*Pt;
    
    
    //Sampling the parameters
    arma::vec ht = h_t[n-i];
    arma::mat Ht = H_t[n-i];
    
    // Numeric Stability;
    if(Ht(0,0) <= 1e-5){
      Ht(0,0) = 1e-5;
    }
    
    // if( n-i<2){
      //   Ht(0,0) = 1000;
      //   Ht(1,1) = 1000;
      //   Ht(0,1) = 0;
      //   Ht(1,0) = 0;
      // }
    H_t[n-i] = Ht;
    
    // Sampling theta
    L = chol(Ht);
    theta[n-i] = L*as<vec>(Rcpp::rnorm(p, 0, 1)) + ht;
    
  }   // FFBS Iterations Finish
  
  
  
  List out;
  // out["H_t"] = H_t;
  // out["h_t"] = h_t;
  out["theta"] = theta;
  return(out);
  
}


