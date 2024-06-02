#include <RcppArmadillo.h>   
#include <myround.cpp>
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]
List ffbsC(arma::mat Y, List filter_output){
  
  // Model specification
  List theta_t = filter_output["theta_t"];
  List P_t = filter_output["P_t"];
  List theta_t_1 = filter_output["theta_t_1"];
  List P_t_1 = filter_output["P_t_1"];
  
  
  arma::mat GG = as<mat>(filter_output["GG"]);
  int k = GG.n_rows;
  
  int n = Y.n_rows;
  List theta(n+1) ;
  arma::vec ht(k) ;
  arma::mat Ht(k,k) ;
  
  // initializing the values
  ht = as<vec>(theta_t[n]);      // m(n)
  Ht = as<mat>(P_t[n]);          // C(n)
  
  // dimension of state vector
  int p = ht.n_rows;
  
  // initializing theta
  // Rcout << Ht << std::endl;
  arma::mat L = chol(Ht);
  
  
  theta[n] = L*as<vec>(Rcpp::rnorm(p, 0, 1))+ht;
  for(int i = 1; i <= n; i++){    // FFBS Iterations Start
    
    arma::vec tt = myround(theta_t[n-i],10);
    // arma::vec tt = theta_t[n-i];
    arma::vec tt1 = myround(theta_t_1[n-i],10);
    // arma::vec tt1 = theta_t_1[n-i];
    arma::vec t = myround(theta[n-i+1],10);
    // arma::vec t = theta[n-i+1];
    
    // Rcout << n-i << std::endl;
    arma::mat Pt = P_t[n-i];
    arma::mat Pt1 = P_t_1[n-i];
    
    arma::mat Pt1I = inv(Pt1);
    ht = tt+Pt*GG.t()*Pt1I*(t-tt1);
    Ht = myround(Pt-(Pt*GG.t())*Pt1I*(GG*Pt),6);
    
    // Numeric Stability;
    // if(Ht(0,0) <= 1e-5){
    //   Ht(0,0) = 1e-5;
    // }
    
    // Sampling theta
    // Rcout << Ht << std::endl;
    // Rcout << i << std::endl;
    L = chol(Ht);
    theta[n-i] = L*as<vec>(Rcpp::rnorm(p, 0, 1)) + ht;
    
  }   // FFBS Iterations Finish
  
  
  
  List out;
  out["theta"] = theta;
  return(out);
  
}


