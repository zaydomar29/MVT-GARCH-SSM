#include <RcppArmadillo.h>   
#include <dmvnormC.cpp>   
#include <diwishart.cpp>   
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]


// When cor = 0 we have that the kalman filter of the joint is giving the same values
// as when running two independent kalman filters



double MVT_lgpostSS(arma::mat par, double cor ,arma::mat W, arma::mat Y, 
                    arma::mat theta, int nt, arma::mat psi, arma::mat Cov){
  
  // arma:: mat psi = eye<mat>(5,5);
  // List out;
  // Initializing the Parameters
  double a00 = par(0,0); double a01 = par(0,1); double b01 = par(0,2);
  double a10 = par(1,0); double a11 = par(1,1); double b11 = par(1,2);
  double a20 = par(2,0); double a21 = par(2,1); double b21 = par(2,2);
  double a30 = par(3,0); double a31 = par(3,1); double b31 = par(3,2);
  double a40 = par(4,0); double a41 = par(4,1); double b41 = par(4,2);
  
  
  // Initializing the Parameters
  arma::vec s_0t(Y.n_rows+1) ;
  arma::vec s_1t(Y.n_rows+1) ;
  arma::vec s_2t(Y.n_rows+1) ;
  arma::vec s_3t(Y.n_rows+1) ;
  arma::vec s_4t(Y.n_rows+1) ;
  
  
  s_0t[0] = a00;
  s_1t[0] = a10;
  s_2t[0] = a20;
  s_3t[0] = a30;
  s_4t[0] = a40;
  
  
  
  
  arma::vec v = NumericVector::create(sqrt(s_0t[0]),sqrt(s_1t[0]),sqrt(s_2t[0]),
                                      sqrt(s_3t[0]),sqrt(s_4t[0]));
  
  arma::mat V = diagmat(v);         // matrix of the garch varainces
  
  // arma::mat R = eye<mat>(5,5);
  // // R(0,1) = R(1,0) = cor;            // matrix of constant correlation
  // R(0,1) = R(1,0) = R(0,2) = R(2,0) = R(0,3) = R(3,0) = R(0,4) = R(4,0) = 0;
  // R(1,2) = R(2,1) = R(1,3) = R(3,1) = R(1,4) = R(4,1) = 0;
  // R(2,3) = R(3,2) = R(2,4) = R(4,2) = 0;
  // R(3,4) = R(4,3) = 0;
  
  // arma::mat st = V*R*V;             // mvt garch matrix
  arma::mat st = V*Cov*V;            // matrix of mvt garch errors
  
  
  // Initializing log-posterior
  double lp = 0;
  double lkhd_obs = 0;
  double lkhd_state= 0;
  
  
  for(int i=0; i<Y.n_rows; i++){
    lkhd_obs = dmvnormC((Y.row(i)).t(),(theta.row(i+1)).t(), st);
    lkhd_state = dmvnormC((theta.row(i+1)).t(),(theta.row(i)).t(),W);
    lp += lkhd_obs+lkhd_state;
    
    
    // Sign of innovation
    arma::vec innov = (Y.row(i)-theta.row(i+1)).t();
    // Assymetry in first series
    if(NumericVector::is_na(Y(i,0)) == false){
      s_0t[i+1] = a00+a01*pow(Y(i,0)-theta(i+1,0),2)+b01*s_0t[i];
      // Rcout << "No NA, pos"  << std::endl;
    }else{
      s_0t[i+1] = a00+b01*s_0t[i];
    }
    
    // Assymetry in second series
    if(NumericVector::is_na(Y(i,1)) == false){
      s_1t[i+1] = a10+a11*pow(Y(i,1)-theta(i+1,1),2)+b11*s_1t[i];
    }else{
      s_1t[i+1] = a10+b11*s_1t[i]; 
    }
    
    // Assymetry in third series
    if(NumericVector::is_na(Y(i,2)) == false){
      s_2t[i+1] = a20+a21*pow(Y(i,2)-theta(i+1,2),2)+b21*s_2t[i];
    }else{
      s_2t[i+1] = a20+b21*s_2t[i]; 
    }
    
    // Assymetry in fourth series
    if(NumericVector::is_na(Y(i,3)) == false){
      s_3t[i+1] = a30+a31*pow(Y(i,3)-theta(i+1,3),2)+b31*s_3t[i];
    }else{
      s_3t[i+1] = a30+b31*s_3t[i]; 
    }
    
    // Assymetry in fifth series
    if(NumericVector::is_na(Y(i,4)) == false){
      s_4t[i+1] = a40+a41*pow(Y(i,4)-theta(i+1,4),2)+b41*s_4t[i];
    }else{
      s_4t[i+1] = a40+b41*s_4t[i]; 
    }
    
    v = NumericVector::create(sqrt(s_0t[i+1]),sqrt(s_1t[i+1]),sqrt(s_2t[i+1]),
                              sqrt(s_3t[i+1]),sqrt(s_4t[i+1]));
    
    V = diagmat(v);                   // matrix of the garch varainces
    // R = eye<mat>(5,5);
    // R(0,1) = R(1,0) = R(0,2) = R(2,0) = R(0,3) = R(3,0) = R(0,4) = R(4,0) = 0;
    // R(1,2) = R(2,1) = R(1,3) = R(3,1) = R(1,4) = R(4,1) = 0;
    // R(2,3) = R(3,2) = R(2,4) = R(4,2) = 0;
    // R(3,4) = R(4,3) = 0;
    
    // st = V*R*V;             // mvt garch matrix
    st = V*Cov*V;            // matrix of mvt garch errors
  }
  double prior = R::dcauchy(a00,0,1,1)+R::dcauchy(a01,0,1,1)+R::dcauchy(b01,0,1,1)+
                 R::dcauchy(a10,0,1,1)+R::dcauchy(a11,0,1,1)+R::dcauchy(b11,0,1,1)+
                 R::dcauchy(a20,0,1,1)+R::dcauchy(a21,0,1,1)+R::dcauchy(b21,0,1,1)+
                 R::dcauchy(a30,0,1,1)+R::dcauchy(a31,0,1,1)+R::dcauchy(b31,0,1,1)+
                 R::dcauchy(a40,0,1,1)+R::dcauchy(a41,0,1,1)+R::dcauchy(b41,0,1,1)+
                 diwishart(W,psi,nt);
  // Rcout << diwishart(W,psi,nt)  << std::endl;
  // Rcout << R::dgamma(1/w1,2,0.5,1)/pow(w1,2) << std::endl; // !!! REMOVE JUST FOR TEST
  
  lp += prior;
  
  // out["s1"] = s_1t;
  // out["s2"] = s_2t;
  // out["lp"] = lp;
  
  return(lp);
}
