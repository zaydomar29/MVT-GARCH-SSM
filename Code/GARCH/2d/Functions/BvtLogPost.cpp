#include <RcppArmadillo.h>   
#include <dmvnormC.cpp>   
#include <diwishart.cpp>   
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]


// When cor = 0 we have that the kalman filter of the joint is giving the same values
// as when running two independent kalman filters



double MVT_lgpostSS(arma::mat par, double cor ,arma::mat W, arma::mat Y, arma::mat theta, int nt, arma::mat psi){
  
  // arma:: mat psi = eye<mat>(2,2);
  // List out;
  // Initializing the Parameters
  double a00 = par(0,0); double a10 = par(1,0);
  double a01 = par(0,1); double a11 = par(1,1);
  double b01 = par(0,2); double b11 = par(1,2);
  
  
  // Initializing the Parameters
  arma::vec s_1t(Y.n_rows+1) ;
  arma::vec s_2t(Y.n_rows+1) ;
  s_1t[0] = a00;
  s_2t[0] = a10;
  
  
  
  
  arma::vec v = NumericVector::create(sqrt(s_1t[0]),sqrt(s_2t[0]));
  arma::mat V = diagmat(v);         // matrix of the garch varainces
  
  arma::mat R = eye<mat>(2,2);
  R(0,1) = R(1,0) = cor;            // matrix of constant correlation
  
  arma::mat st = V*R*V;             // mvt garch matrix
  
  
  
  // Initializing log-posterior
  double lp = 0;
  double lkhd_obs = 0;
  double lkhd_state= 0;
  
  
  for(int i=0; i<Y.n_rows; i++){
    lkhd_obs = dmvnormC((Y.row(i)).t(),(theta.row(i+1)).t(), st);
    lkhd_state = dmvnormC((theta.row(i+1)).t(),(theta.row(i)).t(),W);
    lp += lkhd_obs+lkhd_state;
    
    
    // Innovations
    arma::vec innov = (Y.row(i)-theta.row(i+1)).t();
    
    // First Series
    if(NumericVector::is_na(Y(i,0)) == false){
      s_1t[i+1] = a00+a01*pow(Y(i,0)-theta(i+1,0),2)+b01*s_1t[i];
    }else{
      s_1t[i+1] = a00+b01*s_1t[i];
    }
    
    // Second Series
    if(NumericVector::is_na(Y(i,1)) == false){
      s_2t[i+1] = a10+a11*pow(Y(i,1)-theta(i+1,1),2)+b11*s_2t[i];
    }else{
      s_2t[i+1] = a10+b11*s_2t[i]; 
    }
    
    
    v = NumericVector::create(sqrt(s_1t[i+1]),sqrt(s_2t[i+1]));
    V = diagmat(v);                   // matrix of the garch varainces
    R = eye<mat>(2,2);
    R(0,1) = R(1,0) = cor;            // matrix of constant correlation
    
    st = V*R*V;             // mvt garch matrix
    
  }
  double prior = R::dcauchy(a00,0,1,1)+R::dcauchy(a10,0,1,1)+R::dcauchy(a01,0,1,1)+
    R::dcauchy(a11,0,1,1)+R::dcauchy(b01,0,1,1)+
    R::dcauchy(b11,0,1,1)+diwishart(W,psi,nt);
  // Rcout << diwishart(W,psi,nt)  << std::endl;
  // Rcout << R::dgamma(1/w1,2,0.5,1)/pow(w1,2) << std::endl; // !!! REMOVE JUST FOR TEST
  
  lp += prior;
  
  // out["s1"] = s_1t;
  // out["s2"] = s_2t;
  // out["lp"] = lp;
  
  return(lp);
}
