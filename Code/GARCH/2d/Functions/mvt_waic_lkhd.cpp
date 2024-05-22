#include <RcppArmadillo.h>
#include <MvtKalFiltGarch.cpp>
#include <dmvnormC.cpp>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]

List MVT_waic(arma::mat par,arma::mat Y, std::string mode, Rcpp::Nullable<List> post_theta = R_NilValue){
  
  List out ;
  
  // Model specification
  
  int S = par.n_rows;
  
  // nlk_mat and lk_mat contain the lkhd and lg-lkhd for Y_i given the 
  // j^th posterior samples
  arma::mat nlk_mat(S,Y.n_rows);
  arma::mat lk_mat(S,Y.n_rows);
  
  double nlk;
  double lk;
  double waic = 0;
  
  arma::mat ST_mat(S,Y.size()+1);
  
  if(post_theta.isNotNull() & mode == "GJR"){
    Rcout << "GJR" << std::endl;
    List post_(post_theta.get());         // post_ is the matrix containing Post theta
    for(int j=0; j<S; j++){  
      
      // Initializing the Posterior Parameters using the j^th posterior values
      double a00 = par(j,0); double a10 = par(j,4);
      double a01 = par(j,1); double a11 = par(j,5);
      double g01 = par(j,2); double g11 = par(j,6);
      double b01 = par(j,3); double b11 = par(j,7);
      double cor = par(j,8);
      
      
      arma::mat theta(Y.n_rows+1,Y.n_cols);
      theta = as<arma::mat>(post_[j]);

      // Starting Variance recursion
      arma::vec s_1t(Y.n_rows+1) ;
      arma::vec s_2t(Y.n_rows+1) ;
      s_1t[0] = a00;
      s_2t[0] = a10;
      
      
      arma::vec v = NumericVector::create(sqrt(s_1t[0]),sqrt(s_2t[0]));
      arma::mat V = diagmat(v);         // matrix of the garch varainces
      
      arma::mat R = eye<mat>(2,2);
      R(0,1) = R(1,0) = cor;            // matrix of constant correlation
      
      arma::mat st = V*R*V;             // mvt garch matrix
      
      
      
      for(int i=0; i<Y.n_rows; i++){              
        // This section calculates the lkhd and the variances for Y's given 
        // j^th posterior samples
        
        // Below is the setup for GJR-garch if desired
        nlk_mat(j,i) = exp(dmvnormC((Y.row(i)).t(),(theta.row(i)).t(),st));
        lk_mat(j,i) = dmvnormC((Y.row(i)).t(),(theta.row(i)).t(),st);
        
        double innov1= (Y(i,0)-theta(i+1,0));
        double innov2= (Y(i,1)-theta(i+1,1));
        if(innov1 >= 0){
          s_1t[i+1] = a00+a01*pow(Y(i,0)-theta(i,0),2)+b01*s_1t[i];
        }
        if(innov1 < 0){
          s_1t[i+1] = a00+g01*pow(Y(i,0)-theta(i,0),2)+b01*s_1t[i];
        }
        if(innov2 >= 0){
          s_2t[i+1] = a10+a11*pow(Y(i,1)-theta(i,1),2)+b11*s_2t[i];
        }
        if(innov2 < 0){
          // ST_mat(j,i+1) = a0+(g)*pow(Y[i]-theta[i+1],2)+b1*ST_mat(j,i);
          s_2t[i+1] = a10+g11*pow(Y(i,1)-theta(i,1),2)+b11*s_2t[i];
        }
        
        
        v = NumericVector::create(sqrt(s_1t[i+1]),sqrt(s_2t[i+1]));
        V = diagmat(v);                   // matrix of the garch varainces
        R = eye<mat>(2,2);
        R(0,1) = R(1,0) = cor;            // matrix of constant correlation
        st = V*R*V;                       // mvt garch matrix
      }
    }
  
  // Rcout << nlk_mat << std::endl;
  
    for(int i=0; i<Y.n_rows; i++){
      lk = 0;
      nlk = 0;
        for(int j=0; j<S; j++){
          nlk += nlk_mat(j,i);
          lk += lk_mat(j,i);
          }
      waic += -log(nlk/S)+2*lk/S;           // Calculating the WAIC in this section
      }
  } 
  
  
  // out["Qt"] = KF_mat_Qt;
  // out["ft"] = KF_mat_ft;
  out["waic"] = 2*waic;
  
  return(out);
}
  
  // /***R
  // LG_lgpostC(par, Y)
  // */
  
  
