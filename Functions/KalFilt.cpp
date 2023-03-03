#include <RcppArmadillo.h>   
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]
List KalFiltC1(arma::vec Y, arma::mat GG, arma::mat FF, double V, arma::mat W, arma::vec m0, arma::mat C0){
  
  int n = Y.size();
  
  List theta_t(n+1) ;       // theta_t|t
  List P_t(n+1)   ;         // P_t|t
  arma::vec t1 ;       // theta_t|t-1
  arma::mat Pt1 ;                    // P_t|t-1
  double k_t ;             // Kalman gain
  double yt1;
  
  theta_t[0] = m0;                               // m0
  t1 = GG*as<vec>(theta_t[0]);         // a1
  
  P_t[0] = C0;                                   // C0
  Pt1 = GG*as<mat>(P_t[0])*GG.t()+W;        // R1  
  
  
  
  for(int i=0; i<n; i++){

    arma::vec t = theta_t[i];
    arma::mat Pt = P_t[i];
    double Qt;
    double et;
    
    // Forecast Y
    // yt1 = (FF*t1)[0];             // 1-step obsevartion
    // double yt1 = y_t_1[i];


    // Update (Filter) state param and MSE
    if(NumericVector::is_na(Y[i]) != TRUE){
      // Rcout << Y[i]-FF*t1 << std::endl;

      theta_t[i+1] = t1 + Pt1*FF.t()/(FF*Pt1*FF.t() + V)*(Y[i]-FF*t1);    // Filter state param
      P_t[i+1] = Pt1 - (Pt1*FF.t()/(FF*Pt1*FF.t() + V))*(FF * Pt1);  // MSE filter

      // Kalman Gain matrix
      // k_t[i] = GG*Pt1*FF.t()/Qt ;

    }else{
      et = 0;
      Qt = arma::datum::inf;
      theta_t[i+1] = t1;
      P_t[i+1] = Pt1;

    }
    //Rcout << Pt1 << std::endl; //!!!! THIS LINE JUST FOR CHECK
    // Forecast state param
    if(i < (n-1)){
      t1 = GG * as<vec>(theta_t[i+1]);
      Pt1 = GG * as<mat>(P_t[i+1]) * GG.t() + W;
    }

  }
  
  
  List out;
  out["theta_t"] = theta_t;
  out["P_t"] = P_t;
  // out["theta_t_1"] = theta_t_1;
  // out["P_t_1"] = P_t_1;
  // out["y_t_1"] = y_t_1;
  // out["Q_t"] = Q_t;
  // out["e_t"] = e_t;
  // out["FF"] = FF;
  // out["k_t"] = k_t;
  // out["V"] = V;
  // out["GG"] = GG;
  // out["FF"] = FF;
  // out["Y"] = Y;
  // out["W"] = W;
  
  return(out);
}











