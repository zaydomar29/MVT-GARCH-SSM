#include <RcppArmadillo.h>   
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]
List MVTKalFiltC(arma::mat Y, arma::mat GG, arma::mat FF, arma::mat V, arma::mat W, arma::vec m0, arma::mat C0){
  
  int n = Y.n_rows;
  
  List theta_t(n+1) ;       // theta_t|t
  List P_t(n+1)   ;         // P_t|t
  List theta_t_1(n) ;       // theta_t|t-1
  List P_t_1(n) ;           // P_t|t-1
  List Q_t(n)   ;
  List y_t_1(n) ;           // hat Y_t|t-1
  List e_t(n)  ;            // Y_t-Y_t|t-1
  List k_t(n) ;             // Kalman gain

  theta_t[0] = m0;                               // m0
  theta_t_1[0] = GG*as<vec>(theta_t[0]);         // a1

  P_t[0] = C0;                                   // C0
  P_t_1[0] = GG*as<mat>(P_t[0])*GG.t()+W;        // R1
  
  
  

  for(int i=0; i<n; i++){

    arma::vec t = theta_t[i];
    arma::vec t1 = theta_t_1[i];
    arma::mat Pt = P_t[i];
    arma::mat Pt1 = P_t_1[i];
    arma::mat FF_temp = FF;

    // Forecast Y
    y_t_1[i] = FF*t1;             // 1-step obsevartion
    arma::vec yt1 = y_t_1[i];
    
    if((Y.row(i)).has_nan() == TRUE){
      arma::vec et_temp = (Y.row(i)).t();
      for(int j = 0; j<Y.n_cols; j++ ){
        if(NumericVector::is_na(Y(i,j)) == TRUE){
          et_temp[j] = 0;
          FF_temp(j,j) = 0;
        }else{
          et_temp[j] = Y(i,j)-yt1[j];
        }
      }
      e_t[i] = et_temp;
    }else{
      // Rcout << "FALSE" << std::endl;
      e_t[i] = (Y.row(i)).t()-yt1;
    }
    
    arma::vec et = e_t[i];
    
    
    
    Q_t[i] = FF_temp*Pt1*FF_temp.t() + V;       //  MSE of 1-step obs, innovations errors
    arma::mat Qt = Q_t[i];
    
    theta_t[i+1] = t1+Pt1*FF_temp.t()*Qt.i()*et;    // Filter state param
    P_t[i+1] = Pt1 - Pt1*FF_temp.t()*Qt.i()*FF_temp*Pt1;  // MSE filter
    
    // Kalman Gain matrix
    k_t[i] = GG*Pt1*FF_temp.t()*Qt.i() ;
    

    
    // Forecast state param
    if(i < (n-1)){
      theta_t_1[i+1] = GG * as<vec>(theta_t[i+1]);
      P_t_1[i+1] = GG * as<mat>(P_t[i+1]) * GG.t() + W;
    }

  }


  List out;
  out["theta_t"] = theta_t;
  out["P_t"] = P_t;
  out["theta_t_1"] = theta_t_1;
  out["P_t_1"] = P_t_1;
  out["y_t_1"] = y_t_1;
  out["Q_t"] = Q_t;
  out["e_t"] = e_t;
  out["FF"] = FF;
  out["k_t"] = k_t;
  out["V"] = V;
  out["GG"] = GG;
  out["FF"] = FF;
  out["Y"] = Y;
  out["W"] = W;

  return(out);
}

