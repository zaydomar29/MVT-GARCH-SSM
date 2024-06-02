#include <RcppArmadillo.h>   
#include <myround.cpp>
#include <makeSymmetric.cpp>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]


// When cor = 0 we have that the kalman filter of the joint is giving the same values
// as when running two independent kalman filters



List MVTKalFiltGarch(arma::mat Y, arma::mat GG, arma::mat FF, arma::mat W, arma::vec m0, arma::mat C0, arma::mat garch,
                      arma::mat cor){ //  ,Rcpp::Nullable<int> forecast
  
  int n = Y.n_rows;
  int p = Y.n_cols;
  
  List theta_t(n+1) ;        // theta_t|t
  List P_t(n+1)   ;          // P_t|t
  List theta_t_1(n) ;        // theta_t|t-1
  List P_t_1(n) ;            // P_t|t-1
  List V_t(n);               // Obs-variance at time t        
  
  arma::vec t(p);            // theta_t|t
  arma::mat Pt(p,p);         // P_t|t
  arma::vec t1(p);           // theta_t|t-1
  arma::mat Pt1(p,p);        // P_t|t-1
  arma::vec e_t(p);          // Y_t-Y_t|t-1
  arma::mat S_t(n+1,p);
  
  
  List Q_t(n);
  List y_t_1(n);             // hat Y_t|t-1
  
  
  // // // // Forecast variables
  // if(forecast.isNotNull()){
  //   List a_t_k(n-forecast+1);
  //   List f_t_k(n-forecast+1);
  //   List Q_t_k(n-forecast+1);
  //   List R_t_k(n-forecast+1);
  // }
  
  
  
  
  t = m0;                    // m0
  theta_t[0] = t;
  t1 = GG*t;                 // a1
  theta_t_1[0] = t1;

  Pt = C0;                   // C0
  P_t[0] = Pt;
  
  
  Pt1 = GG*Pt*GG.t()+W;      // R1
  P_t_1[0] = Pt1;
  
  S_t.row(0) = (garch.col(0)).t();
  
  for(int i=0; i<n; i++){
    arma::mat FF_temp = FF;
    
    // Forecast Y
    y_t_1[i] = FF*t1;             // 1-step obsevartion
    arma::vec yt1 = y_t_1[i];
    // HERE HERE
    // Update (Filter) state param and MSE
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
      e_t = et_temp;
    }else{
      // Rcout << "FALSE" << std::endl;
      e_t = (Y.row(i)).t()-yt1;
    }

    
    arma::vec v = (S_t.row(i)).t();
    
    arma::mat V = pow(diagmat(v),0.5);         // matrix of the garch varainces
    arma::mat R = eye<mat>(p,p);
    if(p == 2){
      R(0,1) = R(1,0) = cor(0,0);            // matrix of constant correlation
    }
    if(p>2){
      R = cor;                         // matrix of constant correlation
    }
    
    arma::mat st = V*R*V;                      // matrix of mvt garch errors / Observation variance
    // Rcout << st << std::endl;
    V_t[i] = st;                            
    arma::mat Qt(p,p);
    Qt = makeSymmetric(FF_temp*Pt1*FF_temp.t()+st);           //  MSE of 1-step obs, innovations errors
    Q_t[i] = Qt;
    
    for(int j = 0; j<p; j++ ){
      if(NumericVector::is_na(Y(i,j)) != TRUE){

        arma::vec a = FF_temp*t;
        arma::mat b = (FF_temp*(Pt+t*t.t())*FF_temp.t());
        
        S_t(i+1,j) = garch(j,0)+garch(j,1)*(pow(Y(i,j),2)-2*Y(i,j)*a[j]+b(j,j))+garch(j,2)*S_t(i,j);
        
      }else{
        S_t(i+1,j) = garch(j,0)+garch(j,2)*S_t(i,j);
      }
    }
    
    // Filtering
    arma::mat Pt1F = Pt1*FF_temp.t();
    arma::mat QtI = makeSymmetric(inv(Qt));
    t = myround(t1+Pt1F*QtI*e_t,5);
    theta_t[i+1] = t;
    Pt = makeSymmetric(Pt1 - Pt1F*QtI*Pt1F.t());
    
    P_t[i+1] = Pt;
    
    // Kalman Gain matrix
    // k_t[i] = GG*Pt1*FF_temp.t()*Qt.i() ;
    //

    // Forecast state param
    if(i < (n-1)){
      
      t1 = GG*t;
      theta_t_1[i+1] = t1;
      // Rcout << GG * as<vec>(theta_t[i+1]) << std::endl;
      Pt1 = GG*Pt*GG.t()+W;
      P_t_1[i+1] = Pt1;
      
      // Rcout << GG * as<mat>(P_t[i+1]) * GG.t() + W << std::endl;
    }
    
    



  }

  
  
  List out;
  out["GG"] = GG;
  out["theta_t"] = theta_t;
  out["P_t"] = P_t;
  out["theta_t_1"] = theta_t_1;
  out["P_t_1"] = P_t_1;
  out["y_t_1"] = y_t_1;
  out["Q_t"] = Q_t;
  out["S_t"] = S_t;
  out["V_t"] = V_t;
  
  
  
  
  
  return(out);
}











