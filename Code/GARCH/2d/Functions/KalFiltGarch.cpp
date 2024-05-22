
#include <RcppArmadillo.h>   


using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]


// When cor = 0 we have that the kalman filter of the joint is giving the same values
// as when running two independent kalman filters



List MVTKalFiltGarchC(arma::mat Y, arma::mat GG, arma::mat FF, arma::mat W, arma::vec m0, arma::mat C0, arma::mat garch, double cor){
  
  int n = Y.n_rows;
  int nt = 10; arma::mat psi = eye<mat>(2,2);
  
  List theta_t(n+1) ;       // theta_t|t
  List P_t(n+1)   ;         // P_t|t
  List theta_t_1(n) ;       // theta_t|t-1
  List P_t_1(n) ;           // P_t|t-1
  List Q_t(n)   ;
  List y_t_1(n) ;           // hat Y_t|t-1
  List e_t(n)  ;            // Y_t-Y_t|t-1
  List k_t(n) ;             // Kalman gain
  List s_1t(n+1) ;             // Garch Errors for first series
  List s_2t(n+1) ;             // Garch Errors for second series
  
  theta_t[0] = m0;                               // m0
  theta_t_1[0] = GG*as<vec>(theta_t[0]);         // a1
  
  // Rcout << m0 << std::endl;
  
  P_t[0] = C0;                                   // C0
  P_t_1[0] = GG*as<mat>(P_t[0])*GG.t()+W;        // R1
  
  double a00 = garch(0,0); double a10 = garch(1,0);
  double a01 = garch(0,1); double a11 = garch(1,1);
  double b01 = garch(0,2); double b11 = garch(1,2);
  
  s_1t[0] = a00;
  s_2t[0] = a10;
  
  for(int i=0; i<n; i++){
    
    arma::vec t = theta_t[i];
    arma::vec t1 = theta_t_1[i];
    arma::mat Pt = P_t[i];
    arma::mat Pt1 = P_t_1[i];
    arma::mat FF_temp = FF;
    // Forecast Y
    y_t_1[i] = FF*t1;             // 1-step obsevartion
    arma::vec yt1 = y_t_1[i];
 
    
    
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
      e_t[i] = et_temp;
    }else{
      // Rcout << "FALSE" << std::endl;
      e_t[i] = (Y.row(i)).t()-yt1;
    }
    
    arma::vec et = e_t[i];
    double s1t = s_1t[i];
    double s2t = s_2t[i];
    
    arma::vec v = NumericVector::create(sqrt(s1t),sqrt(s2t));
    arma::mat V = diagmat(v);         // matrix of the garch varainces
    
    arma::mat R = eye<mat>(2,2);
    R(0,1) = R(1,0) = cor;            // matrix of constant correlation
    
    arma::mat st = V*R*V;            // matrix of mvt garch errors
    // Rcout << st << std::endl;
    Q_t[i] = FF_temp*Pt1*FF_temp.t() + st;       //  MSE of 1-step obs, innovations errors
    arma::mat Qt = Q_t[i];
    
    
    // Garch evolution of the series
    arma::vec innov = ((Y.row(i)).t()-t);
    // Rcout << innov << std::endl;
    if(NumericVector::is_na(Y(i,0)) != TRUE){
      s_1t[i+1] = a00+a01*(pow(Y(i,0),2)-2*Y(i,0)*1*t[0]+1*(Pt(0,0)+t[0]*t[0])*1)+b01*s1t;
    }else{
      s_1t[i+1] = a00+b01*s1t;
    }
    
    if(NumericVector::is_na(Y(i,1)) != TRUE){
      s_2t[i+1] = a10+a11*(pow(Y(i,1),2)-2*Y(i,1)*1*t[1]+1*(Pt(1,1)+t[1]*t[1])*1)+b11*s2t;
    }else{
      s_2t[i+1] = a10+b11*s2t;
    }
    
    // a0+a1*(pow(Y[i],2)-2*Y[i]*FF*t+FF*(Pt+t*t.t())*FF.t())+b1*st;
    
    
    // Filtering
    theta_t[i+1] = t1+Pt1*FF_temp.t()*Qt.i()*et;    
    P_t[i+1] = Pt1 - Pt1*FF_temp.t()*Qt.i()*FF_temp*Pt1;  
    
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
  out["k_t"] = k_t;
  out["s_1t"] = s_1t;
  out["s_2t"] = s_2t;
  out["GG"] = GG;
  out["FF"] = FF;
  out["Y"] = Y;
  out["W"] = W;
  
  
  return(out);
}










