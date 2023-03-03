#include <RcppArmadillo.h>   
using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]


// When cor = 0 we have that the kalman filter of the joint is giving the same values
// as when running two independent kalman filters



List MVTKalFiltGarchC(arma::mat Y, arma::mat GG, arma::mat FF, arma::mat W, 
                      arma::vec m0, arma::mat C0, arma::mat garch, double cor, arma::mat Cov){
  
  int n = Y.n_rows;
  
  List theta_t(n+1) ;       // theta_t|t
  List P_t(n+1)   ;         // P_t|t
  List theta_t_1(n) ;       // theta_t|t-1
  List P_t_1(n) ;           // P_t|t-1
  List Q_t(n)   ;
  List y_t_1(n) ;           // hat Y_t|t-1
  List e_t(n)  ;            // Y_t-Y_t|t-1
  List k_t(n) ;             // Kalman gain
  List s_0t(n+1) ;             // Garch Errors for first series
  List s_1t(n+1) ;             // Garch Errors for first series
  List s_2t(n+1) ;             // Garch Errors for second series
  List s_3t(n+1) ;             // Garch Errors for second series
  List s_4t(n+1) ;             // Garch Errors for second series
  
  theta_t[0] = m0;                               // m0
  theta_t_1[0] = GG*as<vec>(theta_t[0]);         // a1
  
  // Rcout << m0 << std::endl;
  
  P_t[0] = C0;                                   // C0
  P_t_1[0] = GG*as<mat>(P_t[0])*GG.t()+W;        // R1
  
  double a00 = garch(0,0); double a01 = garch(0,1); double b01 = garch(0,2);
  double a10 = garch(1,0); double a11 = garch(1,1); double b11 = garch(1,2);
  double a20 = garch(2,0); double a21 = garch(2,1); double b21 = garch(2,2);
  double a30 = garch(3,0); double a31 = garch(3,1); double b31 = garch(3,2);
  double a40 = garch(4,0); double a41 = garch(4,1); double b41 = garch(4,2);
  
  
  s_0t[0] = a00;
  s_1t[0] = a10;
  s_2t[0] = a20;
  s_3t[0] = a30;
  s_4t[0] = a40;
  
  
  for(int i=0; i<n; i++){
    
    arma::vec t = theta_t[i];
    arma::vec t1 = theta_t_1[i];
    arma::mat Pt = P_t[i];
    arma::mat Pt1 = P_t_1[i];
    arma::mat FF_temp = FF;
    // Forecast Y
    y_t_1[i] = FF*t1;             // 1-step obsevartion
    arma::vec yt1 = y_t_1[i];
    // Rcout << yt1 << std::endl;
    
    
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
    double s0t = s_0t[i];
    double s1t = s_1t[i];
    double s2t = s_2t[i];
    double s3t = s_3t[i];
    double s4t = s_4t[i];
    
    arma::vec v = NumericVector::create(sqrt(s0t),sqrt(s1t),sqrt(s2t),sqrt(s3t),sqrt(s4t));
    arma::mat V = diagmat(v);         // matrix of the garch varainces
    
    
    arma::mat st = V*Cov*V;            // matrix of mvt garch errors
    
    // arma::mat R = eye<mat>(5,5);
    // // R(0,1) = R(1,0) = cor;            // matrix of constant correlation
    // R(0,1) = R(1,0) = R(0,2) = R(2,0) = R(0,3) = R(3,0) = R(0,4) = R(4,0) = 0;
    // R(1,2) = R(2,1) = R(1,3) = R(3,1) = R(1,4) = R(4,1) = 0;
    // R(2,3) = R(3,2) = R(2,4) = R(4,2) = 0;
    // R(3,4) = R(4,3) = 0;
    // 
    // arma::mat st = V*R*V;            // matrix of mvt garch errors
    // // Rcout << st << std::endl;
    
    
    
    Q_t[i] = FF_temp*Pt1*FF_temp.t() + st;       //  MSE of 1-step obs, innovations errors
    arma::mat Qt = Q_t[i];
    
    
    
    
    // Garch evolution of the series
    arma::vec innov = ((Y.row(i)).t()-t);
    // Rcout << innov << std::endl;
    if(NumericVector::is_na(Y(i,0)) != TRUE){
      s_0t[i+1] = a00+a01*(pow(Y(i,0),2)-2*Y(i,0)*1*t[0]+1*(Pt(0,0)+t[0]*t[0])*1)+b01*s1t;
    }else{
      s_0t[i+1] = a00+b01*s0t;
    }
    
    if(NumericVector::is_na(Y(i,1)) != TRUE){
      s_1t[i+1] = a10+a11*(pow(Y(i,1),2)-2*Y(i,1)*1*t[1]+1*(Pt(1,1)+t[1]*t[1])*1)+b11*s2t;
    }else{
      s_1t[i+1] = a10+b11*s1t;
    }
    
    if(NumericVector::is_na(Y(i,2)) != TRUE){
      s_2t[i+1] = a20+a21*(pow(Y(i,2),2)-2*Y(i,2)*1*t[2]+1*(Pt(2,2)+t[2]*t[2])*1)+b21*s2t;
    }else{
      s_2t[i+1] = a20+b21*s2t;
    }
    
    if(NumericVector::is_na(Y(i,3)) != TRUE){
      s_3t[i+1] = a30+a31*(pow(Y(i,3),2)-2*Y(i,3)*1*t[3]+1*(Pt(3,3)+t[3]*t[3])*1)+b31*s3t;
    }else{
      s_3t[i+1] = a30+b31*s3t;
    }
    
    if(NumericVector::is_na(Y(i,4)) != TRUE){
      s_4t[i+1] = a40+a41*(pow(Y(i,4),2)-2*Y(i,4)*1*t[4]+1*(Pt(4,4)+t[4]*t[4])*1)+b41*s4t;
    }else{
      s_4t[i+1] = a40+b41*s4t;
    }
    
    
    
    
    // Filtering
    theta_t[i+1] = t1+Pt1*FF_temp.t()*Qt.i()*et;    // Filter state param
    P_t[i+1] = Pt1 - Pt1*FF_temp.t()*Qt.i()*FF_temp*Pt1;  // MSE filter
    
    // Kalman Gain matrix
    k_t[i] = GG*Pt1*FF_temp.t()*Qt.i() ;
    
    //Rcout << Pt1 << std::endl; //!!!! THIS LINE JUST FOR CHECK
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












