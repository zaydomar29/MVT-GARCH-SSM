#include <RcppArmadillo.h>
#include <MvtKalFilt.cpp>
#include <dmvnormC.cpp>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]

List MVT_waic1(List V,List W, arma::mat Y, std::string mode, arma::mat FF, arma::mat GG, arma::vec m0, arma::mat C0){
  
  List out ;
  
  // Model specification
  
  int S = V.length();
  int p = Y.n_cols;
  // nlk_mat and lk_mat contain the lkhd and lg-lkhd for Y_i given the 
  // j^th posterior samples
  arma::mat nlk_mat(S,Y.n_rows);
  arma::mat lk_mat(S,Y.n_rows);
  
  double nlk;
  double lk;
  double waic = 0;
  
  // arma::mat FF = eye(p,p);
  // arma::mat GG = eye(4,4);
  // arma::vec m0 = (Y.row(0)).t();
  // arma::mat C0 = 10*eye(4,4);
  
  
  if(mode == "Std"){
    Rcout << "Std" << std::endl;
    for(int j=0; j<S; j++){  
      
      // Run kalman filter
      List filt_val = MVTKalFiltC(Y,GG,FF,as<arma::mat>(V[j]),as<arma::mat>(W[j]),m0,C0);
      List y_t_1 = filt_val["y_t_1"];
      List Q_t = filt_val["Q_t"];
      
      for(int i=0; i<Y.n_rows; i++){              
        // This section calculates the lkhd and the variances for Y's given 
        // j^th posterior samples
        
        arma::vec mean_y = as<arma::vec>(y_t_1[i]);
        arma::mat var_y = as<arma::mat>(Q_t[i]);
        
        // Below is the setup for GJR-garch if desired
        if(Y.row(i).has_nan() == TRUE){
        nlk_mat(j,i) = 0;
        lk_mat(j,i) = 0;
        }else{
        // Rcout << var_y << std::endl;
        // Rcout << exp(dmvnormC((Y.row(i)).t(),mean_y,var_y)) << std::endl;
        nlk_mat(j,i) = exp(dmvnormC((Y.row(i)).t(),mean_y,var_y));
        lk_mat(j,i) = dmvnormC((Y.row(i)).t(),mean_y,var_y);
        }
      } 
    }
        // Rcout << nlk_mat << std::endl;
        
        for(int i=0; i<Y.n_rows; i++){
        lk = 0;
        nlk = 0;
        
        if(Y.row(i).has_nan() == TRUE){
        waic += 0;
        }else{
        for(int j=0; j<S; j++){
          nlk += nlk_mat(j,i);
          lk += lk_mat(j,i);
        }
          waic += -log(nlk/S)+2*lk/S;           // Calculating the WAIC in this section
        }
        
        
        }
  } 
        
        
        out["waic"] = 2*waic;
        
        return(out);
}
      
     
