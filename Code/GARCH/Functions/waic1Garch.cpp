#include <RcppArmadillo.h>
#include <MvtKalFiltGarchC.cpp>
#include <dmvnormC.cpp>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]

List MVT_Garch_waic1(List par,List cor,List W, arma::mat Y, arma::mat FF, arma::mat GG, arma::vec m0){
  
  List out ;
  
  // Model specification
  
  int S = par.size();
  Rcout << S<<std::endl;
  // nlk_mat and lk_mat contain the lkhd and lg-lkhd for Y_i given the 
  // j^th posterior samples
  arma::mat nlk_mat(S,Y.n_rows);
  arma::mat lk_mat(S,Y.n_rows);
  
  double nlk;
  double lk;
  double waic = 0;
  
  int p = GG.n_cols;
  int k = Y.n_cols;
  arma::mat C0 = 10*eye(p,p);
  
  
  
  
  arma::mat garch(k,3);
  
  
  for(int j=0; j<S; j++){
    
    // Initializing the Posterior Parameters using the j^th posterior values
    List filt_val;
    
    garch = as<arma::mat>(par(j));
    if( k > 2){
      arma::mat rho(k,k);
      rho = as<arma::mat>(cor(j));
      // Run kalman filter
      filt_val = MVTKalFiltGarch(Y,GG,FF,as<arma::mat>(W[j]),m0,C0,garch,rho);
    }else{
      arma::mat rho(1,1);
      rho = as<arma::mat>(cor(j));
      // Run kalman filter
      filt_val = MVTKalFiltGarch(Y,GG,FF,as<arma::mat>(W[j]),m0,C0,garch,rho);
    }
    
    
    
    List y_t_1 = filt_val["y_t_1"];
    List Q_t = filt_val["Q_t"];
    
    for(int i=0; i<Y.n_rows; i++){
      // This section calculates the lkhd and the variances for Y's given
      // j^th posterior samples
      
      arma::vec mean_y = as<arma::vec>(y_t_1[i]);
      arma::mat var_y = as<arma::mat>(Q_t[i]);
      
      if(Y.row(i).has_nan() == TRUE){
        nlk_mat(j,i) = 0;
        lk_mat(j,i) = 0;
      }else{
        // Rcout << var_y << std::endl;
        nlk_mat(j,i) = exp(dmvnormC((Y.row(i)).t(),mean_y,var_y));
        lk_mat(j,i) = dmvnormC((Y.row(i)).t(),mean_y,var_y);
      }
    }
  }
  
  
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
  
  out["waic"] = 2*waic;
  
  return(out);

        
}
      
      
      
      
