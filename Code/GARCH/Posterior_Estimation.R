##################################################################
##################################################################
##############                               #####################
##############     MH-MCMC estimation        #####################
##############                               #####################
##################################################################
##################################################################
Rcpp::sourceCpp("./Functions/MvtKalFiltGarchC.cpp")
Rcpp::sourceCpp("./Functions/FFBS.cpp")
Rcpp::sourceCpp("./Functions/MvtLkhd.cpp")
Rcpp::sourceCpp("./Functions/diwishart.cpp")
library(MCMCpack)


L1 = chol(diag(c(0.5,0.01,0.01)))
S.t = diag(1,dim(GG)[1]); n.t = 10;

burn = 10000
thin = 1
sims = (5000)*thin+burn
sims

## Store posterior samples
S.t = diag(1,dim(GG)[1]); n.t = 10;


S_t.samp = list()

predicted = list()
sampled_theta = list()


garch_par_samp = list(matrix(c(1,0,0,1,0,0),nr = 2,byr=T))
old_param = garch_par_samp[[1]]

W.samp = list()
W.old = diag(1,dim(GG)[2])

cor.samp = list()
if(ncol(Y)==2){
  cor.old = 0
  V.old = matrix(0)
}else{
  cor.old = diag(1,ncol(Y))
  V.old = tcrossprod(cor.old/apply(cor.old, 1,norm, "2"))
}




out = MVTKalFiltGarch(Y=Y,GG=GG,FF=FF,W=W.old,m0=m0,C0=C0,garch=old_param,cor=V.old)
theta = do.call(rbind,lapply(ffbsC(Y, out)$theta,t))




## Start MCMC
for(it in 1:sims){
  if(it %% 100 == 0){print(it)}
  if(it %% 100 == 0){print(old_param)}
  if(it %% 100 == 0){print(cor.old)}
  if(it %% 100 == 0){print(W.old)}
  
  
  n = nrow(Y);
  
  for(j in 1:ncol(Y)){
    prop = t(c(old_param[j,1],old_param[j,2],old_param[j,3])+L1%*%rnorm(3,0,1))
    while(sum(prop[2],prop[3]) > 1 || any(prop<0)){
      prop = t(c(old_param[j,1],old_param[j,2],old_param[j,3])+L1%*%rnorm(3,0,1))
    }
    
    new_param = old_param
    new_param[j,1:3] = c(prop[1],prop[2],prop[3])
    rho = MVT_Lkhd(par=new_param,cor=V.old,Y=Y,W=W.old,theta=theta,FF = FF,GG=GG)-
      MVT_Lkhd(par=old_param,cor=V.old,Y=Y,W=W.old,theta=theta,FF = FF,GG=GG)+
      sum(dcauchy(prop, log = T))-sum(dcauchy(old_param[j,], log = T))
    
    if(log(runif(1))<rho){
      old_param[j,1:3] = c(prop[1],prop[2],prop[3])
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
        garch_par_samp[[(it-burn-1)/thin]] = old_param
      }
    } else {
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 1 ){
        garch_par_samp[[(it-burn-1)/thin]] = old_param
      }
    }
  }
  
  if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
    garch_par_samp[[(it-burn-1)/thin]] = old_param
  }
  
  mu_theta = t(GG %*% t(theta[-n,]))
  SS.t = crossprod(theta[-1,]-mu_theta) + S.t
  df.t = n+n.t
  W.old = riwish(df.t,SS.t)
  if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
    W.samp[[round((it-burn-1)/thin)]] = W.old
  }
  
  
  if(dim(Y)[2]>2){
    
    # Use this reparametrization only when the dimension of Y is greater than 2
    
    j = sample(p-1,1)
    k = sample(j:p,1)
    cor.new = cor.old
    cor.new[j,k] = round(cor.old[j,k]+rnorm(1,0,0.1),5)
    cor.new[j,j:p] = cor.new[j,j:p]/norm(cor.new[j,j:p],"2")
    
    d1=d2=0
    if(j==k){
      # positivity constraint
      cor.new[j,k] = rtruncnorm(1,a=0,mean=cor.old[j,k],sd=0.1)
      cor.new[j,j:p] = cor.new[j,j:p]/norm(cor.new[j,j:p],"2")
      
      d1=dtruncnorm(cor.new[j,k],a=0,mean=cor.old[j,k],sd=0.1)
      d2=dtruncnorm(cor.old[j,k],a=0,mean=cor.new[j,k],sd=0.1)
    }
    
    V.new = tcrossprod(round(cor.new,5))
    rho = MVT_Lkhd(par=old_param,cor=V.new,Y=Y,W=W.old,theta=theta,FF=FF,GG=GG)+d2-
      MVT_Lkhd(par=old_param,cor=V.old,Y=Y,W=W.old,theta=theta,FF=FF,GG=GG)-d1
    
    if(log(runif(1))<rho){
      cor.old = cor.new
      V.old = V.new
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
        cor.samp[[(it-burn-1)/thin]] = cor.old
      }
    } else {
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 1 ){
        cor.samp[[(it-burn-1)/thin]] = cor.old
      }
    }
  }else{
    # When we have the 2-d case the correlation parameter, V, is pretty simple
    V.new = matrix(V.old+runif(1,-0.1,0.1))
    while( V.new > 1 || V.new < -1){
      V.new = matrix(V.old+runif(1,-0.1,0.1))
    }
    rho = MVT_Lkhd(par=old_param,cor=V.new,Y=Y,W=W.old,theta=theta,FF=FF,GG=GG)-
      MVT_Lkhd(par=old_param,cor=V.old,Y=Y,W=W.old,theta=theta,FF=FF,GG=GG)
    
    if(log(runif(1))<rho){
      cor.old = V.new
      V.old = V.new
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
        cor.samp[[(it-burn-1)/thin]] = cor.old
      }
    } else {
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 1 ){
        cor.samp[[(it-burn-1)/thin]] = cor.old
      }
    }
    
    
  }
  
  
  out = MVTKalFiltGarch(Y=Y,GG=GG,FF=FF,W=W.old,m0=m0,C0=C0,garch=old_param,cor=V.old)
  theta = do.call(rbind,lapply(ffbsC(Y, out)$theta,t))
  
  if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
    S_t.samp[[(it-burn-1)/thin]] = out$S_t
    predicted[[(it-burn-1)/thin]] = do.call(rbind, lapply(out$y_t_1,t))
    sampled_theta[[(it-burn-1)/thin]] = theta
  }
  
  
}

