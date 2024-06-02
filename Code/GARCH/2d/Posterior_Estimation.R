##################################################################
##################################################################
##############                               #####################
##############     MH-MCMC estimation        #####################
##############                               #####################
##################################################################
##################################################################
Rcpp::sourceCpp("./Functions/BvtKalFiltGarch.cpp")
Rcpp::sourceCpp("./Functions/BvtKalFFBS.cpp")
Rcpp::sourceCpp("./Functions/BvtLogPost.cpp")
Rcpp::sourceCpp("./Functions/diwishart.cpp")
library(MCMCpack)
library(truncnorm)

nsims = 10e3
burn = 20000
thin = 80
sims = burn+thin*nsims+1


m0 = (Y[1,]); C0 = diag(20,2)

L1 = chol(diag(c(1,0.1,0.1)))
S.t = diag(1,2); n.t = 10;

## Store posterior samples
a0.samp= matrix(c(1,1),nc=2,nr=nsims, byrow = T)
a1.samp= matrix(c(0.1,0.3),nc=2,nr=nsims, byrow = T)
b1.samp= matrix(c(0.85,0.65),nc=2,nr=nsims, byrow = T)
W.samp = list()
rho.samp =c();
sampled_theta = list()
predicted = list()
var1 = list()
var2 = list()


## Old parameters
old_param = rbind(c(1,0.2,0.1),
                  c(1,0.2,0.1))
rho.old = 0.0
W.old = diag(c(1.5,2.5))

out = BvtKalFiltGarchC(Y=Y,GG=GG,FF=FF,W=W.old,m0=Y[1,],C0=C0,garch=old_param,cor=rho.old)
theta = do.call(rbind,lapply(BvtKalffbsC(Y, out)$theta,t))



old_lgpost = BvTLogPost(par = old_param,cor = rho.old,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(10,2))

acc.rate.a0 = c(0,0)
acc.rate.a1 = c(0,0)
acc.rate.b1 = c(0,0)

## Start MCMC
for(it in 1:sims){
  if(it %% 1000 == 0){
    print(it)
    print(old_param)
    print(rho.old)
    print(W.old)
    print(paste("acc.rate.a0:",acc.rate.a0/it))
    print(paste("acc.rate.a1:",acc.rate.a1/it))
    print(paste("acc.rate.b1:",acc.rate.b1/it))
    }
  
  n = nrow(Y);
  
  for(j in 1:2){
    a0.prop = rtruncnorm(1,a=0,b=Inf,mean=old_param[j,1],1)
    new_param = old_param
    new_param[j,1] = a0.prop
    rho = BvTLogPost(par = new_param,cor = rho.old,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(10,2))+
      log(dtruncnorm(old_param[j,1],a=0,b=Inf,mean=a0.prop,0.1))-
      BvTLogPost(par = old_param,cor = rho.old,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(10,2))-
      log(dtruncnorm(a0.prop,a=0,b=Inf,mean=old_param[j,1],0.1))
    
    
    if(log(runif(1))<rho){
      old_param[j,1] = a0.prop
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){a0.samp[(it-burn-1)/thin,j] = a0.prop}
      acc.rate.a0[j] = acc.rate.a0[j]+1
    }else{
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){a0.samp[(it-burn-1)/thin,j] = old_param[j,1]}
    } 
    
    a1.prop = rtruncnorm(1,a=0,b=1-old_param[j,3],mean=old_param[j,2],0.01)
    new_param = old_param
    new_param[j,2] = a1.prop
    rho = BvTLogPost(par = new_param,cor = rho.old,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(10,2))+
      log(dtruncnorm(old_param[j,2],a=0,b=1-old_param[j,3],mean=a1.prop,0.01))-
      BvTLogPost(par = old_param,cor = rho.old,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(10,2))-
      log(dtruncnorm(a1.prop,a=0,b=1-old_param[j,3],mean=old_param[j,2],0.01))
    
    
    if(log(runif(1))<rho){
      old_param[j,2] = a1.prop
      acc.rate.a1[j] = acc.rate.a1[j]+1
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){a1.samp[(it-burn-1)/thin,j] = a1.prop}
    }else{
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){a1.samp[(it-burn-1)/thin,j] = old_param[j,2]}
    }
    
    b1.prop = rtruncnorm(1,a=0,b=1-old_param[j,2],mean=old_param[j,3],0.01)
    new_param = old_param
    new_param[j,3] = b1.prop
    rho = BvTLogPost(par = new_param,cor = rho.old,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(10,2))+
      log(dtruncnorm(old_param[j,3],a=0,b=1-old_param[j,2],mean=b1.prop,0.01))-
      BvTLogPost(par = old_param,cor = rho.old,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(10,2))-
      log(dtruncnorm(b1.prop,a=0,b=1-old_param[j,2],mean=old_param[j,3],0.01))
    
    
    if(log(runif(1))<rho){
      old_param[j,3] = b1.prop
      acc.rate.b1[j] = acc.rate.b1[j]+1
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){b1.samp[(it-burn-1)/thin,j] = b1.prop}
    }else{
      if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){b1.samp[(it-burn-1)/thin,j] = old_param[j,3]}
    }
  }
  
  
  SS.t = crossprod(theta[-1,]-theta[-n,]) + S.t
  df.t = n+n.t
  W.old = riwish(df.t,SS.t)
  if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
    W.samp[[round((it-burn-1)/thin)]] = W.old
  }
  
  

  rho.prop = rho.old+runif(1,-0.1,0.1)
  while( rho.prop > 1 || rho.prop < -1){
    rho.prop = rho.old+runif(1,-0.1,0.1)
  }
  
  rho = BvTLogPost(par = old_param,cor = rho.prop,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(1,2))-
        BvTLogPost(par = old_param,cor = rho.old,W = W.old,Y = Y, theta = theta, nt = 10, psi = diag(1,2))

  if(log(runif(1)) < rho){
    rho.old = rho.prop
    if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
      rho.samp[(it-burn-1)/thin] = rho.prop
    }
  } else {
    if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
      rho.samp[(it-burn-1)/thin] = rho.old
    }
  }

  out = BvtKalFiltGarchC(Y=Y,GG=GG,FF=FF,W=W.old,m0=Y[1,],C0=C0,garch=old_param,cor=rho.old)
  theta = do.call(rbind,lapply(BvtKalffbsC(Y, out)$theta,t))
  
  if((it-burn-1)%%thin == 0 & (it-burn-1)/thin > 0 ){
    var1[[(it-burn-1)/thin]] = do.call(c,lapply(out$Q_t, function(x) x[1,1]))
    var2[[(it-burn-1)/thin]] = do.call(c,lapply(out$Q_t, function(x) x[2,2]))
    predicted[[(it-burn-1)/thin]] = do.call(rbind, lapply(out$y_t_1,t))
    sampled_theta[[(it-burn-1)/thin]] = theta
  }
  
    
}

