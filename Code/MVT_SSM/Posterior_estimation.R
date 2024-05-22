##################################################################
##################################################################
##############                               #####################
##############       A Gibbs Sampler         #####################
##############                               #####################
##################################################################
##################################################################
Rcpp::sourceCpp("./Functions/MvtKalFilt.cpp")
Rcpp::sourceCpp("./Functions/FFBS.cpp")
library(MCMCpack)
library(mvtnorm)


# We use Inverse-Wishart Priors for V and W
k=4; p=4

V.sample = list(diag(50,k))
W.sample = list(diag(0.1,p))
S.y = diag(10,k); n.y = 10;
S.t = diag(10,p); n.t = 10;

nsims = 5000
burn = 1500
thin = 2
est = seq(burn+1,nsims,thin)

m0 = Y[1,]
C0 = diag(100,p)

predicted = list()
sampled_theta = list()

for(i in 2:nsims){
  if(i %% 500 == 0){print(i)} 
  
  
  
  n = length(Y[,1])
  out = MVTKalFiltC(Y,GG,FF,V.sample[[i-1]],W.sample[[i-1]],m0,C0)
  theta = do.call(rbind, lapply(ffbsC(Y,out)$theta,t))
  
  if(round((i-burn-1)/thin)==((i-burn-1)/thin) & (i-burn-1) > 0 ){
    sampled_theta[[round((i-burn-1)/thin)]] = theta
    predicted[[round((i-burn-1)/thin)]] = do.call(rbind, lapply(out$y_t_1,t))
  }
  
  row.na = sort(which(is.na(Y),arr.ind = T)[,1])
  if(length(row.na) != 0){
    SS.y = crossprod((Y-t(FF%*%t(theta[-1,])))[-row.na,]) + S.y
  }else{
    SS.y = crossprod((Y-t(FF%*%t(theta[-1,])))) + S.y
  }
  
  
  # SS.y = crossprod((Y-theta[-1,])) + S.y
  df.y = n+n.y
  V.sample[[i]] = riwish(df.y,SS.y)
  
  SS.t = crossprod(theta[-1,]-t(GG%*%t(theta[-n,])) ) + S.t
  df.t = n+n.t
  W.sample[[i]] = riwish(df.t,SS.t)


}




## Posterior point estimates of the states

state1 = lapply(sampled_theta, function(x) x[,1]) #do.call(sum, lapply(theta, function(x) x[,1]))
state1.post = Reduce('+',state1)/length(est)

state2 = lapply(sampled_theta, function(x) x[,2]) #do.call(sum, lapply(theta, function(x) x[,1]))
state2.post = Reduce('+',state2)/length(est)

state3 = lapply(sampled_theta, function(x) x[,3]) #do.call(sum, lapply(theta, function(x) x[,1]))
state3.post = Reduce('+',state3)/length(est)

state4 = lapply(sampled_theta, function(x) x[,4]) #do.call(sum, lapply(theta, function(x) x[,1]))
state4.post = Reduce('+',state4)/length(est)




##################################################################
##################################################################
##############                               #####################
##############       WAIC for MVT mod        #####################
##############                               #####################
##################################################################
##################################################################
Rcpp::sourceCpp("Functions/waic1_std.cpp")

waic1 = MVT_waic1(V.sample, W.sample, Y,mode = "Std",FF = FF, GG=GG,m0=m0,C0=C0)
waic1

