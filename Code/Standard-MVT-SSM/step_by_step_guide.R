################################
################################
###                          ###
###  STEP-1 DATA GENERATION  ###
###                          ###
################################
################################

rm(list = ls())

library(mvtnorm)
library(MCMCpack)
set.seed(202406)
n = 1000
p = 4
FF = diag(1,p)
GG = diag(1,p)
m0 = c(-10,0,20,50) #rep(50,p)
C0 = diag(10,p)

model_par = rbind(c(1,1,1,1),
                  c(0.1,0.1,0.1,0.1))

sigma.v = diag(model_par[1,])
sigma.w = diag(model_par[2,])

theta = matrix(m0, nc = p, nr = n+1)
Y = matrix(0,nr = n, nc = p)
for(i in 1:n){
  theta[i+1,] = t(GG %*% as.vector(theta[i,])) + rmvnorm(1,rep(0,p),sigma.w)
  Y[i,] = t(FF %*% as.vector(theta[i,])) + rmvnorm(1,rep(0,p),sigma.v)
}

par(mfrow=c(1,1))
plot(Y[,1], type = "l", ylim=c(min(Y),max(Y)))
lines(Y[,2], type = "l", col="red")
lines(Y[,3], type = "l", col="green")
lines(Y[,4], type = "l", col="blue")

par(mfrow=c(1,1))
plot(theta[,1], type = "l")
plot(theta[,2], type = "l", col="red")
plot(theta[,3], type = "l", col="green")
plot(theta[,4], type = "l", col="blue")


#################################################
#################################################
###                                           ###
###  STEP-2 BAYESIAN ESTIMATION OF THE MODEL  ###
###                                           ###
#################################################
#################################################

Rcpp::sourceCpp("./Functions/MvtKalFilt.cpp")
Rcpp::sourceCpp("./Functions/FFBS.cpp")
library(MCMCpack)
library(mvtnorm)

#################################################################
##                 Gibbs sampling and M-H Step                 ##
#################################################################

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
  if(i %% 100 == 0){print(i)} 
  
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
  
  
  df.y = n+n.y
  V.sample[[i]] = riwish(df.y,SS.y)
  
  SS.t = crossprod(theta[-1,]-t(GG%*%t(theta[-n,])) ) + S.t
  df.t = n+n.t
  W.sample[[i]] = riwish(df.t,SS.t)


}

##################################
##################################
###                            ###
###  STEP-3 POSTERIOR SUMMARY  ###
###                            ###
##################################
##################################

##################################################################
##                       GARCH parameters                       ##
##################################################################



state1 = lapply(sampled_theta, function(x) x[,1]) #do.call(sum, lapply(theta, function(x) x[,1]))
state1.post = Reduce('+',state1)/length(est)

state2 = lapply(sampled_theta, function(x) x[,2]) #do.call(sum, lapply(theta, function(x) x[,1]))
state2.post = Reduce('+',state2)/length(est)

state3 = lapply(sampled_theta, function(x) x[,3]) #do.call(sum, lapply(theta, function(x) x[,1]))
state3.post = Reduce('+',state3)/length(est)

state4 = lapply(sampled_theta, function(x) x[,4]) #do.call(sum, lapply(theta, function(x) x[,1]))
state4.post = Reduce('+',state4)/length(est)


y1_forecast = lapply(predicted, function(x) x[,1])
y1_forecast.post = Reduce('+',y1_forecast)/length(est)

y2_forecast = lapply(predicted, function(x) x[,2])
y2_forecast.post = Reduce('+',y2_forecast)/length(est)


y3_forecast = lapply(predicted, function(x) x[,3])
y3_forecast.post = Reduce('+',y3_forecast)/length(est)

y4_forecast = lapply(predicted, function(x) x[,4])
y4_forecast.post = Reduce('+',y4_forecast)/length(est)


#######################################################
#######################################################
###                                                 ###
###  STEP-4 MODEL EVALUATION AND RESIDUAL ANALYSIS  ###
###                                                 ###
#######################################################
#######################################################


##################################################################
##                        WAIC for SSM                          ##
##################################################################

Rcpp::sourceCpp("./Functions/waic1_std.cpp")

est = seq(1501,5e3,2)
waic1 = MVT_waic1(V.sample[est], W.sample[est], Y,mode = "Std",FF = FF, GG=GG,m0=m0,C0=C0)


#################################################################
##                      Residual analysis                      ##
#################################################################

## Residual estimates
e.post = Y-cbind(state1.post,state2.post,state3.post,state4.post)[-1,]
pairs(~e.post[,1]+e.post[,2]+e.post[,3]+e.post[,4])

