library(truncnorm)
Rcpp::sourceCpp("~/MvtKalFiltGarchC.cpp")
Rcpp::sourceCpp("~/FFBS.cpp")
Rcpp::sourceCpp("~/MvtLkhd.cpp")


## MCMC-specifications
burn = 20000
thin = 100
nsims = 5000
sims = (nsims)*thin+burn
sims


## MCMC-initialization
S.t = diag(10,p); n.t = 10;
S_t.samp = list()

Q_t.samp = list()

predicted = list()
sampled_theta = list()


garch_par_samp = list(matrix(c(0.05,0.2,0.5),byrow = T,nr=p,nc=3))
old_param = garch_par_samp[[1]]

W.samp = list()
W.old = diag(1,p)

cor.samp = list()
cor.old = diag(1,p)
V.old = tcrossprod(cor.old/apply(cor.old, 1,norm, "2"))

if(p==2){
  V.old = matrix(0.2)
}


## Kalman Filter-initialization
m0 = Y[1,]
C0 = diag(100,p)

## Random-walk tuning parameter
L1 = chol(diag(c(0.49,0.0025,0.0025)))  


## FFBS-step
out = MVTKalFiltGarch(Y=Y,GG=GG,FF=FF,W=W.old,m0=m0,C0=C0,garch=old_param,cor=V.old)
theta = do.call(rbind,lapply(ffbsC(Y, out)$theta,t))


## Obtain posterior samples
acc.rate = 0
for(it in 1:sims){
  if(it %% 500 == 0){print(it)}
  if(it %% 500 == 0){print(old_param)}
  if(it %% 500 == 0){print(W.old)}
  if(it %% 500 == 0){print(paste("acc.rate: ",acc.rate/it))}
  if(it %% 2000 == 0 & it > burn){
    par(mfrow = c(2,2))
    j = sample(1:4,1)
    print(paste("J=",j,sep = ""))
    plot(unlist(lapply(garch_par_samp, function(X) X[j,1])), type = "l",ylim = c(0,12))
    plot(unlist(lapply(garch_par_samp, function(X) X[j,2])), type = "l", ylim = c(0,0.5))
    plot(unlist(lapply(garch_par_samp, function(X) X[j,3])), type = "l",ylim = c(0.2,1))
  }
  
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
      acc.rate = acc.rate+1
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
    
    V.new = tcrossprod(cor.new)
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
    Q_t.samp[[(it-burn-1)/thin]] = out$Q_t
    predicted[[(it-burn-1)/thin]] = do.call(rbind, lapply(out$y_t_1,t))
    sampled_theta[[(it-burn-1)/thin]] = theta
  }
  
  
}


## Evaluate Posterior Samples

## a0 param
par(mfrow=c(1,1))
par(mar = c(3,5,2,2))
plot(unlist(lapply(garch_par_samp, function(X) X[1,1])), type = "l",ylab=expression(alpha[0]^{(1)}))
abline(h=model_par[1,1], col = "red")
hist(unlist(lapply(garch_par_samp, function(X) X[1,1])))
acf(unlist(lapply(garch_par_samp, function(X) X[1,1])), lag.max = 10)
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[1,1])))/1:(nsims-1), type = "l")
abline(h=model_par[1,1], col = "red")
quantile(unlist(lapply(garch_par_samp, function(X) X[1,1])),probs = c(0.025,0.975))


plot(unlist(lapply(garch_par_samp, function(X) X[2,1])), type = "l",ylab=expression(alpha[0]^{(2)}))
abline(h=model_par[2,1], col = "red")
hist(unlist(lapply(garch_par_samp, function(X) X[2,1])))
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[2,1])))/1:(nsims-1), type = "l")
quantile(unlist(lapply(garch_par_samp, function(X) X[2,1])),probs = c(0.025,0.975))
abline(h=model_par[2,1], col = "red")

plot(unlist(lapply(garch_par_samp, function(X) X[3,1])), type = "l",ylab=expression(alpha[0]^{(3)}))
abline(h=model_par[3,1], col = "red")
hist(unlist(lapply(garch_par_samp, function(X) X[3,1])))
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[3,1])))/1:(nsims-1), type = "l")
abline(h=model_par[3,1], col = "red")
quantile(unlist(lapply(garch_par_samp, function(X) X[3,1])),probs = c(0.025,0.975))


plot(unlist(lapply(garch_par_samp, function(X) X[4,1])), type = "l",ylab=expression(alpha[0]^{(4)}))
abline(h=model_par[4,1], col = "red")
hist(unlist(lapply(garch_par_samp, function(X) X[4,1])))
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[4,1])))/1:(nsims-1), type = "l")
quantile(unlist(lapply(garch_par_samp, function(X) X[4,1])),probs = c(0.025,0.975))










## a1 param
plot(unlist(lapply(garch_par_samp, function(X) X[1,2])), type = "l")
abline(h=model_par[1,2], col = "red")
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[1,2])))/1:(nsims-1), type = "l")
acf(unlist(lapply(garch_par_samp, function(X) X[1,2])), lag.max = 1000)
quantile(unlist(lapply(garch_par_samp, function(X) X[1,2])),probs = c(0.025,0.975))


plot(unlist(lapply(garch_par_samp, function(X) X[2,2])), type = "l")
abline(h=model_par[2,2], col = "red")
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[2,2])))/1:(nsims-1), type = "l")
quantile(unlist(lapply(garch_par_samp, function(X) X[2,2])),probs = c(0.025,0.975))


plot(unlist(lapply(garch_par_samp, function(X) X[3,2])), type = "l")
abline(h=model_par[3,2], col = "red")
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[3,2])))/1:(nsims-1), type = "l")
quantile(unlist(lapply(garch_par_samp, function(X) X[3,2])),probs = c(0.025,0.975))

plot(unlist(lapply(garch_par_samp, function(X) X[4,2])), type = "l")
abline(h=model_par[4,2], col = "red")
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[4,2])))/1:(nsims-1), type = "l")
quantile(unlist(lapply(garch_par_samp, function(X) X[4,2])),probs = c(0.025,0.975))



## b1 param
plot(unlist(lapply(garch_par_samp, function(X) X[1,3])), type = "l")
abline(h=model_par[1,3], col = "red")
acf(unlist(lapply(garch_par_samp, function(X) X[1,3])), lag.max = 1000)
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[1,3])))/1:(nsims-1), type = "l")
quantile(unlist(lapply(garch_par_samp, function(X) X[1,3])),probs = c(0.025,0.975))


plot(unlist(lapply(garch_par_samp, function(X) X[2,3])), type = "l")
abline(h=model_par[2,3], col = "red")
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[2,3])))/1:(nsims-1), type = "l")
quantile(unlist(lapply(garch_par_samp, function(X) X[2,3])),probs = c(0.025,0.975))


plot(unlist(lapply(garch_par_samp, function(X) X[3,3])), type = "l")
abline(h=model_par[3,3], col = "red")
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[3,3])))/1:(nsims-1), type = "l")
quantile(unlist(lapply(garch_par_samp, function(X) X[3,3])),probs = c(0.025,0.975))

plot(unlist(lapply(garch_par_samp, function(X) X[4,3])), type = "l")
abline(h=model_par[4,3], col = "red")
plot(cumsum(unlist(lapply(garch_par_samp, function(X) X[4,3])))/1:(nsims-1), type = "l")

quantile(unlist(lapply(garch_par_samp, function(X) X[4,3])),probs = c(0.025,0.975))




## Effective sample sizes
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[1,1]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[2,1]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[3,1]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[4,1]))[1:4999])

coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[1,2]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[2,2]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[3,2]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[4,2]))[1:4999])

coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[1,3]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[2,3]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[3,3]))[1:4999])
coda::effectiveSize(unlist(lapply(garch_par_samp, function(X) X[4,3]))[1:4999])


## Plotting theta_t
Time=1:n
theta.post = Reduce('+',sampled_theta)/length(sampled_theta)


par(mfrow=c(1,1))
par(mar = c(4,5,1,1))

# theta_t^1
LB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,1]))), 1, function(x){quantile(x,c(0.025))})[-1]
UB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,1]))), 1, function(x){quantile(x,c(0.975))})[-1]
plot(theta.post[,1], type = "l", ylab = expression(theta[1]),xlab = "", main = "", cex.lab = 2.5, cex.axis = 1.75,
     ylim = c(min(LB),max(UB)))
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(theta.post[-1,1])



# theta_t^2
LB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,2]))), 1, function(x){quantile(x,c(0.025))})[-1]
UB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,2]))), 1, function(x){quantile(x,c(0.975))})[-1]
plot(theta.post[,2], type = "l", ylab = expression(theta[2]),xlab = "", main = "", cex.lab = 2.5, cex.axis = 1.75,
     ylim = c(min(LB),max(UB)))
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(theta.post[-1,2])


# theta_t^3
LB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,3]))), 1, function(x){quantile(x,c(0.025))})[-1]
UB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,3]))), 1, function(x){quantile(x,c(0.975))})[-1]
plot(theta.post[,3], type = "l", ylab = expression(theta[3]),xlab = "", main = "", cex.lab = 2.5, cex.axis = 1.75,
     ylim = c(min(LB),max(UB)))
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(theta.post[-1,3])


# theta_t^4
LB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,4]))), 1, function(x){quantile(x,c(0.025))})[-1]
UB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,4]))), 1, function(x){quantile(x,c(0.975))})[-1]
plot(theta.post[,4], type = "l", ylab = expression(theta[4]),xlab = "", main = "", cex.lab = 2.5, cex.axis = 1.75,
     ylim = c(min(LB),max(UB)))
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(theta.post[-1,4])


## Plotting Y_t_1
forecast.post = Reduce('+',predicted)/length(predicted)


par(mfrow=c(1,1))
par(mar = c(3,5,2,2))
# Y_t^1
plot(Y[,1], type = "l", col = "grey", xlab = "", ylab = "HR", main = "", cex.lab = 2, cex.axis = 1.75)
lines(forecast.post[,1], col = "black", lwd = "2")

# Y_t^2
plot(Y[,2], type = "l", col="grey", xlab = "", ylab = "BP-Sys.", main = "", cex.lab = 2, cex.axis = 1.75)
lines(forecast.post[,2], lwd = "2")

# Y_t^3
plot(Y[,3], type = "l",col="grey", xlab = "", ylab = "BP-Dias.", main = "", cex.lab = 2, cex.axis = 1.75)
lines(forecast.post[,3], lwd = "2")

# Y_t^4
plot(Y[,4], type = "l", col="grey", xlab = "", ylab = "Resp.", main = "", cex.lab = 2, cex.axis = 1.75)
lines(forecast.post[,4], lwd = "2")



## Forecast time varying variance
Qt_post = list()
for(i in 1:length(Q_t.samp)){
  Qt_post[[i]] = do.call(rbind,lapply(Q_t.samp[[i]],function(x) diag(x)))
}
Qt_post = Reduce('+',Qt_post)/length(Q_t.samp)

## Plotting forecast time varying variance
plot(forecast.post[-1,1], type = "l", ylab = "Y1",xlab="", main = "", col="black", cex.lab = 2, cex.axis = 1.75,ylim = c(0.95*min(forecast.post[,1]),1.05*max(forecast.post[,1]) ) )
lines(forecast.post[-1,1]-1.96*sqrt(Qt_post[-1,1]))
lines(forecast.post[-1,1]+1.96*sqrt(Qt_post[-1,1]))


plot(forecast.post[-1,2], type = "l", ylab = "Y1",xlab="", main = "", col="black", cex.lab = 2, cex.axis = 1.75,ylim = c(0.5*min(forecast.post[,2]),1.2*max(forecast.post[,2]) ) )
lines(forecast.post[-1,2]-1.96*sqrt(Qt_post[-1,2]))
lines(forecast.post[-1,2]+1.96*sqrt(Qt_post[-1,2]))


plot(forecast.post[-1,3], type = "l", ylab = "Y1",xlab="", main = "", col="black", cex.lab = 2, cex.axis = 1.75,ylim = c(0.5*min(forecast.post[,3]),2*max(forecast.post[,3]) ) )
lines(forecast.post[-1,3]-1.96*sqrt(Qt_post[-1,3]))
lines(forecast.post[-1,3]+1.96*sqrt(Qt_post[-1,3]))


plot(forecast.post[-1,4], type = "l", ylab = "Y1",xlab="", main = "", col="black", cex.lab = 2, cex.axis = 1.75,ylim = c(0.5*min(forecast.post[,4]),2*max(forecast.post[,4]) ) )
lines(forecast.post[-1,4]-1.96*sqrt(Qt_post[-1,4]))
lines(forecast.post[-1,4]+1.96*sqrt(Qt_post[-1,4]))


## Plotting S_t
S_t.post = Reduce('+',S_t.samp)/length(S_t.samp)

# S1
par(mfrow=c(2,2))
par(mar = c(3,5,2,2))
UB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,1])))), 1, function(x){quantile(x,c(0.975))})[-1]
LB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,1])))), 1, function(x){quantile(x,c(0.025))})[-1]
plot(sqrt(S_t.post[-1,1]), type = "l", ylab = "Std. Dev",xlab="", main = "", col="black",ylim = c(min(LB),max(UB)), cex.lab = 2, cex.axis = 1.75)
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(sqrt(S_t.post[-1,1]))




# S2
UB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,2])))), 1, function(x){quantile(x,c(0.975))})[-1]
LB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,2])))), 1, function(x){quantile(x,c(0.025))})[-1]
plot(sqrt(S_t.post[-1,2]), type = "l", ylab = "Std. Dev", main = "",xlab="", col="black",ylim = c(min(LB),max(UB)), cex.lab = 2, cex.axis = 1.75)
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(sqrt(S_t.post[-1,2]))


# S3
UB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,3])))), 1, function(x){quantile(x,c(0.975))})[-1]
LB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,3])))), 1, function(x){quantile(x,c(0.025))})[-1]
plot(sqrt(S_t.post[-1,3]), type = "l", ylab = "Std. Dev",xlab="", main = "", col="black",ylim = c(min(LB),max(UB)), cex.lab = 2, cex.axis = 1.75)
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(sqrt(S_t.post[-1,3]))

# S4
UB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,4])))), 1, function(x){quantile(x,c(0.975))})[-1]
LB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,4])))), 1, function(x){quantile(x,c(0.025))})[-1]
plot(sqrt(S_t.post[-1,4]), type = "l", ylab = "Std. Dev",xlab="", main = "", col="black",ylim = c(min(LB),max(UB)), cex.lab = 2, cex.axis = 1.75)
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
# lines(sqrt(z.t[-1,4]),col="black", lwd = 1, lty = 1)
lines(sqrt(S_t.post[-1,4]))



## Estimates of correlation matrix
par(mfrow = c(1,1))
V.chol = gmat::uchol(V)
boxplot(do.call(rbind,lapply(cor.samp, function(x) c(x[1,1:4],x[2,2:4],x[3,3:4]))),
        main="Box Plot of the Correlation Components",
        names=expression(u[1][1],u[1][2],u[1][3],u[1][4],u[2][2],u[2][3],u[2][4],u[3][3],u[3][4]))
points(c(V.chol[1,1:4],V.chol[2,2:4],V.chol[3,3:4]), pch=4, col="red")







