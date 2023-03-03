##################################################################
##################################################################
##############                               #####################
##############    Reading from the data      #####################
##############                               #####################
##################################################################
##################################################################
rm(list=ls())
setwd("/Users/zaydomar/Dropbox/Zayd/Class_Notes/BrainIT/")
physio = readRDS("Physiological.rds")


## Data 
# table(physio$Patient_Id)
Z1 = physio$HRT[physio$Patient_Id == 84016161]
plot(Z1, type = "l")
Z2 = physio$BPd[physio$Patient_Id == 84016161]
plot(Z2, type = "l")
Z3 = physio$RR[physio$Patient_Id == 84016161]
plot(Z3, type = "l")
Time = physio$Time_Stamp[physio$Patient_Id == 84016161]

# Z1 = Z1[!is.na(Z1)]
# Z2 = Z2[!is.na(Z2)]
HR = Z1[seq(1,length(Z1),1)]
BP = Z2[seq(1,length(Z2),1)]
Y = cbind(HR,BP)

par(mar=c(5,6,6,5))
par(mfrow = c(1,1))
plot(Time, HR, type = "l", main = "Heart Rate, Patient ID = 84016161", bty = "n", cex.lab = 3, cex.axis = 2.00, cex.main = 3)
plot(BP, type = "l", col = "black", main = "Blood Pressure, Patient ID = 84016161", bty = "n", cex.lab = 3, cex.axis = 2.00, cex.main = 2.25)

## DLM specification
FF = matrix(c(1), ncol = 1)
GG = matrix(c(1), nrow = 1, byrow = T)
m0 = c(mean(Y))
C0 = diag(1000,1)
n = length(Y)




##################################################################
##################################################################
##############                               #####################
##############  1)  Conjugate Bayesian       #####################
##############                               #####################
##################################################################
##################################################################
Rcpp::sourceCpp("/Users/zaydomar/Dropbox/Zayd/Class_Notes/Research TSA Books/DLM_Code/Random Walk/Non Garch/BayesConj.cpp")

nsims = 2e3
conj.out = BayesConjC(Y, nsims, c(0.5,0.5))

burn = 500
thin = 1
t = seq(burn, nsims, thin)

plot(conj.out$v[t], type = "l")
plot(conj.out$w1[t], type = "l")

plot(Y, type = "l")
lines(conj.out$PostMean, col = "red", lwd = 2)         ## Posterior Mean State

mean(conj.out$v[t])
mean(conj.out$w1[t])

acf(conj.out$v[t])
acf(conj.out$w1[t])


##################################################################
##################################################################
##############                               #####################
##############    2)  Garch Estimation       #####################
##############                               #####################
##################################################################
##################################################################

Rcpp::sourceCpp("/Users/zaydomar/Dropbox/Zayd/Class_Notes/Research TSA Books/DLM_Code/Random Walk/Garch/KalFiltGarch.cpp")

resid = Y-conj.out$PostMean[-1]
plot(resid, type = "l")            ## Estimated obs error

lkhdC = function(par){
  garch = exp(par[1:3])
  if(garch[2]+garch[3]>0.9999){return(Inf)}
  if(garch[1]<0 || garch[2]<0|| garch[3]<0){return(Inf)}
  W = matrix(exp(par[4]))
  out = KalFiltGarchC(Y,GG,FF,W,m0,C0,garch)
  f.t = unlist(out$y_t_1)
  Q.t = unlist(out$Q_t)
  -sum(dnorm(Y, mean = f.t, sd = sqrt(Q.t), log = T))
  
}

init = log(c(10, 0.1, 0.1,mean(conj.out$w1[t]) ))
lkhdC(init)
param.ll = list(lkhdC(init))
param = list()


for(i in 1:50){
  mle = optim(init, lkhdC, hessian = F)
  param[[i]] = exp(mle$par)
  param.ll[[i]] = mle$value
  init = log(c(2*i,i/(i+1),runif(1,0,1-i/(i+1)),mean(conj.out$w1[t])))
}


which(unlist(param.ll) == min(unlist(param.ll)))
order(unlist(param.ll))
param[[2]]
param[[18]]
param[[17]]
param[[79]]


par = do.call(rbind, param)
plot(unlist(param.ll),par[,1], type = "p")
plot(unlist(param.ll),par[,2], type = "p")
plot(unlist(param.ll),par[,3], type = "p")



##################################################################
##################################################################
##############                               #####################
##############     Innovations MCMC RCPP     #####################
##############                               #####################
##################################################################
##################################################################

Rcpp::sourceCpp("/Users/zaydomar/Dropbox/Zayd/Class_Notes/Research TSA Books/DLM_Code/Random Walk/Garch/Innovations.cpp")


nsims = 50000
burn = 1000
thin = 4
est = seq(from = burn, to = nsims, by = thin)
init_w = mean(conj.out$w1[t])
init_g = param[[18]][1:3]
L = chol(diag(c(1,0.1,0.1)))
t(L) %*% L


out = Innovations(Y, nsims, init_g, init_w, L)

a = out$a
w1 = out$w1

par(mfrow = c(2,2))
plot(a[est,1], type = "l", main = expression(paste("Trace Plot for ", a[0])), ylab =  expression(a[0]))
plot(a[est,2], type = "l", main = expression(paste("Trace Plot for ", a[1])), ylab =  expression(a[1]))
plot(a[est,3], type = "l", main = expression(paste("Trace Plot for ", b[1])), ylab =  expression(b[1]))
plot(w1[est], type = "l", main = expression(paste("Trace Plot for ", w)), ylab =  expression(w))

par(mfrow = c(3,2))
plot(a[est,1],a[est,2], xlab =  expression(a[0]), ylab =  expression(a[1]),
     main = expression(paste("Scatter Plot, ", a[0], " vs ", a[1])) )
plot(a[est,1],a[est,3], xlab =  expression(a[0]), ylab =  expression(b[1]),
     main = expression(paste("Scatter Plot, ", a[0], " vs ", b[1])) )
plot(a[est,2],a[est,3], xlab =  expression(a[1]), ylab =  expression(b[1]),
     main = expression(paste("Scatter Plot, ", a[1], " vs ", b[1])) )
abline(1,-1, col = "green")
plot(a[est,1],w1[est], xlab =  expression(a[0]), ylab =  expression(w),
     main = expression(paste("Scatter Plot, ", a[0], " vs ", w)) )
plot(a[est,2],w1[est], xlab =  expression(a[1]), ylab =  expression(w),
     main = expression(paste("Scatter Plot, ", a[1], " vs ", w)) )
plot(a[est,3],w1[est], xlab =  expression(b[1]), ylab =  expression(w),
     main = expression(paste("Scatter Plot, ", b[1], " vs ", w)) )

par(mfrow = c(2,2))
hist(a[est,1], breaks = 50, probability = T, main = expression(paste("Histogram, ", a[0])), 
     xlab = expression(a[0]) )
curve(dcauchy(x),from = 0, to = 10, col = "red", add = T)
legend(6,0.2,expression("Cauchy(0,1) Prior"),cex=.8,col=c("red"),lty=c(1), bty = "n")

hist(a[est,2], breaks = 50, probability = T, main = expression(paste("Histogram, ", a[1])), 
     xlab = expression(a[1]), xlim = c(0,1) )
curve(dcauchy(x),from = 0, to = 0.9, col = "red", add = T)
legend(0.1,25,expression("Cauchy(0,1) Prior"),cex=.8,col=c("red"),lty=c(1), bty = "n")

hist(a[est,3], breaks = 50, probability = T, main = expression(paste("Histogram, ", b[1])), 
     xlab = expression(b[1]) )
curve(dcauchy(x),from = 0, to = 1, col = "red", add = T)
legend(0.82,26,expression("Cauchy(0,1) Prior"),cex=.8,col=c("red"),lty=c(1), bty = "n")

hist(w1[est], breaks = 100, probability = T, main = expression(paste("Histogram, ", w)), 
     xlab = expression(w), xlim = c(0,2) )
curve(dgamma(x,1,1),from = 0, to = 1, col = "red", add = T)
legend(0.4,3,expression(paste(Gamma(1,1)," Prior")),cex=.8,col=c("red"),lty=c(1), bty = "n")


par(mfrow = c(1,1))
contour(w,v,z,nlevels = 100, xlab = "b1", ylab = "a1")
lines(a[est,3],a[est,2], col = "blue", type = "p", lwd = 1)
lines(mean(a[est,3]),mean(a[est,2]),0.1, col = "brown", type = "p", lwd = 5)
lines(mle.par[3],mle.par[2], col = "red", type = "p", lwd = 5)
lines(0.5,0.1, col = "green", type = "p", lwd = 5)



par(mfrow = c(2,2))
acf(a[est,1], main = expression(a[0]))
acf(a[est,2], main = expression(a[1]))
acf(a[est,3], main = expression(b[1]))
acf(w1[est], main = expression(w))


mean(a[est,1])
mean(a[est,2])
mean(a[est,3])

mean(a[est,2])+mean(a[est,3])
mean(a[est,1])/(1-(mean(a[est,2])+mean(a[est,3])))

mean(w1[est])

var(a[est,1])
var(a[est,2])
var(a[est,3])
var(w1[est])


median(a[est,1])
median(a[est,2])
median(a[est,3])
median(w1[est])




length(unique(a[est,1]))/length(est)
length(unique(w1[est]))/length(est)



##################################################################
##################################################################
##############                               #####################
##############      Generating the data      #####################
##############                               #####################
##################################################################
##################################################################

# Garch(1,1) errors
n = 1000
e.t = rnorm(n,0,1)
z.t = c(0)

a0 = median(a[est,1])
a1 = median(a[est,2])
b1 = median(a[est,3])
v = c(0)
for(i in 2:n){
  z.t[i] = a0 + a1*v[i-1]^2 + b1*z.t[i-1]
  v[i] = sqrt(z.t[i])*e.t[i]
}


par(mfrow = c(1,1))
plot(v, type = "l")
acf(v^2)

## DLM specification
FF = matrix(c(1), ncol = 1)
GG = matrix(c(1), nrow = 1, byrow = T)
theta = matrix(0, nrow = n, ncol = 1)
theta[1,] = matrix(m0, ncol = 1)
Y = c(theta[1,])
for(i in 2:n){
  theta[i,] = GG %*% theta[i-1, ] + as.vector(rnorm(1,mean = 0, sd = sqrt(median(w1[est]))))
  Y[i] = FF%*%theta[i,] + v[i]
}

plot(Y, type = "l")
lines(theta[,1], col = "green")







