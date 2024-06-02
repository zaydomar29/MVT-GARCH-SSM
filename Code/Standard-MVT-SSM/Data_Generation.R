##################################################################
##################################################################
##############                               #####################
##############      Generating the data      #####################
##############                               #####################
##################################################################
##################################################################
rm(list = ls())

library(mvtnorm)
library(MCMCpack)
set.seed(201901)
n = 100
p = 2
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

