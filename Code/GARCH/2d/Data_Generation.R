##################################################################
##################################################################
##############                               #####################
##############      Generating the data      #####################
##############                               #####################
##################################################################
##################################################################
rm(list = ls())
library(mvtnorm)
n = 1000

FF = diag(1,2)
GG = diag(1,2)
m0 = c(-50,10)
C0 = diag(10,2)
# V = diag(1,2)
# W = diag(1,2)

set.seed(1010101)
model_par = rbind(c(1,0.1,0.85,0.1),
                  c(0.5,0.4,0.1,0.1))
rho.pop = -0.2
e.t = rmvnorm(1,c(0,0),diag(1,2))
z1.t = c(model_par[1,1])
z2.t = c(model_par[2,1])

v1 = c(e.t[1])
v2 = c(e.t[2])
ST = list(diag(c(model_par[1,1],model_par[2,1] )))

for(i in 2:n){
  z1.t[i] = model_par[1,1] + model_par[1,2]*v1[i-1]^2 +model_par[1,3]*z1.t[i-1]
  z2.t[i] = model_par[2,1] +model_par[2,2]*v2[i-1]^2 + model_par[2,3]*z2.t[i-1]
  ST[[i]] = diag(sqrt(c(z1.t[i],z2.t[i]))) %*% matrix(c(1,rho.pop,rho.pop,1),nc = 2) %*% diag(sqrt(c(z1.t[i],z2.t[i])))
  v = t(chol(ST[[i]])) %*% rnorm(2,0,1)
  v1[i] = v[1]
  v2[i] = v[2]
  
}

par(mfrow = c(1,1))
plot(v1, type = "l")
lines(v2, type = "l", col = "red")

state = matrix(m0, nc = 2, nr = n+1)
Y = matrix(0,nr = n, nc = 2)
for(i in 1:n){
  state[i+1,] = t(GG %*% as.vector(state[i,]) +
             as.vector(c(rnorm(1,0,sqrt(model_par[1,4])),rnorm(1,0,sqrt(model_par[2,4]))) ))
  Y[i,] = FF %*% as.vector(state[i,]) +as.vector(c(v1[i],v2[i]))
  
}

par(mfrow=c(1,1))
plot(Y[,1], type = "l", ylim = c(min(Y),max(Y)))
lines(Y[,2], type = "l", col = "red")

