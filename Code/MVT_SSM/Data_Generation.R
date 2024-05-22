## 4-d random walk multivariate normal SSM.
## changing FF and GG will give different types of SSMs.
## Here we have a simple independent multivariate SSM.


library(mvtnorm)
library(MCMCpack)
set.seed(201901)
n = 1000
p = 4
FF = diag(1,p)      # Observation matrix
GG = diag(1,p)      # State evolution matrix
m0 = c(-10,0,20,50) # starting point of the series
C0 = diag(10,p)

model_par = rbind(c(1,1,1,1),
                  c(0.1,0.1,0.1,0.1))  # observation level variance and state error variance
rho.obs = -0
rho.state = -0.0
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

