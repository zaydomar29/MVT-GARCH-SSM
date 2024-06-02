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
p = 4
FF = diag(1,p)
GG = diag(1,p)
m0 = rep(0,p)
C0 = diag(10,p)


model_par = rbind(c(2,0.2,0.6,0.1),
                  c(2,0.1,0.8,0.1),
                  c(1,0.1,0.7,0.1),
                  c(1,0.2,0.5,0.1))



u1=c(runif(p,-1,1))
u2=c(0,runif(p-1,0,1))
u3=c(0,0,runif(p-2,0,1))
u4=c(0,0,0,runif(p-3,0,1))

u1 = u1/norm(u1,"2")
u2 = u2/norm(u2,"2")
u3 = u3/norm(u3,"2")
u4 = u4/norm(u4,"2")



U = rbind(u1,u2,u3,u4)
V = tcrossprod(U)

v = matrix(rmvnorm(1,rep(0,p),diag(1,p)), nc=p, nr = n, byr = T)
z.t = matrix(model_par[,1], nrow = n, nc = p,byr = T)
ST = list(diag(c(model_par[1:p,1])))


for(i in 2:n){
  z.t[i,] = model_par[1:p,1] + model_par[1:p,2]*v[i-1,]^2+model_par[1:p,3]*z.t[i-1,]
  ST[[i]] = diag(sqrt(z.t[i-1,])) %*% V %*% diag(sqrt(z.t[i-1,]))
  v[i,] = t(chol(ST[[i]])) %*% rnorm(p,0,1)
  
}

par(mfrow = c(1,1))
plot(z.t[,1], type="l", ylim=c(min(z.t),max(z.t)))
lines(z.t[,2], type="l", col = 2)
lines(z.t[,3], type="l", col = 3)
lines(z.t[,4], type="l", col = 4)


cor(v)-V

state = matrix(c(10.6,-50.5,0,-10), nc = p, nr = n+1,byrow = T)
Y = matrix(0,nr=n,nc=p)
for(i in 1:n){
  state[i+1,] = t(GG %*% as.vector(state[i,]) +
                    as.vector(chol(diag(model_par[1:p,4]))%*%rnorm(p)))
  Y[i,] = FF %*% as.vector(state[i,]) +v[i,]
}

plot(Y[,1], type = "l", ylim = c(min(Y,na.rm = T),max(Y,na.rm = T)))
lines(Y[,2], type = "l", col=2)
lines(Y[,3], type = "l", col=3)
lines(Y[,4], type = "l", col=4)
