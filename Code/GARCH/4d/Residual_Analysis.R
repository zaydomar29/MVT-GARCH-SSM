
##################################################################
##################################################################
##############                               #####################
##############       Residual analysis       #####################
##############                               #####################
##################################################################
##################################################################


garch_par_post = rbind(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[1,])),2,median)[1:4999],
                       apply(do.call(rbind,lapply(garch_par_samp, function(x) x[2,])),2,median)[1:4999],
                       apply(do.call(rbind,lapply(garch_par_samp, function(x) x[3,])),2,median)[1:4999],
                       apply(do.call(rbind,lapply(garch_par_samp, function(x) x[4,])),2,median)[1:4999])

S_t.post = Reduce('+',S_t.samp[1:4999])/4999


## Posterior estimate of the correlation matrix
R.post = lapply(cor.samp, function(x){tcrossprod(x)})[1:4999]
R.post = Reduce('+',R.post)/length(R.post)


## Posterior estimate of the time varying observation variance matrix
V.post = list()
for(i in 1:n){
  V.post[[i]] = diag(sqrt(S_t.post[i,])) %*% R.post %*% diag(sqrt(S_t.post[i,]))
}

## Posterior estimates of the innovations
e.t.post = Y-theta.post[-1,]

# Normalized
par(mfrow = c(2,2))
par(mar = c(5,5,2,2))
qqnorm(e.t.post[,1], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.post[,1])

qqnorm(e.t.post[,2], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.post[,2])

qqnorm(e.t.post[,3], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.post[,3])

qqnorm(e.t.post[,4], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.post[,4])



# Hist Unnormalized
par(mfrow=c(2,2))
par(mar = c(3,5,2,2))
hist(e.t.post[,1], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

hist(e.t.post[,2], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

hist(e.t.post[,3], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

hist(e.t.post[,4], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")



par(mfrow = c(1,1))
pairs(~e.t.post[,1]+e.t.post[,2]+e.t.post[,3]+e.t.post[,4])


## Normalize the estimated innovations
e.t.norm = matrix(nrow=n,ncol=p)
for(i in 1:n){
  e.t.norm[i,] = t(solve(t(chol(V.post[[i]]))) %*% e.t.post[i,])
}


par(mfrow = c(2,2))
par(mar = c(5,5,2,2))
qqnorm(e.t.norm[,1], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.norm[,1])

qqnorm(e.t.norm[,2], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.norm[,2])

qqnorm(e.t.norm[,3], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.norm[,3])

qqnorm(e.t.norm[,4], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.norm[,4])






# Hist normalized
par(mfrow=c(2,2))
par(mar = c(3,5,2,2))
hist(e.t.norm[,1], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

hist(e.t.norm[,2], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

hist(e.t.norm[,3], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

hist(e.t.norm[,4], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")





# Scatter plots of the error terms
par(mfrow = c(1,1))
par(mar = c(5,5,2,2))
pairs(~e.t.norm[,1]+e.t.norm[,2]+e.t.norm[,3]+e.t.norm[,4], cex.lab = 2, cex.axis = 1.5,
      labels = c(expression(hat(e)[1],hat(e)[2],hat(e)[3],hat(e)[4])))
