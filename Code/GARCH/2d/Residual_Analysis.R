## Posterior means of the estimated state vectors

state1 = lapply(sampled_theta, function(x) x[,1])
state1.post = Reduce('+',state1)/(nsims)

state2 = lapply(sampled_theta, function(x) x[,2])
state2.post = Reduce('+',state2)/(nsims)



## Posterior estimates of the innovations
theta.post = cbind(state1.post,state2.post)
e.t.post = Y-theta.post[-1,]

# Normalized
par(mfrow = c(2,1))
par(mar = c(5,5,2,2))
qqnorm(e.t.post[,1], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.post[,1])

qqnorm(e.t.post[,2], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.post[,2])


# Hist of Unnormalized
par(mfrow = c(2,1))
par(mar = c(3,5,2,2))
hist(e.t.post[,1], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

hist(e.t.post[,2], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

par(mfrow = c(1,1))
pairs(~e.t.post[,1]+e.t.post[,2])



## Normalize the estimated innovations
e.t.norm = matrix(nrow=n,ncol=p)
for(i in 1:n){
  e.t.norm[i,] = t(solve(t(chol(V.post[[i]]))) %*% e.t.post[i,])
}


par(mfrow = c(2,1))
par(mar = c(5,5,2,2))
qqnorm(e.t.norm[,1], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.norm[,1])

qqnorm(e.t.norm[,2], main = " ", cex.lab = 2, cex.axis = 1.5)
qqline(e.t.norm[,2])



pairs(~e.t.post[,1]+e.t.post[,2])
pairs(~e.t.norm[,1]+e.t.norm[,2])


plot(e.t.norm[,1])
plot(e.t.norm[,2])


cor(e.t.norm,use = "complete.obs")
cor(e.t.norm[,1],e.t.norm[,2],use = "complete.obs")

                
## Hist normalized
par(mfrow=c(2,1))
par(mar = c(3,5,2,2))
hist(e.t.norm[,1], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")

hist(e.t.norm[,2], freq = F, breaks = 15)
curve(dnorm(x,0,1),-4,4,add = T, col="red")



norm.err = e.t.norm[-which(is.na(e.t.norm),arr.ind = T)[,1], ]
table(is.na(norm.err))

## Test for constant correlation
test.res = rmgarch::DCCtest(norm.err)
test.res







##################################################################
##################################################################
##############                               #####################
##############       WAIC for MVT mod        #####################
##############                               #####################
##################################################################
##################################################################

Rcpp::sourceCpp("./Functions/waic1Garch.cpp")


waic.par = cbind(a0.samp[,1],a1.samp[,1],b1.samp[,1],
                 a0.samp[,2],a1.samp[,2],b1.samp[,2],
                 rho.samp)

waic1 = MVT_Garch_waic1(waic.par,W.samp, Y)
waic1


