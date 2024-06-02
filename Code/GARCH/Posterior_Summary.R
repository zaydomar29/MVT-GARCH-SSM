##################################################################
##################################################################
##############                               #####################
##############       Summary Posterior       #####################
##############                               #####################
##################################################################
##################################################################




## a0 param
par(mfrow=c(1,1))
plot(unlist(lapply(garch_par_samp, function(X) X[1,1])), type = "l")
abline(h=model_par[1,1], col = "red")

plot(unlist(lapply(garch_par_samp, function(X) X[2,1])), type = "l")
abline(h=model_par[2,1], col = "red")


## a1 param
plot(unlist(lapply(garch_par_samp, function(X) X[1,2])), type = "l")
abline(h=model_par[1,2], col = "red")

plot(unlist(lapply(garch_par_samp, function(X) X[2,2])), type = "l")
abline(h=model_par[2,2], col = "red")


## b1 param
plot(unlist(lapply(garch_par_samp, function(X) X[1,3])), type = "l")
abline(h=model_par[1,3], col = "red")

plot(unlist(lapply(garch_par_samp, function(X) X[2,3])), type = "l")
abline(h=model_par[2,3], col = "red")






## Plotting W
W.post = Reduce('+',W.samp)/length(W.samp)

plot(unlist(lapply(W.samp, function(X) X[1,1])), type = "l")
abline(h=model_par[1,4], col = "red")
acf(unlist(lapply(W.samp, function(X) X[1,1])), lag.max = 1000)

plot(unlist(lapply(W.samp, function(X) X[2,2])), type = "l")
abline(h=model_par[2,4], col = "red")


## Plotting theta_t
Time=1:n
theta.post = Reduce('+',sampled_theta)/length(sampled_theta)

plot(theta.post[,1], type = "l",ylab = expression(theta[1]),xlab = "", main = "", cex.lab = 2.5, cex.axis = 1.75)
LB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,1]))), 1, function(x){quantile(x,c(0.025))})[-1]
UB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,1]))), 1, function(x){quantile(x,c(0.975))})[-1]
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(theta.post[-1,1])
# lines(theta[,1], lty=2)

plot(theta.post[,2], type = "l", ylab = expression(theta[2]),xlab = "", main = "", cex.lab = 2.5, cex.axis = 1.75)
LB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,2]))), 1, function(x){quantile(x,c(0.025))})[-1]
UB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,2]))), 1, function(x){quantile(x,c(0.975))})[-1]
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(theta.post[-1,2])
# lines(theta[,2], lty=2)


plot(theta.post[,3], type = "l", ylab = expression(theta[3]),xlab = "", main = "", cex.lab = 2.5, cex.axis = 1.75)
LB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,3]))), 1, function(x){quantile(x,c(0.025))})[-1]
UB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,3]))), 1, function(x){quantile(x,c(0.975))})[-1]
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(theta.post[-1,3])
# lines(theta[,3], lty=2)


plot(theta.post[,4], type = "l", ylab = expression(theta[4]),xlab = "", main = "", cex.lab = 2.5, cex.axis = 1.75)
LB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,4]))), 1, function(x){quantile(x,c(0.025))})[-1]
UB = apply(t(do.call(rbind, lapply(sampled_theta, function(x) x[,4]))), 1, function(x){quantile(x,c(0.975))})[-1]
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
lines(theta.post[-1,4])
# lines(theta[,4], lty=2)




## Plotting Y_t_1
forecast.post = Reduce('+',predicted)/length(predicted)

plot(Y[,1], type = "l", col = "grey", xlab = "", ylab = "HR", main = "", cex.lab = 2, cex.axis = 1.75)
lines(forecast.post[,1], col = "black", lwd = "2")

plot(Y[,2], type = "l", col="grey", xlab = "", ylab = "BP-Sys.", main = "", cex.lab = 2, cex.axis = 1.75)
lines(forecast.post[,2], lwd = "2")


## Plotting S_t
S_t.post = Reduce('+',S_t.samp)/length(S_t.samp)

# S1
UB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,1])))), 1, function(x){quantile(x,c(0.975))})[-1]
LB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,1])))), 1, function(x){quantile(x,c(0.025))})[-1]

plot(sqrt(S_t.post[-1,1]), type = "l", ylab = "HR",xlab="", main = "", col="black",
     cex.lab = 2, cex.axis = 1.75)
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
# lines(sqrt(z.t[-1,1]),col="black", lwd = 0.75, lty = 1)
lines(sqrt(S_t.post[-1,1]))


# S2
UB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,2])))), 1, function(x){quantile(x,c(0.975))})[-1]
LB = apply(t(do.call(rbind, lapply(S_t.samp, function(x) sqrt(x[,2])))), 1, function(x){quantile(x,c(0.025))})[-1]

plot(sqrt(S_t.post[-1,2]), type = "l", ylab = "BP-Sys", main = "",xlab="", col="black", cex.lab = 2, cex.axis = 1.75)
polygon(c(Time,rev(Time)),c(UB,rev(LB)), col = "grey", border = NA)
# lines(sqrt(z.t[-1,2]),col="black", lwd = 0.75, lty = 1)
lines(sqrt(S_t.post[-1,2]))





# Parameter estimates
xtable::xtable(
  rbind(
    # Y1
    # model_par[1,1:3],
    round(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[1,])),2,median), 3),
    round(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[1,])),2,sd), 3),
    # Y2
    # model_par[2,1:3],
    round(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[2,])),2,median), 3),
    round(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[2,])),2,sd), 3),
    # Y3
    # model_par[3,1:3],
    round(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[3,])),2,median), 3),
    round(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[3,])),2,sd), 3),
    # Y4
    # model_par[4,1:3],
    round(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[4,])),2,median), 3),
    round(apply(do.call(rbind,lapply(garch_par_samp, function(x) x[4,])),2,sd), 3))   )


# Credible Intervals
round(quantile(do.call(rbind,lapply(garch_par_samp, function(x) x[1,]))[,1],c(0.0275,0.975)),2)
round(quantile(do.call(rbind,lapply(garch_par_samp, function(x) x[1,]))[,2],c(0.0275,0.975)),2)
round(quantile(do.call(rbind,lapply(garch_par_samp, function(x) x[1,]))[,3],c(0.0275,0.975)),2)

round(quantile(do.call(rbind,lapply(garch_par_samp, function(x) x[2,]))[,1],c(0.0275,0.975)),2)
round(quantile(do.call(rbind,lapply(garch_par_samp, function(x) x[2,]))[,2],c(0.0275,0.975)),2)
round(quantile(do.call(rbind,lapply(garch_par_samp, function(x) x[2,]))[,3],c(0.0275,0.975)),2)


# Rho
round(median(unlist(V.samp)),4)
round(sd(unlist(V.samp)),4)
round(quantile(unlist(V.samp),c(0.0275,0.975)),4)

# Estimates for W
# Y1
round(apply(do.call(rbind,lapply(W.samp, function(X) X[1,1:p])), M=2,FU=median),3)
round(apply(do.call(rbind,lapply(W.samp, function(X) X[1,1:p])), M=2,FU=quantile,c(0.0275,0.975)),3)
round(apply(do.call(rbind,lapply(W.samp, function(X) X[1,1:p])), M=2,FU=sd),3)

# Y2
round(apply(do.call(rbind,lapply(W.samp, function(X) X[2,2:p])), M=2,FU=median),3)
round(apply(do.call(rbind,lapply(W.samp, function(X) X[2,2:p])), M=2,FU=sd),3)
round(apply(do.call(rbind,lapply(W.samp, function(X) X[2,2:p])), M=2,FU=quantile,c(0.0275,0.975)),3)

# Y3
round(apply(do.call(rbind,lapply(W.samp, function(X) X[3,3:p])), M=2,FU=median),3)
round(apply(do.call(rbind,lapply(W.samp, function(X) X[3,3:p])), M=2,FU=sd),3)

# Y4
round(apply(do.call(rbind,lapply(W.samp, function(X) X[4,4:p])), M=2,FU=median),3)
round(apply(do.call(rbind,lapply(W.samp, function(X) X[4,4:p])), M=2,FU=sd),5)
## SD

# Correlation
plot(unlist(cor.samp))
abline(h=-0.2, col = "red")

# Estimates for V
par(mfrow = c(1,1))
V.chol = gmat::uchol(V)
