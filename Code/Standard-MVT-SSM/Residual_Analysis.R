##################################################################
##################################################################
##############                               #####################
##############       Residual analysis       #####################
##############                               #####################
##################################################################
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


##################################################################
##################################################################
##############                               #####################
##############       WAIC for MVT mod        #####################
##############                               #####################
##################################################################
##################################################################
Rcpp::sourceCpp("./Functions/waic1_std.cpp")

est = seq(1501,5e3,2)
waic1 = MVT_waic1(V.sample[est], W.sample[est], Y,mode = "Std",FF = FF, GG=GG,m0=m0,C0=C0)
