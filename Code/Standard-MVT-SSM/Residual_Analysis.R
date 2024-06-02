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


hr_forecast = lapply(predicted, function(x) x[,1])
hr_forecast.post = Reduce('+',hr_forecast)/length(est)

bp_forecast = lapply(predicted, function(x) x[,2])
bp_forecast.post = Reduce('+',bp_forecast)/length(est)


y3_forecast = lapply(predicted, function(x) x[,3])
y3_forecast.post = Reduce('+',y3_forecast)/length(est)

y4_forecast = lapply(predicted, function(x) x[,4])
y4_forecast.post = Reduce('+',y4_forecast)/length(est)


