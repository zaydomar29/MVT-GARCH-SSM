rm(list = ls())
par(mfrow = c(1,1))


load("/Users/zaydomar/Dropbox/Zayd/Class_Notes/Research/DLM Paper/HR_BP.RData")


par.post = matrix(c(median(a0.samp[,1]),median(a1.samp[,1]),median(a1.samp[,1]),median(b1.samp[,1]),
                   median(a0.samp[,2]),median(a1.samp[,2]),median(a1.samp[,2]),median(b1.samp[,2])), nc = 4, byr=T)

W.post = Reduce('+',W.samp)/nsims


MVT_lgpostSS(par.post,mean(rho.samp),W.post,Y,cbind(state1.post,state2.post),n.t,S.t)
# -5809.134

## Truth
MVT_lgpostSS(model_par,mean(rho.samp),diag(c(0.3,0.3)),Y,cbind(state1.post,state2.post),n.t,S.t)
# -5734.977

## Non_Garch
MVT_lgpostSS(matrix(c(10.27,0,0,0,16.57,0,0,0), nr=2,byr = T),mean(rho.samp),diag(c(0.407,0.379)),Y,cbind(state1.post,state2.post),n.t,S.t)
#-5948.385




W11 = do.call(c, lapply(W.samp, function(x) x[1,1]))
W12 = do.call(c, lapply(W.samp, function(x) x[1,2]))
W22 = do.call(c, lapply(W.samp, function(x) x[2,2]))

round(c(median(a0.samp[,1])-1,median(a1.samp[,1])-0.1,median(b1.samp[,1])-0.8,median(W11)-0.3,median(rho.samp)+0.5),3)

round(c(median(a0.samp[,2])-2,median(a1.samp[,2])-0.4,median(b1.samp[,2])-0.5,median(W22)-0.3,median(W12)-0),3)


round(c(sd(a0.samp[,1]),sd(a1.samp[,1]),sd(b1.samp[,1]),sd(W11),sd(rho.samp)),3)

round(c(sd(a0.samp[,2]),sd(a1.samp[,2]),sd(b1.samp[,2]),sd(W22),sd(W12)),3)

