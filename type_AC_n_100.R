# load data

load('type_A_n_100_p_50.RData')

# extract TPR and FPR
mTPR = c(); mFPR = c();
i = 1
for (epsilon in c(seq(0.0000001, 0.0001, length.out = 100), seq(0.0001,0.01,length.out = 50))){ 
  obj = LOCO_roc(results, epsilon = epsilon, iter=250)
  mTPR[i] = obj$TPR
  mFPR[i] = obj$FPR
  i = i + 1
}  



# these are from the GRASS paper
FPR_G = c(0.001,0.003,0.013,0.103,0.204,0.302,0.5011)
TPR_G = 1-c(0.512,0.449,0.353,0.23,0.186,0.148,0.094)

# draw ROC curve
par(mfrow=c(1,2))

plot(mFPR, mTPR, type = 'b', pch = 19, lty = 1,
     xlim=c(min(mFPR,FPR_G),max(mFPR,FPR_G)),
     ylim=c(min(mTPR,TPR_G),max(mTPR,TPR_G)),
     xlab = 'FPR', ylab = 'TPR'
)
#lines(mFPR_g, mTPR_g,  type = 'b', pch = 1, col = 'red')
lines(FPR_G, TPR_G,  type = 'b', pch = 1, col = 'red')

#abline(a=0, b=1, lty = 2)  
legend("bottomright", legend=c("LOCO path",'GRASS'), bty = 'n',
       col=c("black", "red"), pch = c(19,1), lty=1, cex=0.8)




load('type_C_n_100_p_50.RData')

mTPR = c(); mFPR = c();
i = 1
for (epsilon in c(seq(0.0000001, 0.0001, length.out = 100), seq(0.0001,0.01,length.out = 50))){ 
  obj = LOCO_roc(results, epsilon = epsilon, iter=250)
  mTPR[i] = obj$TPR
  mFPR[i] = obj$FPR
  i = i + 1
}  
# 
# mTPR_g = c(); mFPR_g = c();
# i = 1
# for (epsilon in seq(0.000001, 0.01, length.out = 50)){ 
#   obj = glasso_roc(results, epsilon = epsilon, iter = 250)
#   mTPR_g[i] = obj$TPR
#   mFPR_g[i] = obj$FPR
#   i = i + 1
# }  
#par(mfrow=c(1,2))


# these are from the GRASS paper
FPR_G = c(0.008,0.013,0.03,0.13,0.231,0.328,0.52)
TPR_G = 1-c(0.674, 0.592,0.485,0.322,0.255,0.207,0.135)

# draw ROC curve
plot(mFPR, mTPR, type = 'b', pch = 19, lty = 1,
     xlim=c(min(mFPR,FPR_G),max(mFPR,FPR_G)),
     ylim=c(min(mTPR,TPR_G),max(mTPR,TPR_G)),
     xlab = 'FPR', ylab = 'TPR'
)
#lines(mFPR_g, mTPR_g,  type = 'b', pch = 1, col = 'red')
lines(FPR_G, TPR_G,  type = 'b', pch = 1, col = 'red')

#abline(a=0, b=1, lty = 2)  
legend("bottomright", legend=c("LOCO path", 'GRASS'), bty = 'n',
       col=c("black", "red"), pch = c(19,1), lty=1, cex=0.8)

