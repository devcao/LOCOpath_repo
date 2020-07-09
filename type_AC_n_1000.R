
# load data

load('type_A_n_1000_p_50_v1.RData')
results_v1 = results
load('type_A_n_1000_p_50_v2.RData')
results_v2 = results
load('type_A_n_1000_p_50_v3.RData')
results_v3 = results

results=list()
results$Theta_list = append(append(results_v1$Theta_list, results_v2$Theta_list), results_v3$Theta_list)
results$glasso_results = append(append(results_v1$glasso_results, results_v2$glasso_results), results_v3$glasso_results)
results$LOCO_results = append(append(results_v1$LOCO_results, results_v2$LOCO_results), results_v3$LOCO_results)




# extract TPR and FPR
mTPR = c(); mFPR = c();
i = 1
for (epsilon in c(seq(0.0000001, 0.0001, length.out = 100), seq(0.0001,0.05,length.out = 100))){ 
  obj = LOCO_roc(results, epsilon = epsilon, iter=250)
  mTPR[i] = obj$TPR
  mFPR[i] = obj$FPR
  i = i + 1
}  

#mTPR_g = c(); mFPR_g = c();
#i = 1
#for (epsilon in seq(0.000001, 0.01, length.out = 50)){ 
#  obj = glasso_roc(results, epsilon = epsilon, iter = 250)
#  mTPR_g[i] = obj$TPR
#  mFPR_g[i] = obj$FPR
#  i = i + 1
#}  

# these from the GRASS paper
FPR_G = c(0.008,0.011,0.023,0.118,0.218,0.317,0.513)
TPR_G = 1-c(0.276, 0.233,0.184,0.118,0.093,0.076,0.05)

par(mfrow=c(1,2))

# draw ROC curve
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



load('type_C_n_1000_p_50_v1.RData')
results_v1 = results
load('type_C_n_1000_p_50_v2.RData')
results_v2 = results
load('type_C_n_1000_p_50_v3.RData')
results_v3 = results

results=list()
results$Theta_list = append(append(results_v1$Theta_list, results_v2$Theta_list), results_v3$Theta_list)
results$glasso_results = append(append(results_v1$glasso_results, results_v2$glasso_results), results_v3$glasso_results)
results$LOCO_results = append(append(results_v1$LOCO_results, results_v2$LOCO_results), results_v3$LOCO_results)



mTPR = c(); mFPR = c();
i = 1
for (epsilon in c(seq(0.0000001, 0.0001, length.out = 100), seq(0.0001,0.05,length.out = 100))){ 
  obj = LOCO_roc(results, epsilon = epsilon, iter=250)
  mTPR[i] = obj$TPR
  mFPR[i] = obj$FPR
  i = i + 1
}  

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
FPR_G = c(0.007,0.01,0.022,0.115,0.215,0.314,0.51)
TPR_G = 1-c(0.279, 0.237,0.186,0.118,0.092,0.074,0.049)

# draw ROC curve
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

