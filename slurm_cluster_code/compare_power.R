##################################################################################

# # Includes Power simulation of Projection method, etc.

##################################################################################




#######
##### loading required packages
#setwd('/hdi_simu')
require(LOCOpath)
source("path_related.R")
source("path_resample.R")
source("path_power.R")
len = length
# library(scalreg)
# source(file.path("~/SI/Sim-CI.R"))
# source(file.path("~/SI/ST.R"))

######


#ST <- function(X.f, Y.f, sub.size, test.set, M=500, alpha=c(0.2,0.1,0.05,0.01))

################ de-sparsified lasso #########################
ST.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate){
#Return the power of de-sparsified under different settings    
# Args:
# setting: different settings, check 'pathwise_simu_setting.R for details'  
#   rho: related to dependent design setting
# n,p,beta : sample size, features, coefficients
#   iter : # of iterations 
#
# Return:
# A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
#  

  ST.power.nst = matrix(0,iter,4)
  ST.power.st = matrix(0,iter,4)
  
  #pval = matrix(NA, iter, p)
  
  for(s in 1:iter){

    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)



    X_sp = data$X
    Y_sp = data$Y

    sub.size <- n*0.3
    
    results = ST(X.f = X_sp, Y.f = Y_sp, sub.size = sub.size, test.set = which.covariate)

    ST.power.nst[s,] <- results$nst
    ST.power.st[s,] <- results$st
    
    
    if(s %% 10 == 0){  cat("Now computing:", s, "\n")  }
    
  }
 
  simCI.power.nst=apply(ST.power.nst,2,mean)  
  
  simCI.power.st=apply(ST.power.st,2,mean)  
  
  return(list(nst = ST.power.nst, st = ST.power.st))

}  

##############################################################




######
################ de-sparsified lasso #########################
simCI.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate, betaNull){
#Return the power of de-sparsified under different settings    
# Args:
# setting: different settings, check 'pathwise_simu_setting.R for details'  
#   rho: related to dependent design setting
# n,p,beta : sample size, features, coefficients
#   iter : # of iterations 
#
# Return:
# A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
#  

  simCI.power.nst = matrix(0,iter,4)
  simCI.power.st = matrix(0,iter,4)
  
  #pval = matrix(NA, iter, p)
  
  for(s in 1:iter){

    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)



    X_sp = data$X
    Y_sp = data$Y
    
    results = Sim.CI(X = X_sp, Y = Y_sp, set = which.covariate, betaNull = betaNull)

    simCI.power.nst[s,] <- results$nst
    simCI.power.st[s,] <- results$st
    
    
    if(s %% 10 == 0){  cat("Now computing:", s, "\n")  }
    
  }
 
  simCI.power.nst=apply(simCI.power.nst,2,mean)  
  
  simCI.power.st=apply(simCI.power.st,2,mean)  
  
  return(list(nst = simCI.power.nst, st = simCI.power.st))

}  

##############################################################




################ de-sparsified lasso #########################
desparse.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate, betaNull){
#Return the power of de-sparsified under different settings    
# Args:
#	setting: different settings, check 'pathwise_simu_setting.R for details'	
# 	rho: related to dependent design setting
#	n,p,beta : sample size, features, coefficients
# 	iter : # of iterations 
#
# Return:
#	A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
#  
  
  if(len(betaNull) > 1){stop("now only support compute power of 1 coefficients")}

  proj.power = matrix(0,len(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){

    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)



   	X_sp = data$X
    Y_sp = data$Y
    
    fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
      
    pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }

  count = 1
  for(j in which.covariate){

    proj.power[count,1] = mean(pval[,j] < 0.2)
    proj.power[count,2] = mean(pval[,j] < 0.1)
    proj.power[count,3] = mean(pval[,j] < 0.05)
    proj.power[count,4] = mean(pval[,j] < 0.01)
  

    count = count + 1
  }
  
  
  return(proj.power)

}  

##############################################################


################ de-sparsified lasso #########################
normal.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate, betaNull){
#Return the power of de-sparsified under different settings    
# Args:
# setting: different settings, check 'pathwise_simu_setting.R for details'  
#   rho: related to dependent design setting
# n,p,beta : sample size, features, coefficients
#   iter : # of iterations 
#
# Return:
# A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
#  
  
  if(len(betaNull) > 1){stop("now only support compute power of 1 coefficients")}

  proj.power = matrix(0,len(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){

    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)



    X_sp = data$X
    Y_sp = data$Y
    

    pval[s,] = (summary(lm(y~., 
                         data = data.frame(y = data$Y, x = data$X)))$coefficients)[whichCov+1, 4] 

    #fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
    #pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }

  count = 1
  for(j in which.covariate){

    proj.power[count,1] = mean(pval[,j] < 0.2)
    proj.power[count,2] = mean(pval[,j] < 0.1)
    proj.power[count,3] = mean(pval[,j] < 0.05)
    proj.power[count,4] = mean(pval[,j] < 0.01)
  

    count = count + 1
  }
  
  
  return(proj.power)

}  



################ de-sparsified lasso #########################
TTest.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate, betaNull){
  #Return the power of de-sparsified under different settings    
  # Args:
  # setting: different settings, check 'pathwise_simu_setting.R for details'  
  #   rho: related to dependent design setting
  # n,p,beta : sample size, features, coefficients
  #   iter : # of iterations 
  #
  # Return:
  # A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  
  
  if(length(betaNull) > 1){stop("now only support compute power of 1 coefficients")}
  
  proj.power = matrix(0,length(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){
    
    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)
    
    
    
    X_sp = data$X
    Y_sp = data$Y
    
    
    pval[s,] = (summary(lm(y~., 
                           data = data.frame(y = Y_sp-betaNull*X_sp[,which.covariate], x = X_sp)))$coefficients)[which.covariate+1, 4] 
    
    #fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
    #pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }
  
  count = 1
  for(j in which.covariate){
    
    proj.power[count,1] = mean(pval[,j] < 0.2)
    proj.power[count,2] = mean(pval[,j] < 0.1)
    proj.power[count,3] = mean(pval[,j] < 0.05)
    proj.power[count,4] = mean(pval[,j] < 0.01)
    
    
    count = count + 1
  }
  
  
  return(proj.power)
  
}  
##########################################



################ de-sparsified lasso #########################
General.Test.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep'){
  #Return the power of de-sparsified under different settings    
  # Args:
  # setting: different settings, check 'pathwise_simu_setting.R for details'  
  #   rho: related to dependent design setting
  # n,p,beta : sample size, features, coefficients
  #   iter : # of iterations 
  #
  # Return:
  # A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  

  proj.power = c()#matrix(0,length(which.covariate),4)
  pvals = c()
  
  for(s in 1:iter){
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)

    X_sp = data$X
    Y_sp = data$Y
    
    lm_fit = lm(y~., data = data.frame(y = Y_sp, x = X_sp))
    pvals[s] = summary( glht(lm_fit, linfct = "x.1 - x.2 = 0") )$test$pvalues[1]
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }
  #print(pvals)
  proj.power[1] = mean(pvals < 0.2)
  proj.power[2] = mean(pvals < 0.1)
  proj.power[3] = mean(pvals < 0.05)
  proj.power[4] = mean(pvals < 0.01)
    
  return(proj.power)
  
}  
##########################################


##########################################
Logistic.Wald.Test.Power = function(n = 100, p = 80, beta, rho, intercept = 0.5, iter = 500, which.covariate, betaNull = 0){
  #Return the power of de-sparsified under different settings    
  # Args:
  # setting: different settings, check 'pathwise_simu_setting.R for details'  
  #   rho: related to dependent design setting
  # n,p,beta : sample size, features, coefficients
  #   iter : # of iterations 
  #
  # Return:
  # A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  
  
  if(length(betaNull) > 1){stop("now only support compute power of 1 coefficients")}
  proj.power = c()
  
  pval = c()
  
  for(s in 1:iter){
    
    data = binomDesign(n = n, p = p, beta = beta, rho = rho, intercept = intercept)
  
    X_sp = data$X
    Y_sp = data$Y
    if(betaNull==0){
    #if(FALSE){  
      glm_obj = glm(y~., data.frame(y = data$Y, x = data$X), family = binomial(), maxit = 5000)
      pval[s] = coef(summary(glm_obj))[2,4]
    }else{
      glm_obj = glm(y~., data.frame(y = data$Y, x = data$X), family = binomial(), maxit = 5000)
      pval[s] = 2*pnorm(-abs( (coef(summary(glm_obj))[which.covariate+1, 1] - betaNull)/coef(summary(glm_obj))[which.covariate+1, 2]) )
      
    }
    #pval[s,] = (summary(lm(y~., 
    #                       data = data.frame(y = Y_sp-betaNull*X_sp[,which.covariate], x = X_sp)))$coefficients)[which.covariate+1, 4] 
    
    #fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
    #pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }
  
  
    
    proj.power[1] = mean(pval < 0.2)
    proj.power[2] = mean(pval < 0.1)
    proj.power[3] = mean(pval < 0.05)
    proj.power[4] = mean(pval < 0.01)
    
    
  
  
  return(proj.power)
  
}  
##########################################



Poisson.Wald.Test.Power = function(n = 100, p = 80, beta, rho, intercept = 0.5, iter = 500, which.covariate, betaNull = 0){
  #Return the power of de-sparsified under different settings    
  # Args:
  # setting: different settings, check 'pathwise_simu_setting.R for details'  
  #   rho: related to dependent design setting
  # n,p,beta : sample size, features, coefficients
  #   iter : # of iterations 
  #
  # Return:
  # A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  
  
  if(length(betaNull) > 1){stop("now only support compute power of 1 coefficients")}
  proj.power = c()
  
  pval = c()
  
  for(s in 1:iter){
    
    data = poissonDesign(n = n, p = p, beta = beta, rho = rho, intercept = intercept)
    
    X_sp = data$X
    Y_sp = data$Y
    
    glm_obj = glm(y~., data.frame(y = data$Y, x = data$X), family = poisson(), maxit = 5000)
    pval[s] = coef(summary(glm_obj))[2,4]
    
    #pval[s,] = (summary(lm(y~., 
    #                       data = data.frame(y = Y_sp-betaNull*X_sp[,which.covariate], x = X_sp)))$coefficients)[which.covariate+1, 4] 
    
    #fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
    #pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }
  
  
  
  proj.power[1] = mean(pval < 0.2)
  proj.power[2] = mean(pval < 0.1)
  proj.power[3] = mean(pval < 0.05)
  proj.power[4] = mean(pval < 0.01)
  
  
  
  
  return(proj.power)
  
}  
##########################################





FTest.Power = function(n = 100, p = 1000, beta, rho, iter = 500, setting = 'dep', which.covariate = 1, betaNull = 1){
  #Return the power of de-sparsified under different settings    
  # Args:
  # setting: different settings, check 'pathwise_simu_setting.R for details'  
  #   rho: related to dependent design setting
  # n,p,beta : sample size, features, coefficients
  #   iter : # of iterations 
  #
  # Return:
  # A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  
  
  #### specified here
  betaNull = 1; which.covariate=1;
  ####
  
  if(length(betaNull) > 1){stop("now only support compute power of 1 coefficients")}
  
  proj.power = matrix(0,length(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){
    
    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)
    
    
    
    X_sp = data$X
    Y_sp = data$Y
    
    
    lm_model = lm(y~., data = data.frame(y=Y_sp-betaNull*X_sp[,which.covariate], x=X_sp))
    
    pval[s,] = linearHypothesis(lm_model, c('x.1=0', 'x.11=0', 'x.12=0'))[[6]][2]
    
    
    #pval[s,] = (summary(lm(y~., 
    #                       data = data.frame(y = Y_sp-betaNull*X_sp[,which.covariate], x = X_sp)))$coefficients)[which.covariate+1, 4] 
    #fit.proj <- lasso.proj(X_sp, Y_sp-betaNull*X_sp[,which.covariate], standardize = TRUE, parallel = TRUE, ncores = 40)
    #pval[s,] = fit.proj$pval
    
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  }
  
  count = 1
  for(j in which.covariate){
    
    proj.power[count,1] = mean(pval[,j] < 0.2)
    proj.power[count,2] = mean(pval[,j] < 0.1)
    proj.power[count,3] = mean(pval[,j] < 0.05)
    proj.power[count,4] = mean(pval[,j] < 0.01)
    
    
    count = count + 1
  }
  
  
  return(proj.power)
  
}  
# ###############
# require(car)
# 
# data = dataGen(setting = 'dep', n = 100, p = 80, beta = c(rep(1,10),rep(0,70)), rho = 0)
# X_sp = data$X; Y_sp = data$Y;
# a = lm(y~., data = data.frame(y = Y_sp-1*X_sp[,1], x = X_sp))
# 
# 
# b = linearHypothesis(a, c('x.1=0', 'x.11=0', 'x.12=0'))
# linearHypothesis(a, c('x.1=0', 'x.11=0', 'x.12=0'))[[6]][2]
# 
# 1-pf(112, df1=2,df2=87)
# 
# 









