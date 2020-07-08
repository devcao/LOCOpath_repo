
# This file contains the bootstrap code and power simulation code for LOCO statistic in linear regression case


##### Linear regression bootstrap the null distribution of LOCO statistic, return statistcal decesion #####
Net.Resample = function(X, Y, which.covariate, betaNull, multiTest, B = 500, parallel = FALSE, beta.init = 'adaptive', family = 'gaussian', beta.true = NULL, ...){
  # Linear regression bootstrap the null distribution of Path-based statistic, return statistcal decesion
  #
  # Args:
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do simultanoues testing.
  #                   E.g, which.covariate = 1, we test for beta_1
  #                   E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  #                   E.g, if which.covariate = 1, we test for beta_1
  #                         then specify betaNull=0 or any other value  
  #                   E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  #                         then specify betaNull = list(c(0,0)) or any other value.
  #                   The data type and length must match with which.covarite
  # multiTest: boolean, FALSE then do single test, TRUE then do simultanoues test
  #                     must match with which.covarite and betaNull
  #
  #	B: int, # of bootstrap replications
  #	parallel: boolean, run in parallel or nor. Note: This may nor be supported on Windows machine.
  # beta.init: ('adaptive', 'Truth')
  #             If 'adaptive', use adaptive LASSO. If 'Truth', use True value, and beta.true must be specified. 
  #	beta.true : vector, for simulation only, the true value of beta would be used as initial estimates. Default NULL. 
  #             If beta.init='Truth', beta.true must be specified. 
  #	
  # Return:  A list of 
  #	            rej: boolean, reject or not under alpha = 0.2,0.1,0.05,0.01
  # 	          pval: pvalue of the test
  #             TS_null: the boostrapped null distribution under H0
  #             TS: the test statistic based on the LOCO path
  
  
  n = nrow(X)
  p = ncol(X)
  rej = matrix(0,length(which.covariate),4)  # 1st : which cov, 2nd: na, 3rd: which alpha, 0.2,0.1,0.05,0.01
  pval = numeric()
  
  TS = ExactNet.TS(X = X, Y = Y, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest, family = family, ...)
  
  if(beta.init == "adaptive"){
      bhat = adalasso(X = X, y = Y, k = 10, use.Gram = FALSE, both = TRUE, intercept = FALSE)$coefficients.adalasso
  }else if (beta.init == "Truth"){
    bhat = beta.true
  }else{
    stop("wrong initial beta")
  } 
  
  
  residual = Y - X%*%bhat
  
 
  ################### HERE WE GO ! ! ! ###########################################
  ###################### This part could be parallelized ##################	
  count = 1
  
  for(wc_cov in which.covariate){   
    
    b.Null = bhat
    #b.Null[wc_cov] = 0
    
    if(multiTest) { 
      to.which.covariate = list(wc_cov)
      to.betaNull = list(betaNull[[count]])
      
      b.Null[wc_cov] = betaNull[[count]]
      
      
    }else{
      to.which.covariate = wc_cov
      to.betaNull = betaNull[count]
      
      b.Null[wc_cov] = betaNull[count]
      
    } # then run multiple testing
    
    
    TS_null = Net.Resample.Process(X = X, Y = Y, multiTest = multiTest, residual = residual, b.Null = b.Null, betaNull = to.betaNull, 
                                    beta.index = to.which.covariate, B = B, parallel = parallel, ...)
    
    rej[count,1] = TS[count] > quantile(TS_null,0.8)
    rej[count,2] = TS[count] > quantile(TS_null,0.9)
    rej[count,3] = TS[count] > quantile(TS_null,0.95)
    rej[count,4] = TS[count] > quantile(TS_null,0.99)
    pval[count] = mean(TS_null >= TS[count])
    
    count = count + 1
    
  }
  
  ##########################################################
  
  return(list(rej = rej, pval = pval, TS_null = TS_null, TS = TS))
  
} 



############## A backend boostrapping function, won't be used alone #######################
Net.Resample.Process = function(X, Y, multiTest, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  
  # A backend boostrapping function, won't be used alone
  # Args:
  # X, Y, multiTest: same argument in function Net.Resample, check Net.Resample for more details
  # residual: vector, residual
  # b.Null: vector, the working initial as in paper LCOC path 
  # beta.index: same as which.covariate, just different name
  # betaNull: same argument in function Net.Resample, check Net.Resample for more details
  # B: int, # of bootstrap replications
  # parallel: boolean, run in parallel or not. This may nor be supported on Windows machine.
  # Returns: null distribtion of LOCO path statistic under H0
  #
  #
  n = nrow(X)
  p = ncol(X)
  
  TS_null = numeric()
  
  if(parallel){ # running in parallel
    mat = list()
    
    
    for(bs in 1:B){
      ind = sample(1:n,replace = TRUE)
      boot_residual = residual[ind]
      #b_null = bhat
      #b_null[beta.index] = 0
      Y = X %*% b.Null + boot_residual  ## change this part in logistic
      mat[[bs]] = cbind(X,Y) 
    }
    
    #rgs = list(...)
    #Args = c(which.covariate = beta.index, betaNull = betaNull, exact = exact, multiTest = multiTest, args)
    
    # On a cluster, just use
    
    no_cores <- detectCores() 
    cat("n_cores detected:", no_cores, "\n")
    # Initiate cluster
    #cl <- makeCluster(no_cores)
    cl <- makeCluster(no_cores, type = "FORK")
    ###### if using window ###### 
    # load special packages
    #clusterEvalQ(cl, .libPaths("~/R"))
    #clusterEvalQ(cl, library(glmnet))
    #clusterEvalQ(cl, library(lars))
    #clusterEvalQ(cl, library(MASS))
    #clusterEvalQ(cl,library(pryr))
    #clusterEvalQ(cl,library(plus))
    #clusterEvalQ(cl,source("~/hdi_path/bin/pathwise_ts.R"))
    
    #clusterExport(cl, varlist = c("beta.index", "exact", "betaNull", "multiTest",...), envir = environment())
    #clusterExport(cl, varlist = 'Args', envir = environment())
    ###### if using window ######
    
    re_list = parLapply(cl, mat, ExactNet.TS.Para, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull, ...)
    #re_list = parLapply(cl, mat, Path.TS.Para, list = Args)
    
    ######## in case run out of MEMORY
    print("Cluster MEM:")
    print(mem_used())
    ########
    stopCluster(cl)  # END parallel bootsrap
    
    
    for(bss in 1:B){
      
      TS_null[bss] = re_list[[bss]]
      
    }
    
    return(TS_null)
    
  }else{ # not parallel, could be slow
    
    for(bs in 1:B){   
      ind = sample(1:n,replace = TRUE)
      boot_residual = residual[ind]
      #b_null = bhat 
      #b_null[beta.index] = 0
      Y = X %*% b.Null + boot_residual
      
      TS_null[bs] = ExactNet.TS(X = X, Y = Y, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull,...)
      
    }
    
    return(TS_null)
  }
  
  
  
}

######## Simulate power of the LOCO path statistic for linear regression ############
Net.Resample.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, 
                              iter = 500, B = 500, setting = 'dep',
                              which.covariate = 1, betaNull = 1, multiTest = FALSE, ...){
  # Simulate power of the LOCO path statistic for linear regression 
  # Args:
  #	setting: use 'dep',  currently	didn't add other options
  # rho: related to dependent design setting
  #	n,p,beta : sample size, number of features, regression oefficients
  # iter : int, # of iterations 
  # B: int, # of bootstrap replications
  # which.covariate, betaNull, multiTest: same argument in function Net.Resample, check Net.Resample for more details
  # Returns:
  #	A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  
  
  path.power = array(0,dim = c(iter,length(which.covariate),4))
  
  for(s in 1:iter){
    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)
    
    X_sp = data$X
    Y_sp = data$Y
    
    results = Net.Resample(X = X_sp, Y = Y_sp, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest, B = B, beta.true = beta, ...)
    
    ######## keep track of MEMORY
    print("After Bootstrap:")
    print(mem_used())  
    ########
    
    path.power[s,,] = results$rej
    
    if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
    
    
  }    
  
  
  power = apply(path.power,c(2,3),mean)
  
  return( list(path.power = path.power, power = power) )
  
  
}    

######## Bootstrap the null distribution of the LOCO path statistic for testing H0: beta1-beta2 = 0 in linear regression ############
Path.Resample.Equality = function(X, Y, which.covariate.equal=1, betaNull, ...){
  # Bootstrap the null distribution of the LOCO path statistic for testing H0: beta1-beta2 = 0 in linear regression
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate.equal: ignore, or use default, this is hardcoding the covariate to test after transformation 
  # betaNull: specifying the null hypothesis H0: beta1-beta2 = betaNull.
  #           if betaNull = 0, we test H0: beta1-beta2 = 0  
  n = nrow(X)
  p = ncol(X)
  X_adj = X
  X_adj[,2] = X[,1] + X[,2]
  
  return(Path.Resample(X=X_adj, Y=Y, which.covariate = 1, betaNull = betaNull, multiTest = FALSE, ...))
}


######## Simulate power of the LOCO path statistic for testing H0: beta1-beta2 = 0 in linear regression ############
Path.Resample.Equality.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, 
                              iter = 500, B = 500, setting = 'dep',
                              which.covariate.equal = 1, betaNull = 0, ...){
  # Simulate power of the LOCO path statistic for testing beta1=beta2 in linear regression 
  # Args:
  #	setting: use 'dep',  currently	didn't add other options
  # rho: related to dependent design setting
  #	n,p,beta : sample size, number of features, regression oefficients
  # iter : int, # of iterations 
  # B: int, # of bootstrap replications
  # which.covariate.equal: use default, this is hardcoding the covariate to test after transformation 
  # betaNull, specifying the null hypothesis H0: beta1-beta2 = betaNull.
  #           if betaNull = 0, we test H0: beta1-beta2 = 0  
  # Returns: A list of:
  #             path.power: matrix, statistical decision in all the iterations 
  #             power: the simulated power 
  #  
  
  path.power = matrix(NA, iter, 4)#array(0,dim = c(iter,length(which.covariate.equal),4))
  
  for(s in 1:iter){
    
    data = dataGen(setting = setting, n = n, p = p, beta = beta, rho = rho)
    
    X_sp = data$X
    Y_sp = data$Y
    
    results = Path.Resample.Equality(X = X_sp, Y = Y_sp, 
                                    which.covariate.equal = which.covariate.equal, betaNull = betaNull, 
                                    B = B, beta.true = beta, ...)
    
    ######## keep track of MEMORY
    print("After Bootstrap:")
    print(mem_used())  
    ########
    
    #path.power[s,,] = results$rej
    path.power[s,] = results$rej
    if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
    
    
  }    
  
  
  power = apply(path.power,2,mean)
  
  return( list(path.power = path.power, power = power) )
  
  
}    



#data = binomDesign(n = 100, p = 80, beta = c(rep(1,3), rep(0,77)), rho = 0)
#glm(y~., data.frame(y = data$Y, x = data$X), family = binomial)

#require('glmaag')
#glmaag(y = data$Y, x = data$X, fam = 'Logistic')
#data = depenDesign(n = 100, p = 1000, beta = c(0,rep(1,9),rep(0,990)), rho = 0)
#Path.Resample(X=data$X, Y=data$Y, which.covariate = 1, betaNull = 0, multiTest = FALSE, B = 500, parallel = FALSE, beta.init = 'adaptive', beta.true = 0)


#Net.Resample(X=data$X, Y=data$Y, which.covariate = 1, betaNull = 0, multiTest = FALSE, B = 500, parallel = FALSE, beta.init = 'adaptive', beta.true = 0)




#Net.Resample.Power(n = 100, p = 80, beta = c(0,1,1,rep(0,77)), rho = 0, iter = 500, B = 500, 
#                   setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE)
  

