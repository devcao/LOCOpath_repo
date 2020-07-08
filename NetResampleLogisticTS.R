


# This file contains the bootstrap code and power simulation code for LOCO statistic in Logistic and Poisson regression case

##### data generataion for logistic regression #########
binomDesign = function(n, p, beta, rho, intercept = 0.5){
  # data generataion for logistic regression
  # Args:
  # n,p,beta : sample size, number of features, regression oefficients (no intercept here)
  # rho: can be 'equl': compound symmetry with correlation = 0.8
  #             'weak_equl': compound symmetry with correlation = 0.5
  #             positive value:  toeplitz matrix with correlation = rho, the specified value
  #             0: independent case
  # intercept: intercept for logistic regression
  # Returns: A list of randomly genereated:
  #             X: design matrix
  #             Y: response vector
  #
  #
  if(rho == 'equl'){  # equi corr
    Sigma = matrix(rep(0.8,p*p),p,p)
    diag(Sigma) = rep(1,p)
    Mu=rep(0,p)
    X=rmvn(n, mu = Mu,sigma = Sigma)
    beta = c(intercept, beta)
    logit = cbind(1,X) %*% beta
    pr = exp(logit) / (1 + exp(logit))
    y = rbinom(n = n, size = 1, prob = pr)
  }else if (rho == "weak_equl"){
    Sigma = matrix(rep(0.5,p*p),p,p)
    diag(Sigma) = rep(1,p)
    Mu=rep(0,p)
    X=rmvn(n, mu = Mu, sigma = Sigma)
    beta = c(intercept, beta)
    logit = cbind(1,X) %*% beta
    pr = exp(logit) / (1 + exp(logit))
    y = rbinom(n = n, size = 1, prob = pr) 
  }else if (rho > 0){  # toeplitz matrix
    Sigma = toeplitz( rho^(0:(p - 1)) )
    Mu = rep(0,p)
    #X = rmvnorm(n,mean = Mu,sigma = Sigma)
    X = rmvn(n, mu = Mu,sigma = Sigma)
    beta = c(intercept, beta)
    logit = cbind(1,X) %*% beta
    pr = exp(logit) / (1 + exp(logit))
    y = rbinom(n = n, size = 1, prob = pr) 
  }else if (rho == 0){  # independent
    Sigma = diag(rep(1,p))
    Mu=rep(0,p)
    X = matrix(rnorm(n*p), n, p)
    beta = c(intercept, beta)
    logit = cbind(1,X) %*% beta
    pr = exp(logit) / (1 + exp(logit))
    y = rbinom(n = n, size = 1, prob = pr) 
  }
  return(list(X = X, Y = y))   
  
}

##### data generataion for Poisson regression #########
poissonDesign = function(n, p, beta, rho, intercept = 0.5){
  # data generataion for Poisson regression
  # Args:
  # n,p,beta : sample size, number of features, regression oefficients (no intercept here)
  # rho: can be 'equl': compound symmetry with correlation = 0.8
  #             'weak_equl': compound symmetry with correlation = 0.5
  #             positive value:  toeplitz matrix with correlation = rho, the specified value
  #             0: independent case
  # intercept: intercept for logistic regression
  # Returns: A list of randomly genereated:
  #             X: design matrix
  #             Y: response vector
  #
  #
  if(rho == 'equl'){  # equi corr
    Sigma = matrix(rep(0.8,p*p),p,p)
    diag(Sigma) = rep(1,p)
    Mu=rep(0,p)
    X=rmvn(n, mu = Mu,sigma = Sigma)
    beta = c(intercept, beta)
    log_mu = cbind(1,X) %*% beta
    mu = exp(log_mu) 
    y = rpois(n = n, lambda = mu)
  }else if (rho == "weak_equl"){
    Sigma = matrix(rep(0.5,p*p),p,p)
    diag(Sigma) = rep(1,p)
    Mu=rep(0,p)
    X=rmvn(n, mu = Mu, sigma = Sigma)
    beta = c(intercept, beta)
    log_mu = cbind(1,X) %*% beta
    mu = exp(log_mu) 
    y = rpois(n = n, lambda = mu)
  }else if (rho > 0){  # toeplitz matrix
    Sigma = toeplitz( rho^(0:(p - 1)) )
    Mu = rep(0,p)
    #X = rmvnorm(n,mean = Mu,sigma = Sigma)
    X = rmvn(n, mu = Mu,sigma = Sigma)
    beta = c(intercept, beta)
    log_mu = cbind(1,X) %*% beta
    mu = exp(log_mu) 
    y = rpois(n = n, lambda = mu)
  }else if (rho == 0){  # independent
    Sigma = diag(rep(1,p))
    Mu=rep(0,p)
    X = matrix(rnorm(n*p), n, p)
    beta = c(intercept, beta)
    log_mu = cbind(1,X) %*% beta
    mu = exp(log_mu) 
    y = rpois(n = n, lambda = mu)
  }
  return(list(X = X, Y = y))   
  
}




####### adpative LASSO/elastic net for GLM based on glmnet #################################
adaptive_glm = function(X, Y, alpha = 0.5, family = 'binomial'){
  # adpative LASSO/elastic net for GLM based on glmnet 
  # Args:
  # X, Y: design matrix and response vector
  # alpha: float between 0 and 1, argument in glmnet, controls the balance between L1 and L2 norm
  # family: can be ('binomial', 'poisson'). Other case didn't test.
  # Returns: vector, adaptive LASSO/elastic net estimator 
  #
  cv_lambda = cv.glmnet(X, Y, alpha = alpha, family=family)
  
  beta_lasso = as.numeric(coef(cv_lambda, s = cv_lambda$lambda.min))
  #print("beta_lasso")
  #print(beta_lasso[1:10])
  if(all(beta_lasso[-1]==0)){
    return(beta_lasso)
  }else{
  cv_adp_lambda = cv.glmnet(X, Y, penalty.factor = 1 / abs(beta_lasso[-1]), alpha = alpha, family=family)
  beta_adp_lasso = as.numeric( coef(cv_adp_lambda, s = cv_adp_lambda$lambda.min) )
  #print("beta_adp_lasso")
  #print(as.numeric(beta_adp_lasso)[1:10])
  }
  return(beta_adp_lasso)
}



# Logistic regression bootstrap the null distribution of Path-based statistic, return statistcal decesion
Net.Resample.Logistic = function(X, Y, which.covariate, betaNull, 
                                 multiTest, B = 500, parallel = FALSE, 
                                 beta.init.alpha = 1,
                                 beta.init = 'adaptive', beta.null.estimate = FALSE, beta.true = NULL, ...){
  # Logistic regression bootstrap the null distribution of Path-based statistic, return statistcal decesion
  #
  # Args:
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do simultanoues testing.
  #                   E.g, which.covariate = 1, we test for beta_1
  #                   E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  #                   E.g, if which.covariate = 1, we test for beta_1
  #                         then specify betaNull=0 
  #                   E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  #                         then specify betaNull = list(c(0,0)) 
  #                   The data type and length must match with which.covarite
  #                   Only support testing betaNull=0 or a list of 0
  #                   For testing non-zero betaNull, check code Net.Resample.Logistic.Con
  # multiTest: boolean, FALSE then do single test, TRUE then do simultanoues test
  #                     must match with which.covarite and betaNull
  #
  #	B: int, # of bootstrap replications
  #	parallel: boolean, run in parallel or nor. Note: This may nor be supported on Windows machine.
  # beta.init.alpha: specify the alpha for the adptive elastic net, default is 1, using adaptive LASSO
  # beta.init: ('adaptive', 'Truth')
  #             If 'adaptive', use adaptive LASSO. If 'Truth', use True value, and beta.true must be specified. 
  # beta.null.estimate: default FALSE. 
  #                     If true, the initial estimator will be estimated by removing beta_{which.covariate} first, 
  #                         then set the estimated beta_{which.covariate}=0.
  #                     If false,the initial estimator will be estimated first (all variables included), then set beta_{which.covariate} = 0
  #                     Simulations shows: if set FALSE: better power, if set TRUE: better size
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
  
  TS = ExactNet.TS.Logistic(X = X, Y = Y, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest, family = 'binomial', ...)
  
  if(beta.init == "adaptive"){
    if(beta.null.estimate){
      if(multiTest){stop('not implemented yet')}
      bhat = rep(0, p+1)
      bhat[-(which.covariate+1)]=adaptive_glm(X = X[, -which.covariate], Y = Y, alpha = beta.init.alpha)
     # cat('bhat:',bhat,'\n')
    }else{
      bhat=adaptive_glm(X = X, Y = Y, alpha = beta.init.alpha)
    }
  }else if (beta.init == "Truth"){
    bhat = beta.true
  }else if (beta.init == "glm"){
    bhat = coef( glm.fit(x=cbind(1,X), y=Y, family = binomial()) )
  }else{
    stop("wrong initial beta")
  } 
  
  #TS_null = matrix(NA, nrow = B, ncol = length(which.covariate))
  
  ################### HERE WE GO ! ! ! ###########################################
  
  
  ###################### This part could be parallelized ##################	
  count = 1
  
  for(wc_cov in which.covariate){   
    
    b.Null = bhat
    #b.Null[wc_cov] = 0
    
    if(multiTest) { 
      to.which.covariate = list(wc_cov)
      to.betaNull = list(betaNull[[count]])
      
      b.Null[wc_cov+1] = betaNull[[count]]
      
      
    }else{
      to.which.covariate = wc_cov
      to.betaNull = betaNull[count]
      
      b.Null[wc_cov+1] = betaNull[count]
      
    } # then run multiple testing
    
    #cat('bNull:',b.Null,'\n')
    
    TS_null = Net.Resample.Logistic.Process(X = X, Y = Y, multiTest = multiTest, b.Null = b.Null, betaNull = to.betaNull, 
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


############## A backend boostrapping function for Net.Resample.Logistic, won't be used alone #######################
Net.Resample.Logistic.Process = function(X, Y, multiTest, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  
  # A backend boostrapping function for Net.Resample.Logistic, won't be used alone
  # Args:
  # X, Y, multiTest: same argument in function Net.Resample.Logistic, check Net.Resample.Logistic for more details
  # residual: vector, residual
  # b.Null: vector, the working initial as in paper LCOC path 
  # beta.index: same as which.covariate, just different name
  # betaNull: same argument in function Net.Resample.Logistic, check Net.Resample.Logistic for more details
  # B: int, # of bootstrap replications
  # parallel: boolean, run in parallel or not. This may nor be supported on Windows machine.
  # Returns: null distribtion of LOCO path statistic under H0
  #
  #
  
  
  n = nrow(X)
  p = ncol(X)
  #cat('b.Null:',b.Null)
  TS_null = numeric()
  
  if(parallel){ # running in parallel
    mat = list()
    
    
    for(bs in 1:B){
    
      #Y = X %*% b.Null + boot_residual  ## change this part in logistic
      logit = cbind(1, X) %*% b.Null
      pr = exp(logit) / (1 + exp(logit))
      Y = rbinom(n = n, size = 1, prob = pr) 
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
    
    re_list = parLapply(cl, mat, ExactNet.TS.Logistic.Para, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull, family = 'binomial', ...)
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
      #Y = X %*% b.Null + boot_residual
      #cat('b.Null:', b.Null)
      logit = cbind(1, X) %*% b.Null
      pr = exp(logit) / (1 + exp(logit))
      Y = rbinom(n = n, size = 1, prob = pr) 
      #cat("ss:", sum(Y==0),'\n')
      
      TS_null[bs] = ExactNet.TS.Logistic(X = X, Y = Y, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull, family = 'binomial', ...)
      
    }
    
    return(TS_null)
  }
  
  
  
}

############# Poisson #################

###### Poisson regression bootstrap the null distribution of Path-based statistic, return statistcal decesion
Net.Resample.Poisson = function(X, Y, which.covariate, betaNull, multiTest, B = 500, parallel = FALSE, beta.init = 'adaptive', beta.init.alpha = 1, beta.true = beta, ...){
  # Poisson regression bootstrap the null distribution of Path-based statistic, return statistcal decesion
  #
  # Args:
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do simultanoues testing.
  #                   E.g, which.covariate = 1, we test for beta_1
  #                   E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  #                   E.g, if which.covariate = 1, we test for beta_1
  #                         then specify betaNull=0 
  #                   E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  #                         then specify betaNull = list(c(0,0)) 
  #                   The data type and length must match with which.covarite
  #                   Only support testing betaNull=0 or a list of 0
  #                   
  # multiTest: boolean, FALSE then do single test, TRUE then do simultanoues test
  #                     must match with which.covarite and betaNull
  #
  #	B: int, # of bootstrap replications
  #	parallel: boolean, run in parallel or nor. Note: This may nor be supported on Windows machine.
  # beta.init.alpha: specify the alpha for the adptive elastic net, default is 1, using adaptive LASSO
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
  
  TS = ExactNet.TS.Logistic(X = X, Y = Y, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest, family = 'poisson', ...)
  
  if(beta.init == "adaptive"){
    bhat=adaptive_glm(X = X, Y = Y, alpha = beta.init.alpha, family = 'poisson')
    cat('beta.init.alpha', beta.init.alpha, '\n', 'beta_hat', bhat[1:20], '\n')
  }else if (beta.init == "Truth"){
    bhat = beta.true
  }else if (beta.init == "glm"){
    bhat = coef( glm.fit(x=cbind(1,X), y=Y, family = binomial()) )
  }else{
    stop("wrong initial beta")
  } 
  
  #TS_null = matrix(NA, nrow = B, ncol = length(which.covariate))
  
  ################### HERE WE GO ! ! ! ###########################################
  
  
  ###################### This part could be parallelized ##################	
  count = 1
  
  for(wc_cov in which.covariate){   
    
    b.Null = bhat
    #b.Null[wc_cov] = 0
    
    if(multiTest) { 
      to.which.covariate = list(wc_cov)
      to.betaNull = list(betaNull[[count]])
      
      b.Null[wc_cov+1] = betaNull[[count]]
      
      
    }else{
      to.which.covariate = wc_cov
      to.betaNull = betaNull[count]
      
      b.Null[wc_cov+1] = betaNull[count]
      
    } # then run multiple testing
    
    
    TS_null = Net.Resample.Poisson.Process(X = X, Y = Y, multiTest = multiTest, residual = residual, b.Null = b.Null, betaNull = to.betaNull, 
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

############## A backend boostrapping function for Net.Resample.Poisson, won't be used alone #######################
Net.Resample.Poisson.Process = function(X, Y, multiTest, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  # A backend boostrapping function for Net.Resample.Poisson.Process, won't be used alone
  # Args:
  # X, Y, multiTest: same argument in function Net.Resample.Poisson, check Net.Resample.Poisson for more details
  # residual: vector, residual
  # b.Null: vector, the working initial as in paper LCOC path 
  # beta.index: same as which.covariate, just different name
  # betaNull: same argument in function Net.Resample.Poisson, check Net.Resample.Poisson for more details
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
      
      #Y = X %*% b.Null + boot_residual  ## change this part in logistic
      log_mu = cbind(1, X) %*% b.Null
      mu = exp(log_mu) 
      Y = rpois(n = n, lambda = mu) 
      
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
    
    re_list = parLapply(cl, mat, ExactNet.TS.Logistic.Para, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull, family = 'poisson', ...)
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
      #Y = X %*% b.Null + boot_residual
      logit = cbind(1, X) %*% b.Null
      pr = exp(logit) / (1 + exp(logit))
      Y = rbinom(n = n, size = 1, prob = pr) 
      
      TS_null[bs] = ExactNet.TS.Logistic(X = X, Y = Y, multiTest = multiTest, which.covariate = beta.index, betaNull = betaNull, family = 'poisson', ...)
      
    }
    
    return(TS_null)
  }
  
  
  
}




######## Simulate power of the LOCO path statistic for logistic regression  ###########
Net.Resample.Logistic.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, intercept=0.5,
                              iter = 500, B = 500, setting = 'dep',
                              which.covariate = 1, betaNull = 1, multiTest = FALSE, beta.init = 'Truth',...){
  
  # Simulate power of the LOCO path statistic for logistic regression 
  # Args:
  # 
  # rho: related to dependent design setting
  #	n,p,beta, intercept : sample size, number of features, regression oefficients (not include intercept), intercept
  # iter : int, # of iterations 
  # B: int, # of bootstrap replications
  # which.covariate, betaNull, multiTest: same argument in function Net.Resample.Logistic, check Net.Resample.Logistic for more details
  # setting: just use default 'dep', takes no effect in actual code
  # Returns: A list of:
  #             path.power: matrix, statistical decision in all the iterations 
  #             power: the simulated power 
  
  
  path.power = array(0,dim = c(iter,length(which.covariate),4))
    
  for(s in 1:iter){
    
    try({
    data = binomDesign(n = n, p = p, beta = beta, rho = rho, intercept = intercept)
    
    X_sp = data$X
    Y_sp = data$Y
    
    results = Net.Resample.Logistic(X = X_sp, Y = Y_sp, which.covariate = which.covariate, betaNull = betaNull, 
                                    multiTest = multiTest, B = B, beta.true = c(intercept, beta), beta.init = beta.init, ...)
    
    ######## keep track of MEMORY
    print("After Bootstrap:")
    print(mem_used())  
    ########
    
    path.power[s,,] = results$rej
    
    if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
    
    }) 
  }    
  
  
  power = apply(path.power,c(2,3),mean)
  
  return( list(path.power = path.power, power = power) )
  
  
}    


######## Simulate power of the LOCO path statistic for Poisson regression  ###########
Net.Resample.Poisson.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, intercept=0.5,
                                       iter = 500, B = 500, setting = 'dep',
                                       which.covariate = 1, betaNull = 1, multiTest = FALSE, beta.init = 'Truth',...){
  # Simulate power of the LOCO path statistic for Poisson regression 
  # Args:
  # 
  # rho: related to dependent design setting
  #	n,p,beta, intercept : sample size, number of features, regression oefficients (not include intercept), intercept
  # iter : int, # of iterations 
  # B: int, # of bootstrap replications
  # which.covariate, betaNull, multiTest: same argument in function Net.Resample.Poisson, check Net.Resample.Poisson for more details
  # setting: just use default 'dep', takes no effect in actual code
  # Returns: A list of:
  #             path.power: matrix, statistical decision in all the iterations 
  #             power: the simulated power 
  
  
  path.power = array(0,dim = c(iter,length(which.covariate),4))
  
  for(s in 1:iter){
    
    try({
      
    data = poissonDesign(n = n, p = p, beta = beta, rho = rho, intercept = intercept)
    
    X_sp = data$X
    Y_sp = data$Y
    
    results = Net.Resample.Poisson(X = X_sp, Y = Y_sp, which.covariate = which.covariate, betaNull = betaNull, 
                                    multiTest = multiTest, B = B, beta.true = c(intercept, beta), beta.init = beta.init, ...)
    
    ######## keep track of MEMORY
    print("After Bootstrap:")
    print(mem_used())  
    ########
    
    path.power[s,,] = results$rej
    
    if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
    
    })
  }    
  
  
  power = apply(path.power,c(2,3),mean)
  
  return( list(path.power = path.power, power = power) )
  
  
}    






################ de-sparsified logistic lasso power simulation #########################
desparse.Logistic.Power = function(n = 100, p = 1000, beta, rho, intercept = 0.5, iter = 500, which.covariate, betaNull){
  #Return the power of de-sparsified logistic lasso under different settings    
  # Args:
  #	setting: use 'dep',  currently	didn't add other options
  # rho: related to dependent design setting
  #	n,p,beta : sample size, number of features, regression coefficients
  # intercept: intercept of regression coefficients
  # iter : # of iterations 
  # which.covariate: Test which coefficients, if set to be 1, test beta_1=0 or not
  # betaNull: the null hypothesis of beta, default 0, testing beta_i=0
  # Return:
  #	A matrix of len(which.covariate) * 4 : Simulated power for j-th covariate
  #  
  
  
  if(length(betaNull) > 1){stop("now only support compute power of 1 coefficients")}
  
  proj.power = matrix(0,length(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){
    data = binomDesign(n = n, p = p, beta = beta, rho = rho, intercept = intercept)
    X_sp = data$X
    Y_sp = data$Y
    
    fit.proj <- lasso.proj(X_sp, Y_sp, standardize = TRUE, parallel = TRUE, family = 'binomial', ncores = 28)
    
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





##########################################################################################################################################


# Based on group LASSO Logistic regression, bootstrap the null distribution of Path-based statistic, return statistcal decesion
Net.Resample.Logistic.Multi = function(X, Y, which.covariate, betaNull, B = 500, parallel = FALSE, beta.init = 'adaptive', beta.true = NULL, ...){
  # Based on group LASSO Logistic regression, bootstrap the null distribution of Path-based statistic, return statistcal decesion
  #
  # Args:
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate: should be a list of vectors
  #                  E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  # betaNULL: ignroe or just use default, won't be used in the actual code
  #
  #	B: int, # of bootstrap replications
  #	parallel: boolean, run in parallel or nor. Note: This may nor be supported on Windows machine.
  # beta.init: ('adaptive', 'Truth')
  #             If 'adaptive', use adaptive LASSO. If 'Truth', use True value, and beta.true must be specified. 
  #
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
  
  TS = ExactNet.TS.Logistic.Multi(X = X, Y = Y, which.covariate = which.covariate, betaNull = betaNull, ...)
  
  if(beta.init == "adaptive"){
    bhat=adaptive_glm(X = X, Y = Y, alpha = 1)
  }else if (beta.init == "Truth"){
    bhat = beta.true
  }else if (beta.init == "glm"){
    bhat = coef( glm.fit(x=cbind(1,X), y=Y, family = binomial()) )
  }else{
    stop("wrong initial beta")
  } 
  
  #TS_null = matrix(NA, nrow = B, ncol = length(which.covariate))
  
  ################### HERE WE GO ! ! ! ###########################################
  
  
  ###################### This part could be parallelized ##################	
  count = 1
  
  for(wc_cov in which.covariate){   
    
    b.Null = bhat
    #b.Null[wc_cov] = 0
    
    #if(multiTest) { s
      to.which.covariate = list(wc_cov)
      to.betaNull = list(betaNull[[count]])
      
      b.Null[wc_cov+1] = betaNull[[count]]
      
      
    #}else{
      #to.which.covariate = wc_cov
      #to.betaNull = betaNull[count]
      
      #b.Null[wc_cov+1] = betaNull[count]
      
    #} # then run multiple testing
    
    
    TS_null = Net.Resample.Logistic.Multi.Process(X = X, Y = Y, residual = residual, b.Null = b.Null, betaNull = to.betaNull, 
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






############## A backend boostrapping function for Net.Resample.Logistic.Multi, won't be used alone #######################
Net.Resample.Logistic.Multi.Process = function(X, Y, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  # A backend boostrapping function for Net.Resample.Logistic.Multi, won't be used alone
  # Args:
  # X, Y, multiTest: same argument in function Net.Resample.Logistic.Multi, check Net.Resample.Logistic.Multi for more details
  # residual: vector, residual
  # b.Null: vector, the working initial as in paper LCOC path 
  # beta.index: same as which.covariate, just different name
  # betaNull: same argument in function Net.Resample.Logistic.Multi, check Net.Resample.Logistic.Multi for more details
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
      
      #Y = X %*% b.Null + boot_residual  ## change this part in logistic
      logit = cbind(1, X) %*% b.Null
      pr = exp(logit) / (1 + exp(logit))
      Y = rbinom(n = n, size = 1, prob = pr) 
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
    
    re_list = parLapply(cl, mat, ExactNet.TS.Logistic.Multi.Para, which.covariate = beta.index, betaNull = betaNull, ...)
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
      #Y = X %*% b.Null + boot_residual
      logit = cbind(1, X) %*% b.Null
      pr = exp(logit) / (1 + exp(logit))
      Y = rbinom(n = n, size = 1, prob = pr) 
      
      TS_null[bs] = ExactNet.TS.Logistic.Multi(X = X, Y = Y, which.covariate = beta.index, betaNull = betaNull, ...)
      
    }
    
    return(TS_null)
  }
  
  
  
}



######## Simulate power of the LOCO path statistic for logistic regression based on group lasso ###########
Net.Resample.Logistic.Multi.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, intercept=0.5,
                                             iter = 500, B = 500, 
                                             which.covariate, betaNull, beta.init = 'Truth',...){
  # Simulate power of the LOCO path statistic for logistic regression based on group lasso 
  # Args:
  #
  # rho: related to dependent design setting
  #	n,p,beta, intercept : sample size, number of features, regression oefficients (not include intercept), intercept
  # iter : int, # of iterations 
  # B: int, # of bootstrap replications
  # which.covariate, betaNull, multiTest: same argument in function Net.Resample.Logistic.Multi, check Net.Resample.Logistic.Multi for more details
  # Returns: A list of:
  #             path.power: matrix, statistical decision in all the iterations 
  #             power: the simulated power 
  
  
  path.power = array(0,dim = c(iter,length(which.covariate),4))
  
  for(s in 1:iter){
    
    try({
      data = binomDesign(n = n, p = p, beta = beta, rho = rho, intercept = intercept)
      
      X_sp = data$X
      Y_sp = data$Y
      
      results = Net.Resample.Logistic.Multi(X = X_sp, Y = Y_sp, which.covariate = which.covariate, betaNull = betaNull, 
                                            B = B, beta.true = c(intercept, beta), beta.init = beta.init, ...)
      
      ######## keep track of MEMORY
      print("After Bootstrap:")
      print(mem_used())  
      ########
      
      path.power[s,,] = results$rej
      
      if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
      
    }) 
  }    
  
  
  power = apply(path.power,c(2,3),mean)
  
  return( list(path.power = path.power, power = power) )
  
  
}    






#####################################################################
#####################################################################
#####################################################################


# Based my own coordinate descent code for logistic regression, bootstrap the null distribution of Path-based statistic, return statistcal decesion
# support non-zero betaNull
Net.Resample.Logistic.Con = function(X, Y, which.covariate, betaNull, B = 500, parallel = FALSE, beta.init = 'adaptive', beta.true = NULL, ...){
  # Based my own coordinate descent code for logistic regression, bootstrap the null distribution of Path-based statistic, return statistcal decesion
  #
  # Args:
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate: if is a vector, indicating which covariate we will be computing; 
  #                   E.g, which.covariate = 1, we test for beta_1
  #                  
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  #                   E.g, if which.covariate = 1, we test for beta_1
  #                         then specify betaNull=0 or any other value  
  #                   The data type and length must match with which.covarite
  #	B: int, # of bootstrap replications
  #	parallel: boolean, run in parallel or nor. Note: This may nor be supported on Windows machine.
  # beta.init: ('adaptive', 'Truth')
  #             If 'adaptive', use adaptive LASSO. If 'Truth', use True value, and beta.true must be specified. 
  #
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
  
  TS = ExactNet.TS.Logistic.Con(X = X, Y = Y, which.covariate = which.covariate, betaNull = betaNull, ...)
  
  if(beta.init == "adaptive"){
    bhat=adaptive_glm(X = X, Y = Y, alpha = 1)
  }else if (beta.init == "Truth"){
    bhat = beta.true
  }else if (beta.init == "glm"){
    bhat = coef( glm.fit(x=cbind(1,X), y=Y, family = binomial()) )
  }else{
    stop("wrong initial beta")
  } 
  
  #TS_null = matrix(NA, nrow = B, ncol = length(which.covariate))
  
  ################### HERE WE GO ! ! ! ###########################################
  
  
  ###################### This part could be parallelized ##################	
  count = 1
  
  for(wc_cov in which.covariate){   
    
    b.Null = bhat
    #b.Null[wc_cov] = 0
    
    #if(multiTest) { 
    #  to.which.covariate = list(wc_cov)
    #  to.betaNull = list(betaNull[[count]])
      
    #  b.Null[wc_cov+1] = betaNull[[count]]
      
      
    #}else{
      to.which.covariate = wc_cov
      to.betaNull = betaNull[count]
      
      b.Null[wc_cov+1] = betaNull[count]
      
    #} # then run multiple testing
    
    
    TS_null = Net.Resample.Logistic.Con.Process(X = X, Y = Y, residual = residual, b.Null = b.Null, betaNull = to.betaNull, 
                                            beta.index = to.which.covariate, B = B, parallel = parallel, ...)
    
    rej[count,1] = TS[count] > quantile(TS_null,0.8, na.rm = TRUE)
    rej[count,2] = TS[count] > quantile(TS_null,0.9, na.rm = TRUE)
    rej[count,3] = TS[count] > quantile(TS_null,0.95, na.rm = TRUE)
    rej[count,4] = TS[count] > quantile(TS_null,0.99, na.rm = TRUE)
    pval[count] = mean(TS_null >= TS[count], na.rm = TRUE)
    
    count = count + 1
    
  }
  
  ##########################################################
  
  return(list(rej = rej, pval = pval, TS_null = TS_null, TS = TS))
  
} 



############## A backend boostrapping function for Net.Resample.Logistic.Con, won't be used alone #######################
Net.Resample.Logistic.Con.Process = function(X, Y, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  # A backend boostrapping function for Net.Resample.Logistic.Con, won't be used alone
  # Args:
  # X, Y, multiTest: same argument in function Net.Resample.Logistic.Con, check Net.Resample.Logistic.Con for more details
  # residual: vector, residual
  # b.Null: vector, the working initial as in paper LCOC path 
  # beta.index: same as which.covariate, just different name
  # betaNull: same argument in function Net.Resample.Logistic.Con, check Net.Resample.Logistic.Con for more details
  # B: int, # of bootstrap replications
  # parallel: boolean, run in parallel or not. This may nor be supported on Windows machine.
  # Returns: null distribtion of LOCO path statistic under H0
  #
  
  n = nrow(X)
  p = ncol(X)
  
  TS_null = numeric()
  
  if(parallel){ # running in parallel
    mat = list()
    
    
    for(bs in 1:B){
      
      #Y = X %*% b.Null + boot_residual  ## change this part in logistic
      logit = cbind(1, X) %*% b.Null
      pr = exp(logit) / (1 + exp(logit))
      Y = rbinom(n = n, size = 1, prob = pr) 
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
    
    re_list = parLapply(cl, mat, ExactNet.TS.Logistic.Con.Para, which.covariate = beta.index, betaNull = betaNull, ...)
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
      #Y = X %*% b.Null + boot_residual
      logit = cbind(1, X) %*% b.Null
      pr = exp(logit) / (1 + exp(logit))
      Y = rbinom(n = n, size = 1, prob = pr) 
      
      TS_null[bs] = ExactNet.TS.Logistic.Con(X = X, Y = Y, which.covariate = beta.index, betaNull = betaNull, ...)
      
    }
    
    return(TS_null)
  }
  
  
  
}



######## Simulate power of the LOCO path statistic for logistic regression testing non-zero betaNull ############
Net.Resample.Logistic.Con.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, intercept=0.5,
                                             iter = 500, B = 500, 
                                             which.covariate, betaNull, beta.init = 'Truth',...){
  
  # Simulate power of the LOCO path statistic for logistic regression testing non-zero betaNull 
  # Args:
  #
  # rho: related to dependent design setting
  #	n,p,beta, intercept : sample size, number of features, regression oefficients (not include intercept), intercept
  # iter : int, # of iterations 
  # B: int, # of bootstrap replications
  # which.covariate, betaNull, multiTest: same argument in function Net.Resample.Logistic.Con, check Net.Resample.Logistic.Con for more details
  # Returns: A list of:
  #             path.power: matrix, statistical decision in all the iterations 
  #             power: the simulated power 
  
  
  #TS = matrix(NA, iter, len(which.covariate))
  #b_size_la <-  matrix(0,iter,4)
  
  path.power = array(0,dim = c(iter,length(which.covariate),4))
  
  for(s in 1:iter){
    
    try({
      data = binomDesign(n = n, p = p, beta = beta, rho = rho, intercept = intercept)
      
      X_sp = data$X
      Y_sp = data$Y
      
      results = Net.Resample.Logistic.Con(X = X_sp, Y = Y_sp, which.covariate = which.covariate, betaNull = betaNull, 
                                            B = B, beta.true = c(intercept, beta), beta.init = beta.init, ...)
      
      ######## keep track of MEMORY
      print("After Bootstrap:")
      print(mem_used())  
      ########
      
      path.power[s,,] = results$rej
      
      if(s %% 10 == 0){  cat("Now Computing:",s,"\n") }
      
    }) 
  }    
  
  
  power = apply(path.power,c(2,3),mean)
  
  return( list(path.power = path.power, power = power) )
  
  
}    



