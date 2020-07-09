
binomDesign = function(n, p, beta, rho, intercept = 0.5){
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


poissonDesign = function(n, p, beta, rho, intercept = 0.5){
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





adaptive_glm = function(X, Y, alpha = 0.5, family = 'binomial'){
  
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




Net.Resample.Logistic = function(X, Y, which.covariate, betaNull, multiTest, B = 500, parallel = FALSE, 
                                 beta.init = 'adaptive', beta.init.alpha=1, beta.null.estimate = FALSE, beta.true = beta, ...){
  # Bootstrap the null distribution of Path-based statistic, and return reject or not
  #
  # Args:
  #	X, Y, which.covariate : feed in to Path-based TS function
  #	B : # of bootstrap replications
  #	parallel : run in parallel or nor
  #	exact : use exact TS or approx TS
  #	beta.true : for simulation only, the true value of beta would be used as initial estimates.
  #	
  # Return:  
  #	Reject or not under alpha = 0.2,0.1,0.05,0.01
  # 	p.values of the test
  
  
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
    
    #print(b.Null)
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



Net.Resample.Logistic.Process = function(X, Y, multiTest, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  # Bootstrap the null distribution of Path-based statistic of coef beta.index
  #
  # Args:
  #	X, Y, which.covariate : feed in to Path-based TS function
  #	B : # of bootstrap replications
  #	parallel : run in parallel or nor
  #	exact : use exact TS or approx TS
  #	beta.index : which coef 
  #
  # Return:  
  #	A vector of the bootstrapped null 
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


Net.Resample.Poisson = function(X, Y, which.covariate, betaNull, multiTest, B = 500, parallel = FALSE, beta.init = 'adaptive', beta.init.alpha = 1, beta.true = beta, ...){
  # Bootstrap the null distribution of Path-based statistic, and return reject or not
  #
  # Args:
  #	X, Y, which.covariate : feed in to Path-based TS function
  #	B : # of bootstrap replications
  #	parallel : run in parallel or nor
  #	exact : use exact TS or approx TS
  #	beta.true : for simulation only, the true value of beta would be used as initial estimates.
  #	
  # Return:  
  #	Reject or not under alpha = 0.2,0.1,0.05,0.01
  # 	p.values of the test
  
  
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



Net.Resample.Poisson.Process = function(X, Y, multiTest, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  # Bootstrap the null distribution of Path-based statistic of coef beta.index
  #
  # Args:
  #	X, Y, which.covariate : feed in to Path-based TS function
  #	B : # of bootstrap replications
  #	parallel : run in parallel or nor
  #	exact : use exact TS or approx TS
  #	beta.index : which coef 
  #
  # Return:  
  #	A vector of the bootstrapped null 
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






Net.Resample.Logistic.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, intercept=0.5,
                              iter = 500, B = 500, setting = 'dep',
                              which.covariate = 1, betaNull = 1, multiTest = FALSE, beta.init = 'Truth',...){
  #Return the power fo Path-based test under different settings    
  # Args:
  #	setting: different settings, check 'pathwise_simu_setting.R for details'	
  # 	rho: related to dependent design setting
  #	n,p,beta : sample size, features, coefficients
  # 	iter : # of iterations 
  #
  # Return:
  #	Simulated power
  #
  
  #TS = matrix(NA, iter, len(which.covariate))
  #b_size_la <-  matrix(0,iter,4)
  
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




Net.Resample.Poisson.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, intercept=0.5,
                                       iter = 500, B = 500, setting = 'dep',
                                       which.covariate = 1, betaNull = 1, multiTest = FALSE, beta.init = 'Truth',...){
  #Return the power fo Path-based test under different settings    
  # Args:
  #	setting: different settings, check 'pathwise_simu_setting.R for details'	
  # 	rho: related to dependent design setting
  #	n,p,beta : sample size, features, coefficients
  # 	iter : # of iterations 
  #
  # Return:
  #	Simulated power
  #
  
  #TS = matrix(NA, iter, len(which.covariate))
  #b_size_la <-  matrix(0,iter,4)
  
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







desparse.Logistic.Power = function(n = 100, p = 1000, beta, rho, intercept = 0.5, iter = 500, which.covariate, betaNull){
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
  
  if(length(betaNull) > 1){stop("now only support compute power of 1 coefficients")}
  
  proj.power = matrix(0,length(which.covariate),4)
  
  pval = matrix(NA, iter, p)
  
  for(s in 1:iter){
  try({
    data = binomDesign(n = n, p = p, beta = beta, rho = rho, intercept = intercept)
    X_sp = data$X
    Y_sp = data$Y
    
    fit.proj <- lasso.proj(X_sp, Y_sp, standardize = TRUE, parallel = TRUE, family = 'binomial', ncores = 28)
    
    pval[s,] = fit.proj$pval
    
    if(s %% 100 == 0){  cat("Now computing:", s, "\n")  }
    
  })
  }
  count = 1
  for(j in which.covariate){
    
    proj.power[count,1] = mean(pval[,j] < 0.2, na.rm=TRUE)
    proj.power[count,2] = mean(pval[,j] < 0.1, na.rm = TRUE)
    proj.power[count,3] = mean(pval[,j] < 0.05, na.rm = TRUE)
    proj.power[count,4] = mean(pval[,j] < 0.01, na.rm = TRUE)
    
    
    count = count + 1
  }
  
  
  return(list(path.power=pval, power = proj.power))
  
}  





#######################



Net.Resample.Logistic.Multi = function(X, Y, which.covariate, betaNull, B = 500, parallel = FALSE, beta.init = 'adaptive', beta.true = beta, ...){
  # Bootstrap the null distribution of Path-based statistic, and return reject or not
  #
  # Args:
  #	X, Y, which.covariate : feed in to Path-based TS function
  #	B : # of bootstrap replications
  #	parallel : run in parallel or nor
  #	exact : use exact TS or approx TS
  #	beta.true : for simulation only, the true value of beta would be used as initial estimates.
  #	
  # Return:  
  #	Reject or not under alpha = 0.2,0.1,0.05,0.01
  # 	p.values of the test
  
  
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



Net.Resample.Logistic.Multi.Process = function(X, Y, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  # Bootstrap the null distribution of Path-based statistic of coef beta.index
  #
  # Args:
  #	X, Y, which.covariate : feed in to Path-based TS function
  #	B : # of bootstrap replications
  #	parallel : run in parallel or nor
  #	exact : use exact TS or approx TS
  #	beta.index : which coef 
  #
  # Return:  
  #	A vector of the bootstrapped null 
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




Net.Resample.Logistic.Multi.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, intercept=0.5,
                                             iter = 500, B = 500, 
                                             which.covariate, betaNull, beta.init = 'Truth',...){
  #Return the power fo Path-based test under different settings    
  # Args:
  #	setting: different settings, check 'pathwise_simu_setting.R for details'	
  # 	rho: related to dependent design setting
  #	n,p,beta : sample size, features, coefficients
  # 	iter : # of iterations 
  #
  # Return:
  #	Simulated power
  #
  
  #TS = matrix(NA, iter, len(which.covariate))
  #b_size_la <-  matrix(0,iter,4)
  
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

Net.Resample.Logistic.Con = function(X, Y, which.covariate, betaNull, B = 500, parallel = FALSE, beta.init = 'adaptive', beta.true = beta, ...){
  # Bootstrap the null distribution of Path-based statistic, and return reject or not
  #
  # Args:
  #	X, Y, which.covariate : feed in to Path-based TS function
  #	B : # of bootstrap replications
  #	parallel : run in parallel or nor
  #	exact : use exact TS or approx TS
  #	beta.true : for simulation only, the true value of beta would be used as initial estimates.
  #	
  # Return:  
  #	Reject or not under alpha = 0.2,0.1,0.05,0.01
  # 	p.values of the test
  
  
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



Net.Resample.Logistic.Con.Process = function(X, Y, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
  # Bootstrap the null distribution of Path-based statistic of coef beta.index
  #
  # Args:
  #	X, Y, which.covariate : feed in to Path-based TS function
  #	B : # of bootstrap replications
  #	parallel : run in parallel or nor
  #	exact : use exact TS or approx TS
  #	beta.index : which coef 
  #
  # Return:  
  #	A vector of the bootstrapped null 
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






Net.Resample.Logistic.Con.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, intercept=0.5,
                                             iter = 500, B = 500, 
                                             which.covariate, betaNull, beta.init = 'Truth',...){
  #Return the power fo Path-based test under different settings    
  # Args:
  #	setting: different settings, check 'pathwise_simu_setting.R for details'	
  # 	rho: related to dependent design setting
  #	n,p,beta : sample size, features, coefficients
  # 	iter : # of iterations 
  #
  # Return:
  #	Simulated power
  #
  
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



