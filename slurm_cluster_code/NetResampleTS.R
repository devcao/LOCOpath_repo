Net.Resample = function(X, Y, which.covariate, betaNull, multiTest, B = 500, parallel = FALSE, beta.init = 'adaptive', family = 'gaussian', beta.true = beta, ...){
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
  
  TS = ExactNet.TS(X = X, Y = Y, which.covariate = which.covariate, betaNull = betaNull, multiTest = multiTest, family = family, ...)
  
  if(beta.init == "adaptive"){
      bhat = adalasso(X = X, y = Y, k = 10, use.Gram = FALSE, both = TRUE, intercept = FALSE)$coefficients.adalasso
  }else if (beta.init == "Truth"){
    bhat = beta.true
  }else{
    stop("wrong initial beta")
  } 
  
  
  residual = Y - X%*%bhat
  
  #	}else{ # low dimenstion just use LSE
  #  	bhat = ginv(t(X)%*%X)%*%t(X)%*%Y    
  #  	residual = Y - X%*%bhat
  
  #	}
  
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




Net.Resample.Process = function(X, Y, multiTest, residual, b.Null, beta.index, betaNull, B = 500, parallel = FALSE, ...){
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


Net.Resample.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, 
                              iter = 500, B = 500, setting = 'dep',
                              which.covariate = 1, betaNull = 1, multiTest = FALSE, ...){
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


Path.Resample.Equality = function(X, Y, which.covariate.equal, betaNull, ...){
  n = nrow(X)
  p = ncol(X)
  X_adj = X
  X_adj[,2] = X[,1] + X[,2]
  return(Path.Resample(X=X_adj, Y=Y, which.covariate = 1, betaNull = betaNull, multiTest = FALSE, ...))
}



Path.Resample.Equality.Power = function(n = 100, p = 1000, beta=c(rep(1,10),rep(0,990)), rho=0.5, 
                              iter = 500, B = 500, setting = 'dep',
                              which.covariate.equal = 1, betaNull = 0, ...){
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
  

