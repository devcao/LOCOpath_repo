require(SIS)
require(LOCOpath)
require(mvnfast)
require(SIS)
source("path_power.R")
source('NetTS.R')
source('NetResampleTS.R')
source('NetResampleLogisticTS.R')


TS_util_fun = function(x_sp, y_sp, which.covariate = 1, betaNull = 0, multiTest = FALSE, path.method = "lars",
                       norm = "L1", normalize = TRUE, intercept = FALSE){
  return(
    ExactPath.TS(X = x_sp, Y = y_sp, which.covariate = which.covariate, betaNull = betaNull,
                 multiTest = multiTest, path.method = path.method,
                 norm = norm, normalize = normalize, intercept = intercept)
  )
}


logistic_TS_util_fun = function(x_sp, y_sp, which.covariate = 1, betaNull = 0, multiTest = FALSE, 
                                nlambda = 100,  alpha = 1, family = "binomial", lambda.min.ratio = 0.001,
                                norm = 'L2.squared'){
  return(
    ExactNet.TS.Logistic(X = x_sp, Y = y_sp, which.covariate = which.covariate, betaNull = betaNull,
                 multiTest = multiTest, nlambda = nlambda, alpha = alpha, family = family, 
                 lambda.min.ratio = lambda.min.ratio, norm = norm)
  )
}





logistic_screen_simu = function(n, p, signal, rho, iter=200, norm = 'L1'){
  
result_L1 = c(); result_L2 = c()
for (i in 1:iter){
  
  #n = 20;p = 100;rho = 0.5
  data = binomDesign(n = n, p = p, beta = c(rep(signal,3), rep(0, p-3)), rho, intercept = 0.5)

  n_threads = detectCores()
  cl = makeCluster(n_threads, type = "FORK")
  
  TS=unlist(parLapply(cl, X=1:p, logistic_TS_util_fun, 
                    x_sp = data$X, y_sp = data$Y,
                    betaNull = 0, multiTest = FALSE, nlambda = 100, 
                    alpha = 1,family = "binomial", lambda.min.ratio = 0.001,
                    norm = 'L1'))
  stopCluster(cl)
  
  result_L1[i] = all( TS[1:3] %in% sort(TS, decreasing = TRUE)[1:n-1] )
    
  cl = makeCluster(n_threads, type = "FORK") 
  
  TS=unlist(parLapply(cl, X=1:p, logistic_TS_util_fun, 
                      x_sp = data$X, y_sp = data$Y,
                      betaNull = 0, multiTest = FALSE, nlambda = 100, 
                      alpha = 1,family = "binomial", lambda.min.ratio = 0.001,
                      norm = 'L2.squared'))
  
  stopCluster(cl)
  
  result_L2[i] = all( TS[1:3] %in% sort(TS, decreasing = TRUE)[1:n-1] )
  
 
  cat("Now running:", i, '\n')
} 
  
  return( list(L1 = mean(result_L1), L2 = mean(result_L2)) )
}



screen_simu = function(n, p, signal, rho, iter=200, norm = 'L1'){
  
  result_L1 = c(); result_L2 = c()
  for (i in 1:iter){
    
    #n = 20;p = 100;rho = 0.5
    data = depenDesign(n = n, p = p, beta = c(rep(signal,3), rep(0, p-3)), rho)
    
    n_threads = detectCores()
    cl = makeCluster(n_threads, type = "FORK")
    
    TS=unlist(parLapply(cl, X=1:p, TS_util_fun, 
                        x_sp = data$X, y_sp = data$Y,
                        betaNull = 0, multiTest = FALSE, path.method = "lars",
                        norm = 'L1', normalize = TRUE, intercept = FALSE))
    stopCluster(cl)
    
    result_L1[i] = all( TS[1:3] %in% sort(TS, decreasing = TRUE)[1:n-1] )
    
    cl = makeCluster(n_threads, type = "FORK") 
    TS=unlist(parLapply(cl, X=1:p, TS_util_fun, 
                        x_sp = data$X, y_sp = data$Y,
                        betaNull = 0, multiTest = FALSE, path.method = "lars",
                        norm = 'L2.squared', normalize = TRUE, intercept = FALSE))
    stopCluster(cl)
    
    result_L2[i] = all( TS[1:3] %in% sort(TS, decreasing = TRUE)[1:n-1] )
    
    
    cat("Now running:", i, '\n')
  } 
  
  return( list(L1 = mean(result_L1), L2 = mean(result_L2)) )
}







sis_simu = function(n, p, signal, rho, iter=200){
  
  result_sis = c(); result_isis = c()
  for (i in 1:iter){
    
    #n = 20;p = 100;rho = 0.5
    data = depenDesign(n = n, p = p, beta = c(rep(signal,3), rep(0, p-3)), rho)
    sis_index = SIS(x = data$X, y = data$Y, family = 'gaussian', iter = FALSE, nsis = n-1)$sis.ix0
    isis_index = SIS(x = data$X, y = data$Y, family='gaussian', tune='bic', nsis = n-1)$ix0
    
    result_sis[i] = all( 1:3 %in% sis_index)
    result_isis[i] = all( 1:3 %in%  isis_index)
    
    
    cat("Now running:", i, '\n')
  }
  
  return( list(sis = mean(result_sis), isis = mean(result_isis)) )
}



logistic_sis_simu = function(n, p, signal, rho, iter=200){
  
  result_sis = c(); result_isis = c()
  for (i in 1:iter){
    
    #n = 20;p = 100;rho = 0.5
    data = binomDesign(n = n, p = p, beta = c(rep(signal,3), rep(0, p-3)), rho, intercept = 0.5)
    sis_index = SIS(x = data$X, y = data$Y, family = 'binomial', iter = FALSE, nsis = n-1)$sis.ix0
    isis_index = SIS(x = data$X, y = data$Y, family='binomial', tune='bic', nsis = n-1)$ix0
    
    result_sis[i] = all( 1:3 %in% sis_index)
    result_isis[i] = all( 1:3 %in%  isis_index)
    
    
    cat("Now running:", i, '\n')
  }
  
  return( list(sis = mean(result_sis), isis = mean(result_isis)) )
}


#rate_n_20_p_1000_l2 = c(); rate_n_20_p_1000_l1 = c();
#i = 1
#for(rho in c(0, 0.1, 0.5, 0.9)){
#  rate_n_20_p_1000_l2[i] = screen_simu(n = 20, p = 1000, rho = rho, norm = "L2.squared")
#  rate_n_20_p_1000_l1[i] = screen_simu(n = 20, p = 1000, rho = rho, norm = "L1")
#  i = i + 1
#}

