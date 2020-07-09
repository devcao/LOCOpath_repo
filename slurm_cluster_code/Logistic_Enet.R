
Rcpp::sourceCpp('lognet.cpp')

Logistic_Enet = function(X, Y, alpha = 1, nlambda = 100, constraint = FALSE,
                         which.covariate = 1, betaNull = 0,
                         lambda.min.ratio = 0.001, epsilon = 0.0001){
  n = dim(X)[1]; p = dim(X)[2]
  X = scale(X); 
  
  max_lam = max(abs( t(X) %*% Y  / n)) / alpha * 1.10 * n
  
  lambda = seq(lambda.min.ratio*max_lam, max_lam, length=nlambda)
  beta_cd = matrix(0,nlambda,p)
  beta0_cd = c()
  
  if (constraint) {
    for(l in 1:nlambda){
      #if(l %% 10 == 0){cat('Now computing', l, lambda[l], ' ')}
      bhat = logistic_enet_constraint(Yr = Y, Xr = cbind(rep(1,n),X), 
                                      which_cov = which.covariate, betaNull = betaNull,
                                      lambda = lambda[l], theta = alpha, 
                                      binitr = rep(0,p+1), delta = epsilon)$b
      beta_cd[l, ] = bhat[-1]
      beta0_cd = bhat[1]
    }
  }else{
    for(l in 1:nlambda){
      #if(l %% 10 == 0){cat('Now computing', l, lambda[l], ' ')}
      bhat = logistic_enet(Yr = Y, Xr = cbind(rep(1,n),X), 
                           which_cov = which.covariate, betaNull = betaNull,
                                      lambda = lambda[l], theta = alpha, 
                                      binitr = rep(0,p+1), delta = epsilon)$b
      beta_cd[l, ] = bhat[-1]
      beta0_cd = bhat[1]
    }
  }  
  return(list(lambda = lambda/n, beta_hat = beta_cd, intercept = beta0_cd))  
}
