
# This file includes the R code wrapper for our c++ logistic regression coordinate descent code.

Rcpp::sourceCpp('lognet.cpp')


# R code wrapper for our c++ logistic regression coordinate descent code
Logistic_Enet = function(X, Y, alpha = 1, nlambda = 100, constraint = FALSE,
                         which.covariate = 1, betaNull = 0,
                         lambda.min.ratio = 0.001, epsilon = 0.0001){
  # R code wrapper for our c++ logistic regression coordinate descent code
  # Args:
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate:  numeric, indicating which covariate we will be computing; 
  #                   E.g, which.covariate = 1, we test for beta_1
  #
  # betaNULL: numeric, specify the null hypothesis H0: beta = betaNULL. 
  #                   E.g, if which.covariate = 1, we test for beta_1
  #                   
  # constraint: boolean, if FALSE, fitting coordinate descent with no constriant
  #                      if TRUE, fitting coordinate descent with constriant beta_{which.covariate} = betaNULL
  #                      E.g, if which.covariate = 1, betaNull = 1, itting coordinate descent with constriant beta_1=1.
  # nlambda: int, argument for glmnet function, gives the length of lambda sequence
  # alpha: float, between 0 and 1, argument for glmnet function, 0 for L2 and 1 for L1 norm
  # lambda.min.ratio: argument for glmnet function, 
  #                   specify the lambda_min = lambda.min.ratio * lambda_max
  # epsilon: convergence criterion
  # Returns: a list of:
  #             lambda: lambda sequence
  #             beta_hat: estiamted beta for each lambda
  #             intercept: estimated intercept for each lambda
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
