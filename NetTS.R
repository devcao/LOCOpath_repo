
# This file contains all the functions calculating our LOCO test statistuc


######## Compute LOCO statistic for linear regression using glmnet ##########
ExactNet.TS <- function(X, Y, which.covariate, betaNull = 0, multiTest = FALSE,
                        nlambda = 100,  alpha = 1, family = "gaussian", 
                        norm = 'L2.squared'){
  # Compute LOCO statistic for linear regression using glmnet
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
  # nlambda: int, argument for glmnet function, gives the length of lambda sequence
  # alpha: float, between 0 and 1, argument for glmnet function, 0 for L2 and 1 for L1 norm
  # family: argument for glmnet function, specify the regression type of GLM, 
  #         you can ignore if you are doing linear regression
  # norm: the norm of LOCO statistic, shoule be ('L2.squared', 'L1', 'L_inf')
  # Returns:
  #       Numeric/Vector/List, data type matches with which.covaraite 
  #       returns the value of our LOCO statistic
  #
  #
  n <- nrow(X)
  p <- ncol(X)
  
  ################# Here we go ! ! ! #######################################################  
  
  ##### extract the path #######
  
  if (length(which.covariate) != length(betaNull)){
    
    stop("Length of variables being tested must equal the length of their Null hypothesis")
    
  }#else{
  #  multiTest <- length(which.covariate) > 1
  #}
  l = 1
  TS = c()
  
  for (j in which.covariate){  
  # simultanoues testing 
  if(multiTest){
    
    #### ajust Y by using Y - X_j * beta_j_NULL - X_k * beta_k_NULL - ...
    adjust.X <- rowSums( t( t(X) * betaNull[[l]]  ) )     ## ? what the fuck is this?
    adjust.Y <- as.vector( Y - adjust.X - mean(Y))
    X_scaled = scale(X)
    
    ## TODO
    max_lam = NA
    
    net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = FALSE, standardize = TRUE)
    lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
    
    net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = FALSE, standardize = TRUE)
    lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    
    
    # indivdual testing   
  }else if(!multiTest){  #indivdual test
    
    adjust.Y <- as.vector( Y - betaNull[l] * X[,j] - mean(Y) )
    X_scaled <- scale(X)
    
    ## TODO
    net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = FALSE, standardize = TRUE)
    lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
    
    net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = FALSE, standardize = TRUE)
    lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    
    
  }else{
    stop("wrong input of multiTest, must be boolean")
  } 
  
  
  minLam <- min(lambda.hat, lambda.j.hat)
  
  lambda.hat = lambda.hat[lambda.hat > minLam]
  
  lambda.j.hat = lambda.hat[lambda.hat > minLam]
  
  union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
  
  beta.j.hat <- Matrix(0, length(union.lambda), p)  # sparse Matrix
  
  beta.hat <- predict(net.out, s=union.lambda, type="coef")
  beta.hat <- t( beta.hat[-1, ] )  # we don't need intercept
  
  
  beta.j.tmp <- predict(net.j.out, s=union.lambda, type="coef")
  beta.j.hat[, -j] <- t( beta.j.tmp[-1, ] ) # we don't need intercep
  
  TS.k <- numeric()
  M <- length(union.lambda)
  
  for (k in 1:p){
    
    # get beta.hat and beta.j.hat at all values of lambda in union lambda:
    #beta.hat.union.lambda[,k] <- 
    #  approx(x = lambda.hat, y = beta.hat[,k], xout = union.lambda,yright=0)$y 
    
    #beta.j.hat.union.lambda[,k] <- 
    #  approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = union.lambda,yright=0)$y
    # get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
    
    delta <- (beta.hat[,k] - beta.j.hat[,k])
    
    if (norm == "L2.squared"){
      TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
    }else if (norm == "L1"){
      TS.k[k] <- 0.5*sum(diff(union.lambda)*(abs(delta[-M]) + abs(delta[-1])))
    }else if (norm == "L_inf"){
      TS.k[k] <- max(abs(delta))
    }
  }
  
  if (norm == "L_inf"){
    TS[l] <- max(TS.k)
  }else{
    TS[l] <- sum(TS.k)
    
  }
  
  l=l+1
  }
  return(TS)
  #return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}

######## Parallel computing version of ExactNet.TS, won't be used alone ##########
ExactNet.TS.Para <- function(mat, ...){
  # Parallel computing version of ExactNet.TS
  # mat: matrix, the last column gives response Y, 1st to the 2nd to last columns gives X
  # Other arguments same as ExactNet.TS
  n = nrow(mat)
  p = ncol(mat) - 1
  X = mat[,1:p] 
  Y = mat[,p+1] 
  
  return(ExactNet.TS(X,Y,...))
  
}




######## Compute LOCO statistic for Logistic regression using glmnet ##########
ExactNet.TS.Logistic <- function(X, Y, which.covariate, betaNull = 0, multiTest = FALSE,
                        nlambda = 100,  alpha = 1, family = "binomial", lambda.min.ratio = 0.001,
                        norm = 'L2.squared'){
  # Compute LOCO statistic for Logistic regression using glmnet
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
  #                   For testing non-zero betaNull, check code ExactNet.TS.Logistic.Con
  # multiTest: boolean, FALSE then do single test, TRUE then do simultanoues test
  #                     must match with which.covarite and betaNull
  # nlambda: int, argument for glmnet function, gives the length of lambda sequence
  # alpha: float, between 0 and 1, argument for glmnet function, 0 for L2 and 1 for L1 norm
  # family: argument for glmnet function, specify the regression type of GLM, 
  #         you can ignore if you are doing Logistic regression
  # lambda.min.ratio: argument for glmnet function, 
  #                   specify the lambda_min = lambda.min.ratio * lambda_max
  # norm: the norm of LOCO statistic, shoule be ('L2.squared', 'L1', 'L_inf')
  # Returns:
  #       Numeric/Vector/List, data type matches with which.covaraite 
  #       returns the value of our LOCO statistic
  #
  #
  
  n <- nrow(X)
  p <- ncol(X)
  
  ################# Here we go ! ! ! #######################################################  
  
  ##### extract the path #######
  
  if (length(which.covariate) != length(betaNull)){
    stop("Length of variables being tested must equal the length of their Null hypothesis")
  }#else{
  #  multiTest <- length(which.covariate) > 1
  #}
  
  l = 1
  TS = c()
  
  for (j in which.covariate){
    
  # simultanoues testing 
  if(multiTest){
    
    #### ajust Y by using Y - X_j * beta_j_NULL - X_k * beta_k_NULL - ...
    adjust.X <- rowSums( t( t(X) * betaNull[[l]]  ) )     ## ? what the fuck is this?
    adjust.Y <- as.vector( Y - adjust.X )
    X_scaled = scale(X)
    
    net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
    
    net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    
    
    # indivdual testing   
  }else if(!multiTest){  #indivdual test
    
    #adjust.Y <- as.vector( Y - betaNull * X[,which.covariate] - mean(Y) )
    adjust.Y <- as.vector( Y - betaNull[l] * X[,j] )
    
    X_scaled <- scale(X)
    

    net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
    
    net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    
    
  }else{
    stop("wrong input of multiTest, must be boolean")
  } 
  
  
  minLam <- min(lambda.hat, lambda.j.hat)
  
  lambda.hat = lambda.hat[lambda.hat > minLam]
  
  lambda.j.hat = lambda.hat[lambda.hat > minLam]
  
  union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
  
  beta.j.hat <- Matrix(0, length(union.lambda), p)  # sparse Matrix
  
  beta.hat <- predict(net.out, s=union.lambda, type="coef")
  beta.hat <- t( beta.hat[-1, ] )  # we don't need intercept
  
  
  beta.j.tmp <- predict(net.j.out, s=union.lambda, type="coef")
  beta.j.hat[, -j] <- t( beta.j.tmp[-1, ] ) # we don't need intercep
  
  TS.k <- numeric()
  M <- length(union.lambda)
  
  for (k in 1:p){
    
    # get beta.hat and beta.j.hat at all values of lambda in union lambda:
    #beta.hat.union.lambda[,k] <- 
    #  approx(x = lambda.hat, y = beta.hat[,k], xout = union.lambda,yright=0)$y 
    
    #beta.j.hat.union.lambda[,k] <- 
    #  approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = union.lambda,yright=0)$y
    # get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
    
    delta <- (beta.hat[,k] - beta.j.hat[,k])
    
    if (norm == "L2.squared"){
      TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
    }else if (norm == "L1"){
      TS.k[k] <- 0.5*sum(diff(union.lambda)*(abs(delta[-M]) + abs(delta[-1])))
    }else if (norm == "L_inf"){
      TS.k[k] <- max(abs(delta))
    }
  }
  
  if (norm == "L_inf"){
    TS[l] <- max(TS.k)
  }else{
    TS[l] <- sum(TS.k)
    
  }
  
  l=l+1
  ## finish for loop
  }
  return(TS)
  #return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}

######## Parallel computing version of ExactNet.TS.Logistic, won't be used alone ##########
ExactNet.TS.Logistic.Para <- function(mat, ...){
  # Parallel computing version of ExactNet.TS.Logistic
  # mat: matrix, the last column gives response Y, 1st to the 2nd to last columns gives X
  # Other arguments same as ExactNet.TS.Logistic
  
  
  n = nrow(mat)
  p = ncol(mat) - 1
  X = mat[,1:p] 
  Y = mat[,p+1] 
  
  return(ExactNet.TS.Logistic(X,Y,...))
  
}






# Compute LOCO statistic for testing beta_j = non_zero_value in Logistic regression 
# Based on our own coordinate descent code
ExactNet.TS.Logistic.Con <- function(X, Y, which.covariate, betaNull = 0, 
                                 nlambda = 100,  alpha = 1, lambda.min.ratio = 0.01,
                                 norm = 'L2.squared'){
  # Compute LOCO statistic for testing beta_j = non_zero_value in Logistic regression 
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
  # nlambda: int, argument for coordinate descent, gives the length of lambda sequence
  # alpha: float, between 0 and 1, argument for coordinate descent, 0 for L2 and 1 for L1 norm
  # lambda.min.ratio: argument for coordinate descent,
  #                   specify the lambda_min = lambda.min.ratio * lambda_max
  # norm: the norm of LOCO statistic, shoule be ('L2.squared', 'L1', 'L_inf')
  # Returns:
  #       Numeric/Vector/List, data type matches with which.covaraite 
  #       returns the value of our LOCO statistic
  #
  #
  n <- nrow(X)
  p <- ncol(X)
  
  ################# Here we go ! ! ! #######################################################  
  
  ##### extract the path #######
  
  if (length(which.covariate) != length(betaNull)){
    stop("Length of variables being tested must equal the length of their Null hypothesis")
  }#else{
  #  multiTest <- length(which.covariate) > 1
  #}
  
  l = 1
  TS = c()
  
  for (j in which.covariate){
    
    
      
      #adjust.Y <- as.vector( Y - betaNull * X[,which.covariate] - mean(Y) )
      adjust.Y <- as.vector( Y )#- betaNull[l] * X[,j] )
      #X_scaled <- scale(X)
      
      net.j.out = Logistic_Enet(X = X, Y = adjust.Y, constraint = TRUE,
                                which.covariate = j, betaNull = betaNull[l],
                                nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, epsilon = 0.001)
      
      net.out = Logistic_Enet(X = X, Y = adjust.Y, constraint = FALSE, which.covariate = j, betaNull = betaNull[l],
                              nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, epsilon = 0.001)
      
      #net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
      lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
      
      #net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
      lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
      union.lambda = lambda.hat
      
      beta.hat = net.out$beta_hat
      beta.j.hat = net.j.out$beta_hat  
      
      
    #}else{
    #  stop("wrong input of multiTest, must be boolean")
    #} 
    
    
    #minLam <- min(lambda.hat, lambda.j.hat)
    #lambda.hat = lambda.hat[lambda.hat > minLam]
    #lambda.j.hat = lambda.hat[lambda.hat > minLam]
    #union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
    
    #beta.j.hat <- Matrix(0, length(union.lambda), p)  # sparse Matrix
    
    #beta.hat <- predict(net.out, s=union.lambda, type="coef")
    #beta.hat <- t( beta.hat[-1, ] )  # we don't need intercept
    
    
    #beta.j.tmp <- predict(net.j.out, s=union.lambda, type="coef")
    #beta.j.hat[, -j] <- t( beta.j.tmp[-1, ] ) # we don't need intercep
    
    TS.k <- numeric()
    M <- length(union.lambda)
    
    for (k in 1:p){
      
      # get beta.hat and beta.j.hat at all values of lambda in union lambda:
      #beta.hat.union.lambda[,k] <- 
      #  approx(x = lambda.hat, y = beta.hat[,k], xout = union.lambda,yright=0)$y 
      
      #beta.j.hat.union.lambda[,k] <- 
      #  approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = union.lambda,yright=0)$y
      # get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
      
      delta <- (beta.hat[,k] - beta.j.hat[,k])
      
      if (norm == "L2.squared"){
        TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
      }else if (norm == "L1"){
        TS.k[k] <- 0.5*sum(diff(union.lambda)*(abs(delta[-M]) + abs(delta[-1])))
      }else if (norm == "L_inf"){
        TS.k[k] <- max(abs(delta))
      }
    }
    
    if (norm == "L_inf"){
      TS[l] <- max(TS.k)
    }else{
      TS[l] <- sum(TS.k)
      
    }
    
    l=l+1
    ## finish for loop
  }
  return(TS)
  #return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}


######## Parallel computing version of ExactNet.TS.Logistic.Con, won't be used alone ##########
ExactNet.TS.Logistic.Con.Para <- function(mat, ...){
  # Parallel computing version of ExactNet.TS.Logistic.Con
  # mat: matrix, the last column gives response Y, 1st to the 2nd to last columns gives X
  # Other arguments same as ExactNet.TS.Logistic.Con
  
  n = nrow(mat)
  p = ncol(mat) - 1
  X = mat[,1:p] 
  Y = mat[,p+1] 
  
  return(ExactNet.TS.Logistic.Con(X,Y,...))
  
}






# Compute LOCO statistic for simultanoues testing in Logistic regression 
# Based on group LASSO 
ExactNet.TS.Logistic.Multi <- function(X, Y, which.covariate, betaNull = 0,
                                     nlambda = 100,  lambda.min.ratio = 0.01,
                                     norm = 'L2.squared'){
  # Compute LOCO statistic for simultanoues testing in Logistic regression 
  #
  # Args:
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate: should be a list of vectors
  #                  E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  # betaNULL: ignroe or just use default, won't be used in the actual code
  # nlambda: int, argument for gglasso, gives the length of lambda sequence
  # lambda.min.ratio: argument for gglasso,
  #                   specify the lambda_min = lambda.min.ratio * lambda_max
  # norm: the norm of LOCO statistic, shoule be ('L2.squared', 'L1', 'L_inf')
  # Returns: returns the value of our LOCO statistic
  #
  # 
  n <- nrow(X)
  p <- ncol(X)
  
  ################# Here we go ! ! ! #######################################################  
  
  ##### extract the path #######
  
  if (length(which.covariate) != length(betaNull)){
    stop("Length of variables being tested must equal the length of their Null hypothesis")
  }#else{
  #  multiTest <- length(which.covariate) > 1
  #}
  
  l = 1
  TS = c()
  
  for (j in which.covariate){
    
    # simultanoues testing 
    #if(multiTest){
    #### ajust Y by using Y - X_j * beta_j_NULL - X_k * beta_k_NULL - ...
    # adjust.X <- rowSums( t( t(X) * betaNull[[l]]  ) )     ## ? what the fuck is this?
    # adjust.Y <- as.vector( Y - adjust.X )
    # X_scaled = scale(X)
    # net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    # lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
    # net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    # lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    
    # indivdual testing   
    #}else if(!multiTest){  #indivdual test
    
    #adjust.Y <- as.vector( Y - betaNull * X[,which.covariate] - mean(Y) )
    adjust.Y <- 2*as.vector(Y) - 1#- betaNull[l] * X[,j] )
    X_scaled <- scale(X)
    
    #m1 <- gglasso(x=data$X, y=2*data$Y-1,group=group1,loss="logit")
    
    group = rep(0, p)
    
    group[j] = 1
    group[-j] = 2:(p-length(j)+1)
    print(group[1:15])
    #group.j = rep(1, p-length(j))
    
    gglasso.out = gglasso(x = X_scaled, y = adjust.Y, loss = 'logit', group = group, 
                          nlambda = nlambda, lambda.factor = lambda.min.ratio)
    beta.hat = t(gglasso.out$beta)
    lambda.hat = gglasso.out$lambda
    
    gglasso.j.out = gglasso(x = X_scaled[,-j], y = adjust.Y, loss = 'logit', 
                            group = NULL, lambda = gglasso.out$lambda)
    
    beta.j.hat = matrix(0, nrow = dim(beta.hat)[1], ncol = dim(beta.hat)[2])
    beta.j.hat[,-j] = t(gglasso.j.out$beta)
    #net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    lambda.hat <- sort(lambda.hat, decreasing = FALSE)
    #net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    #lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    union.lambda = lambda.hat
    
    #}else{
    #  stop("wrong input of multiTest, must be boolean")
    #} 
    
    
    #minLam <- min(lambda.hat, lambda.j.hat)
    #lambda.hat = lambda.hat[lambda.hat > minLam]
    #lambda.j.hat = lambda.hat[lambda.hat > minLam]
    #union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
    
    #beta.j.hat <- Matrix(0, length(union.lambda), p)  # sparse Matrix
    
    #beta.hat <- predict(net.out, s=union.lambda, type="coef")
    #beta.hat <- t( beta.hat[-1, ] )  # we don't need intercept
    
    
    #beta.j.tmp <- predict(net.j.out, s=union.lambda, type="coef")
    #beta.j.hat[, -j] <- t( beta.j.tmp[-1, ] ) # we don't need intercep
    
    TS.k <- numeric()
    M <- length(union.lambda)
    
    for (k in 1:p){
      
      # get beta.hat and beta.j.hat at all values of lambda in union lambda:
      #beta.hat.union.lambda[,k] <- 
      #  approx(x = lambda.hat, y = beta.hat[,k], xout = union.lambda,yright=0)$y 
      
      #beta.j.hat.union.lambda[,k] <- 
      #  approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = union.lambda,yright=0)$y
      # get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
      
      delta <- (beta.hat[,k] - beta.j.hat[,k])
      
      if (norm == "L2.squared"){
        TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
      }else if (norm == "L1"){
        TS.k[k] <- 0.5*sum(diff(union.lambda)*(abs(delta[-M]) + abs(delta[-1])))
      }else if (norm == "L_inf"){
        TS.k[k] <- max(abs(delta))
      }
    }
    
    if (norm == "L_inf"){
      TS[l] <- max(TS.k)
    }else{
      TS[l] <- sum(TS.k)
      
    }
    
    l=l+1
    ## finish for loop
  }
  return(TS)
  #return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}


######## Parallel computing version of ExactNet.TS.Logistic.Multi, won't be used alone ##########
ExactNet.TS.Logistic.Multi.Para <- function(mat, ...){
  # Parallel computing version of ExactNet.TS.Logistic.Multi
  # mat: matrix, the last column gives response Y, 1st to the 2nd to last columns gives X
  # Other arguments same as ExactNet.TS.Logistic.Multi
  
  n = nrow(mat)
  p = ncol(mat) - 1
  X = mat[,1:p] 
  Y = mat[,p+1] 
  
  return(ExactNet.TS.Logistic.Multi(X,Y,...))
  
}


  # Compute LOCO statistic for Logistic regression using glmnet
  #
  # Args:
  # X,Y: design matrix (matrix) and response vector (vecotr or matrix of 1 column)
  # which.covariate: should be a list of vectors
  #                  E.g, which.covariate = list(c(1,2)), we test simultaneously for beta_1=.. and beta_2= ..
  # betaNULL: ignroe or just use default, won't be used in the actual code
  # nlambda: int, argument for gglasso, gives the length of lambda sequence
  # lambda.min.ratio: argument for gglasso,
  #                   specify the lambda_min = lambda.min.ratio * lambda_max
  # norm: the norm of LOCO statistic, shoule be ('L2.squared', 'L1', 'L_inf')
  # Returns: returns the value of our LOCO statistic
  #



# Compute LOCO statistic for covariance matrix in graphical models
# Based on our own graphical LASSO code
ExactNet.TS.Graph.Sigma <- function(S, c1, c2, n_rho = 10, norm = 'L2.squared',...){
  
  # Compute LOCO statistic for covariance matrix in graphical models
  #
  # Args:
  # S: sample matrix, covariance matrix
  # c1, c2: int, specify the row and column entry of S
  # n_rho: int, gives the length of lambda sequence
  # norm: the norm of LOCO statistic, shoule be ('L2.squared', 'L1', 'L_inf')
  # Returns: returns the value of our LOCO statistic
  #
  #
  p <- ncol(S)
  
  ################# Here we go ! ! ! #######################################################  
  
  ##### extract the path #######
  
  if (length(c1) != length(c2)){
    stop("Length of c1 not equal to length of c2")
  }
  
  l = 1
  TS = c()
  n_iter = length(c1)
  
  for (j in n_iter){
    
    graph.j.out = graphLASSO.path.c(S, n_rho = n_rho, c1=c1[j], c=c2[j])#, ...)
    
    graph.out = graphLASSO.path(S, n_rho = n_rho)#, ...)
    
    #net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    lambda.hat <- sort(graph.out$rholist, decreasing = FALSE)
    
    #net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    
    union.lambda = lambda.hat
    
    beta.hat = trans_lars_matrix(graph.out)
    beta.j.hat = trans_lars_matrix(graph.j.out)
    
    
    #}else{
    #  stop("wrong input of multiTest, must be boolean")
    #} 
    
    
    #minLam <- min(lambda.hat, lambda.j.hat)
    #lambda.hat = lambda.hat[lambda.hat > minLam]
    #lambda.j.hat = lambda.hat[lambda.hat > minLam]
    #union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
    
    #beta.j.hat <- Matrix(0, length(union.lambda), p)  # sparse Matrix
    
    #beta.hat <- predict(net.out, s=union.lambda, type="coef")
    #beta.hat <- t( beta.hat[-1, ] )  # we don't need intercept
    
    
    #beta.j.tmp <- predict(net.j.out, s=union.lambda, type="coef")
    #beta.j.hat[, -j] <- t( beta.j.tmp[-1, ] ) # we don't need intercep
    
    TS.k <- numeric()
    M <- length(union.lambda)
    num_coef = p*(p-1)/2
    
    for (k in 1:num_coef){
      
      # get beta.hat and beta.j.hat at all values of lambda in union lambda:
      #beta.hat.union.lambda[,k] <- 
      #  approx(x = lambda.hat, y = beta.hat[,k], xout = union.lambda,yright=0)$y 
      
      #beta.j.hat.union.lambda[,k] <- 
      #  approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = union.lambda,yright=0)$y
      # get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
      
      delta <- (beta.hat[,k] - beta.j.hat[,k])
      
      if (norm == "L2.squared"){
        TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
      }else if (norm == "L1"){
        TS.k[k] <- 0.5*sum(diff(union.lambda)*(abs(delta[-M]) + abs(delta[-1])))
      }else if (norm == "L_inf"){
        TS.k[k] <- max(abs(delta))
      }
    }
    
    if (norm == "L_inf"){
      TS[l] <- max(TS.k)
    }else{
      TS[l] <- sum(TS.k)
      
    }
    
    l=l+1
    ## finish for loop
  }
  return(TS)
  #return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}


# Compute LOCO statistic for precision matrix in graphical models
# Based on glasso code
ExactNet.TS.Graph <- function(S, c1, c2, n_rho = 10, norm = 'L2.squared', large_pen=10e4, ...){
  
  # Compute LOCO statistic for precision matrix in graphical models
  #
  # Args:
  # S: sample matrix, covariance matrix
  # c1, c2: int, specify the row and column entry of S
  # n_rho: int, gives the length of lambda sequence
  # norm: the norm of LOCO statistic, shoule be ('L2.squared', 'L1', 'L_inf')
  # large_pen: the penalty on the (c1,c2) entry of S, should be large number
  #            just use the default unless you are experimenting 
  # ..., other arguments for glassopath.c and glassopath.g function
  # Returns: returns the value of our LOCO statistic
  #
  #
  
  p <- ncol(S)
  
  ################# Here we go ! ! ! #######################################################  
  
  ##### extract the path #######
  
  if (length(c1) != length(c2)){
    stop("Length of c1 not equal to length of c2")
  }
  
  l = 1
  TS = c()
  n_iter = length(c1)
  
  
  rholist = seq(max(abs(S))/n_rho, max(abs(S)), length = n_rho)
  
  for (j in n_iter){
    
    # graph.j.out = graphLASSO.path.c(S, n_rho = n_rho, c1=c1[j], c=c2[j])#, ...)
    # graph.out = graphLASSO.path(S, n_rho = n_rho)#, ...)
    
    graph.out = glassopath.g(S = S, rholist = rholist,...)
    
    graph.j.out = glassopath.c(S = S, c1 = c1, c2 = c2, rholist = rholist, large_pen=large_pen, ...)
    #net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    lambda.hat <- sort(graph.out$rholist, decreasing = FALSE)
    
    #net.j.out <- glmnet(X_scaled[, -j], adjust.Y, nlambda = nlambda, alpha = alpha, lambda.min.ratio = lambda.min.ratio, family = family, intercept = TRUE, standardize = TRUE)
    
    union.lambda = lambda.hat
    
    beta.hat = trans_lars_matrix(graph.out, use.glasso = FALSE)
    beta.j.hat = trans_lars_matrix(graph.j.out, use.glasso = FALSE)
    
    
    #}else{
    #  stop("wrong input of multiTest, must be boolean")
    #} 
    
    
    #minLam <- min(lambda.hat, lambda.j.hat)
    #lambda.hat = lambda.hat[lambda.hat > minLam]
    #lambda.j.hat = lambda.hat[lambda.hat > minLam]
    #union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
    
    #beta.j.hat <- Matrix(0, length(union.lambda), p)  # sparse Matrix
    
    #beta.hat <- predict(net.out, s=union.lambda, type="coef")
    #beta.hat <- t( beta.hat[-1, ] )  # we don't need intercept
    
    
    #beta.j.tmp <- predict(net.j.out, s=union.lambda, type="coef")
    #beta.j.hat[, -j] <- t( beta.j.tmp[-1, ] ) # we don't need intercep
    
    TS.k <- numeric()
    M <- length(union.lambda)
    num_coef = p*(p-1)/2
    
    for (k in 1:num_coef){
      
      # get beta.hat and beta.j.hat at all values of lambda in union lambda:
      #beta.hat.union.lambda[,k] <- 
      #  approx(x = lambda.hat, y = beta.hat[,k], xout = union.lambda,yright=0)$y 
      
      #beta.j.hat.union.lambda[,k] <- 
      #  approx(x = lambda.j.hat, y = beta.j.hat[,k], xout = union.lambda,yright=0)$y
      # get absolute difference between beta.hat and beta.j.hat at all values of lambda in union lambda
      
      delta <- (beta.hat[,k] - beta.j.hat[,k])
      
      if (norm == "L2.squared"){
        TS.k[k] <- sum(diff(union.lambda)*(delta[-M]^2 + diff(delta) * delta[-M] + (1/3) * diff(delta)^2 ))
      }else if (norm == "L1"){
        TS.k[k] <- 0.5*sum(diff(union.lambda)*(abs(delta[-M]) + abs(delta[-1])))
      }else if (norm == "L_inf"){
        TS.k[k] <- max(abs(delta))
      }
    }
    
    if (norm == "L_inf"){
      TS[l] <- max(TS.k)
    }else{
      TS[l] <- sum(TS.k)
      
    }
    
    l=l+1
    ## finish for loop
  }
  return(TS)
  #return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}



######## Parallel computing version of ExactNet.TS.Graph, won't be used alone ##########
ExactNet.TS.Graph.Para <- function(mat, ...){
  
  # Parallel computing version of ExactNet.TS.Graph
  # mat: matrix of the covariance matrix
  # Other arguments same as ExactNet.TS.Graph
  
  S = mat
  return(ExactNet.TS.Graph(S,...))
  
}






