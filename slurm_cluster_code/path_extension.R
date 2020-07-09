ExactPath.TS_v1=function (X, Y, which.covariate, betaNull, multiTest, path.method = "lars", 
          norm = "L2.squared", normalize = TRUE, intercept = FALSE, full_path=FALSE) 
{
  n = nrow(X)
  p = ncol(X)
  l <- 1
  TS <- numeric(length(which.covariate))
  for (j in which.covariate) {
    if (path.method == "lars") {
      if (multiTest & is.list(which.covariate) & is.list(betaNull)) {
        adjust.X = rowSums(t(apply(X[, j], 1, function(x) {
          return(x * betaNull[[l]])
        })))
        newY = Y - adjust.X
        X.sc = scale(X)
        lars.out <- lars(X.sc, newY, type = "lasso", 
                         intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.hat <- sort(c(lars.out$lambda, 0), decreasing = FALSE)
        beta.hat <- coef(lars.out)[seq(length(lambda.hat), 
                                       1, -1), ]
        lars.j.out <- lars(X.sc[, -j], newY, type = "lasso", 
                           intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.j.hat <- sort(c(lars.j.out$lambda, 0), 
                             decreasing = FALSE)
        beta_val <- coef(lars.j.out)[seq(length(lambda.j.hat), 
                                         1, -1), ]
      }
      else if ((!multiTest)) {
        newY = Y - betaNull[l] * X[, j]
        X.sc = scale(X)
        lars.out <- lars(X.sc, newY, type = "lasso", 
                         intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.hat <- sort(c(lars.out$lambda, 0), decreasing = FALSE)
        beta.hat <- coef(lars.out)[seq(length(lambda.hat), 
                                       1, -1), ]
        lars.j.out <- lars(X.sc[, -j], newY, type = "lasso", 
                           intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.j.hat <- sort(c(lars.j.out$lambda, 0), 
                             decreasing = FALSE)
        beta_val <- coef(lars.j.out)[seq(length(lambda.j.hat), 
                                         1, -1), ]
      }
      else {
        stop("wrong input")
      }
    }
    else if (path.method == "plus.lasso") {
      if (multiTest & is.list(which.covariate) & is.list(betaNull)) {
        adjust.X = rowSums(t(apply(X[, j], 1, function(x) {
          return(x * betaNull[[l]])
        })))
        newY = Y - adjust.X
        X.sc = scale(X)
        lars.out <- plus(X.sc, as.vector(newY), method = "lasso", 
                         intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.hat <- sort(lars.out$lam.path, decreasing = FALSE)
        beta.hat <- coef(lars.out, lam = lambda.hat)
        lars.j.out <- plus(X.sc[, -j], as.vector(newY), 
                           method = "lasso", intercept = intercept, use.Gram = FALSE, 
                           normalize = normalize)
        lambda.j.hat <- sort(lars.j.out$lam.path, decreasing = FALSE)
        beta_val <- coef(lars.j.out, lam = lambda.j.hat)
      }
      else if ((!multiTest)) {
        newY = Y - betaNull[l] * X[, j]
        X.sc = scale(X)
        lars.out <- plus(X.sc, as.vector(newY), method = "lasso", 
                         intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.hat <- sort(lars.out$lam.path, decreasing = FALSE)
        beta.hat <- coef(lars.out, lam = lambda.hat)
        lars.j.out <- plus(X.sc[, -j], as.vector(newY), 
                           method = "lasso", intercept = intercept, use.Gram = FALSE, 
                           normalize = normalize)
        lambda.j.hat <- sort(lars.j.out$lam.path, decreasing = FALSE)
        beta_val <- coef(lars.j.out, lam = lambda.j.hat)
      }
      else {
        stop("wrong input")
      }
    }
    else if (path.method == "plus.mc+") {
      if (multiTest & is.list(which.covariate) & is.list(betaNull)) {
        adjust.X = rowSums(t(apply(X[, j], 1, function(x) {
          return(x * betaNull[[l]])
        })))
        newY = Y - adjust.X
        X.sc = scale(X)
        lars.out <- plus(X.sc, as.vector(newY), method = "mc+", 
                         intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.hat <- sort(lars.out$lam.path, decreasing = FALSE)
        beta.hat <- coef(lars.out, lam = lambda.hat)
        lars.j.out <- plus(X.sc[, -j], as.vector(newY), 
                           method = "mc+", intercept = intercept, use.Gram = FALSE, 
                           normalize = normalize)
        lambda.j.hat <- sort(lars.j.out$lam.path, decreasing = FALSE)
        beta_val <- coef(lars.j.out, lam = lambda.j.hat)
      }
      else if ((!multiTest)) {
        newY = Y - betaNull[l] * X[, j]
        X.sc = scale(X)
        lars.out <- plus(X.sc, as.vector(newY), method = "mc+", 
                         intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.hat <- sort(lars.out$lam.path, decreasing = FALSE)
        beta.hat <- coef(lars.out, lam = lambda.hat)
        lars.j.out <- plus(X.sc[, -j], as.vector(newY), 
                           method = "mc+", intercept = intercept, use.Gram = FALSE, 
                           normalize = normalize)
        lambda.j.hat <- sort(lars.j.out$lam.path, decreasing = FALSE)
        beta_val <- coef(lars.j.out, lam = lambda.j.hat)
      }
      else {
        stop("wrong input")
      }
    }
    else if (path.method == "plus.scad") {
      if (multiTest & is.list(which.covariate) & is.list(betaNull)) {
        adjust.X = rowSums(t(apply(X[, j], 1, function(x) {
          return(x * betaNull[[l]])
        })))
        newY = Y - adjust.X
        X.sc = scale(X)
        lars.out <- plus(X.sc, as.vector(newY), method = "scad", 
                         intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.hat <- sort(lars.out$lam.path, decreasing = FALSE)
        beta.hat <- coef(lars.out, lam = lambda.hat)
        lars.j.out <- plus(X.sc[, -j], as.vector(newY), 
                           method = "scad", intercept = intercept, use.Gram = FALSE, 
                           normalize = normalize)
        lambda.j.hat <- sort(lars.j.out$lam.path, decreasing = FALSE)
        beta_val <- coef(lars.j.out, lam = lambda.j.hat)
      }
      else if ((!multiTest)) {
        newY = Y - betaNull[l] * X[, j]
        X.sc = scale(X)
        lars.out <- plus(X.sc, as.vector(newY), method = "scad", 
                         intercept = intercept, use.Gram = FALSE, normalize = normalize)
        lambda.hat <- sort(lars.out$lam.path, decreasing = FALSE)
        beta.hat <- coef(lars.out, lam = lambda.hat)
        lars.j.out <- plus(X.sc[, -j], as.vector(newY), 
                           method = "scad", intercept = intercept, use.Gram = FALSE, 
                           normalize = normalize)
        lambda.j.hat <- sort(lars.j.out$lam.path, decreasing = FALSE)
        beta_val <- coef(lars.j.out, lam = lambda.j.hat)
      }
      else {
        stop("wrong input")
      }
    }
    if (lambda.hat[1] != lambda.j.hat[1]) {
      leftmost = c(lambda.hat[1], lambda.j.hat[1])
      whichone = which.max(leftmost)
      if (whichone == 1) {
        lambda.j.hat <- lambda.j.hat[lambda.j.hat >= 
                                       lambda.hat[1]]
        lambda.j.hat <- c(lambda.hat[1], lambda.j.hat)
      }
      else {
        lambda.hat <- lambda.hat[lambda.hat >= lambda.j.hat[1]]
        lambda.hat <- c(lambda.j.hat[1], lambda.hat)
      }
      beta.hat <- coef(lars.out, lam = lambda.hat)
      beta_val <- coef(lars.j.out, lam = lambda.j.hat)
    }
    new_beta <- matrix(0, dim(beta_val)[1], p)
    new_beta[, -j] <- beta_val
    beta.j.hat <- new_beta
    union.lambda <- sort(unique(c(lambda.hat, lambda.j.hat)), 
                         decreasing = FALSE)
    M <- length(union.lambda)
    beta.hat.union.lambda <- beta.j.hat.union.lambda <- matrix(NA, 
                                                               length(union.lambda), p)
    TS.k <- numeric()
    if(full_path){
    for (k in 1:p) {
      beta.hat.union.lambda[, k] <- approx(x = lambda.hat, 
                                           y = beta.hat[, k], xout = union.lambda, yright = 0)$y
      beta.j.hat.union.lambda[, k] <- approx(x = lambda.j.hat, 
                                             y = beta.j.hat[, k], xout = union.lambda, yright = 0)$y
      delta <- (beta.hat.union.lambda[, k] - beta.j.hat.union.lambda[, 
                                                                     k])
      if (norm == "L2.squared") {
        TS.k[k] <- sum(diff(union.lambda) * (delta[-M]^2 + 
                                               diff(delta) * delta[-M] + (1/3) * diff(delta)^2))
      }
      else if (norm == "L1") {
        TS.k[k] <- 0.5 * sum(diff(union.lambda) * (abs(delta[-M]) + 
                                                     abs(delta[-1])))
      }
      else if (norm == "L2") {
        TS.k[k] <- sum(diff(union.lambda) * (delta[-M]^2 + 
                                               diff(delta) * delta[-M] + (1/3) * diff(delta)^2))
      }
      else if (norm == "L_inf") {
        TS.k[k] <- max(abs(delta))
      }
    }
    }else{
      k = j
      beta.hat.union.lambda[, k] <- approx(x = lambda.hat, 
                                           y = beta.hat[, k], xout = union.lambda, yright = 0)$y
      beta.j.hat.union.lambda[, k] <- approx(x = lambda.j.hat, 
                                             y = beta.j.hat[, k], xout = union.lambda, yright = 0)$y
      delta <- (beta.hat.union.lambda[, k] - beta.j.hat.union.lambda[, 
                                                                     k])
      if (norm == "L2.squared") {
        TS.k[1] <- sum(diff(union.lambda) * (delta[-M]^2 + 
                                               diff(delta) * delta[-M] + (1/3) * diff(delta)^2))
      }
      else if (norm == "L1") {
        TS.k[1] <- 0.5 * sum(diff(union.lambda) * (abs(delta[-M]) + 
                                                     abs(delta[-1])))
      }
      else if (norm == "L2") {
        TS.k[1] <- sum(diff(union.lambda) * (delta[-M]^2 + 
                                               diff(delta) * delta[-M] + (1/3) * diff(delta)^2))
      }
      else if (norm == "L_inf") {
        TS.k[1] <- max(abs(delta))
      }
      
    }
    if (norm == "L_inf") {
      TS[l] <- max(TS.k)
    }
    else if (norm == "L2") {
      TS[l] <- sqrt(sum(TS.k))
    }
    else {
      TS[l] <- sum(TS.k)
    }
    l <- l + 1
  }
  return(TS)
}

Path.TS_v1 = function (exact = TRUE, ...) 
{
  if (exact) {
    return(ExactPath.TS_v1(...))
  }
  else {
    return(ApproxPath.TS(...))
  }
}

ExactPath.TS.Para_v1=function (mat, ...) 
{
  n = nrow(mat)
  p = ncol(mat) - 1
  X = mat[, 1:p]
  Y = mat[, p + 1]
  return(ExactPath.TS_v1(X, Y, ...))
}

Path.TS.Para_v1 = function (exact = TRUE, ...) 
{
  if (exact) {
    return(ExactPath.TS.Para_v1(...))
  }
  else {
    return(ApproxPath.TS.Para(...))
  }
}


Path.Resample_v1=function (X, Y, which.covariate, betaNull, multiTest, B = 500, 
          parallel = FALSE, exact = TRUE, beta.init = "adaptive", beta.true = beta, 
          ...) 
{
  n = nrow(X)
  p = ncol(X)
  rej = matrix(0, length(which.covariate), 4)
  pval = numeric()
  TS = Path.TS_v1(exact = exact, X = X, Y = Y, which.covariate = which.covariate, 
               betaNull = betaNull, multiTest = multiTest, ...)
  if (beta.init == "adaptive") {
    bhat = adalasso(X = X, y = Y, k = 10, use.Gram = FALSE, 
                    both = TRUE, intercept = FALSE)$coefficients.adalasso
  }
  else if (beta.init == "de-sparse") {
    bhat = as.vector(lasso.proj(X, Y, standardize = TRUE, 
                                parallel = TRUE, ncores = 40)$bhat)
  }
  else if (beta.init == "MC+") {
    bhat = coef(cv.ncvreg(X = X, y = Y, penalty = "MCP", 
                          family = "gaussian", nfold = 10))[-1]
  }
  else if (beta.init == "SCAD") {
    bhat = coef(cv.ncvreg(X = X, y = Y, penalty = "SCAD", 
                          family = "gaussian", nfold = 10))[-1]
  }
  else if (beta.init == "Truth") {
    bhat = beta.true
  }
  residual = Y - X %*% bhat
  count = 1
  for (wc_cov in which.covariate) {
    b.Null = bhat
    if (multiTest) {
      to.which.covariate = list(wc_cov)
      to.betaNull = list(betaNull[[count]])
      b.Null[wc_cov] = betaNull[[count]]
    }
    else {
      to.which.covariate = wc_cov
      to.betaNull = betaNull[count]
      b.Null[wc_cov] = betaNull[count]
    }
    TS_null = Path.Resample.Process_v1(X = X, Y = Y, multiTest = multiTest, 
                                    residual = residual, b.Null = b.Null, betaNull = to.betaNull, 
                                    beta.index = to.which.covariate, B = B, exact = exact, 
                                    parallel = parallel, ...)
    rej[count, 1] = TS[count] > quantile(TS_null, 0.8)
    rej[count, 2] = TS[count] > quantile(TS_null, 0.9)
    rej[count, 3] = TS[count] > quantile(TS_null, 0.95)
    rej[count, 4] = TS[count] > quantile(TS_null, 0.99)
    pval[count] = mean(TS_null > TS[count])
    count = count + 1
  }
  return(list(rej = rej, pval = pval, TS_null = TS_null, TS = TS))
}




Path.Resample.Process_v1=function (X, Y, multiTest, residual, b.Null, beta.index, betaNull, 
          B = 500, exact, parallel = FALSE, ...) 
{
  n = nrow(X)
  p = ncol(X)
  TS_null = numeric()
  if (parallel) {
    mat = list()
    for (bs in 1:B) {
      ind = sample(1:n, replace = TRUE)
      boot_residual = residual[ind]
      Y = X %*% b.Null + boot_residual
      mat[[bs]] = cbind(X, Y)
    }
    no_cores <- detectCores()
    cat("n_cores detected:", no_cores, "\n")
    cl <- makeCluster(no_cores, type = "FORK")
    re_list = parLapply(cl, mat, Path.TS.Para_v1, exact = exact, 
                        multiTest = multiTest, which.covariate = beta.index, 
                        betaNull = betaNull, ...)
    print("Cluster MEM:")
    print(mem_used())
    stopCluster(cl)
    for (bss in 1:B) {
      TS_null[bss] = re_list[[bss]]
    }
    return(TS_null)
  }
  else {
    for (bs in 1:B) {
      ind = sample(1:n, replace = TRUE)
      boot_residual = residual[ind]
      Y = X %*% b.Null + boot_residual
      TS_null[bs] = Path.TS_v1(exact = exact, X = X, Y = Y, 
                            multiTest = multiTest, which.covariate = beta.index, 
                            betaNull = betaNull, ...)
    }
    return(TS_null)
  }
}



Path.Resample.Power_v1=function (n = 100, p = 1000, beta = c(rep(1, 10), rep(0, 990)), 
          rho = 0.5, iter = 500, B = 500, setting = "dep", which.covariate = 1, 
          betaNull = 1, multiTest = FALSE, ...) 
{
  path.power = array(0, dim = c(iter, length(which.covariate), 
                                4))
  for (s in 1:iter) {
    data = dataGen(setting = setting, n = n, p = p, beta = beta, 
                   rho = rho)
    X_sp = data$X
    Y_sp = data$Y
    results = Path.Resample_v1(X = X_sp, Y = Y_sp, which.covariate = which.covariate, 
                            betaNull = betaNull, multiTest = multiTest, B = B, 
                            beta.true = beta, ...)
    print("After Bootstrap:")
    print(mem_used())
    path.power[s, , ] = results$rej
    if (s%%10 == 0) {
      cat("Now Computing:", s, "\n")
    }
  }
  path.power = apply(path.power, c(2, 3), mean)
  return(path.power = path.power)
}



