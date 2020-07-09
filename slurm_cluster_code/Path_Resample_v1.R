
Path.Resample = function (X, Y, which.covariate, betaNull, multiTest, B = 500, 
          parallel = FALSE, exact = TRUE, beta.init = "adaptive", beta.true = beta, 
          ...) 
{
  n = nrow(X)
  p = ncol(X)
  rej = matrix(0, length(which.covariate), 4)
  pval = numeric()
  TS = Path.TS(exact = exact, X = X, Y = Y, which.covariate = which.covariate, 
               betaNull = betaNull, multiTest = multiTest, ...)
  # if (p >= n) {
    if (beta.init == "adaptive") {
	print('all used adaptive')
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
  # }
  # else {
  #  bhat = ginv(t(X) %*% X) %*% t(X) %*% Y
  #  residual = Y - X %*% bhat
  # }
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
    TS_null = Path.Resample.Process(X = X, Y = Y, multiTest = multiTest, 
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
