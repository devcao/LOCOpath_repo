.libPaths("~/R")
require(hdi)
require(LOCOpath)

data(riboflavin)

TS = c()
# for (i in 1:length(riboflavin$y)){
#   TS[i] = 
#     TS_util_fun(x = riboflavin$x, y = riboflavin$y, which.covariate = 1, betaNull = 0, multiTest = FALSE, path.method = "lars",
#                norm = "L1", normalize = TRUE, intercept = FALSE)
#   print(TS[i])  
# }

require(doMC)


TS_util_fun = function(x_sp, y_sp, which.covariate = 1, betaNull = 0, multiTest = FALSE, path.method = "lars",
                       norm = "L1", normalize = TRUE, intercept = FALSE){
  return(
    ExactPath.TS(X = x_sp, Y = y_sp, which.covariate = which.covariate, betaNull = betaNull, 
                 multiTest = multiTest, path.method = path.method,
                 norm = norm, normalize = normalize, intercept = intercept)
  ) 
}

n_threads = -1
if(n_threads == -1){n_threads = detectCores()}
cl = makeCluster(n_threads, type = "FORK")


TS=unlist(parLapply(cl, X=1:4088, TS_util_fun, x_sp = riboflavin$x, y_sp = riboflavin$y, betaNull = 0, multiTest = FALSE, path.method = "lars",
                         norm = "L1", normalize = TRUE, intercept = FALSE))


save(TS, file = "TS_ribo.RData")

