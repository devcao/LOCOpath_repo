.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_simu")
source("path_power.R")
require(mvnfast)


source('~/hdi_simu/NetResampleLogisticTS.R')
source('~/hdi_simu/NetTS.R')
source('~/hdi_simu/NetResampleTS.R')
source('~/hdi_simu/Logistic_Enet.R')

require(LOCOpath)
require(gglasso)

#prostate = readRDS('Singh2002.rds')


#proj_prostate = lasso.proj(prostate$X, ifelse(prostate$y == 'normal', 1, 0), 
  #                         standardize = TRUE, parallel = TRUE, family = 'binomial', ncores = 40)

#save(proj_prostate, file = 'lasso_prostate.RData')



Logistic_TS_util_fun = function(x_sp, y_sp, which.covariate = 1, betaNull = 0, multiTest = FALSE, norm = 'L1', ...){
  return(
    ExactNet.TS.Logistic(X = x_sp, Y = y_sp, which.covariate = which.covariate, betaNull = betaNull,
                 multiTest = multiTest, ...)
  )
}



glc_amd = readRDS('glc-amd.rds')
str(glc_amd)



#n_threads = detectCores()
#cl = makeCluster(n_threads, type = "FORK")

#TS=unlist(parLapply(cl, X=1:dim(glc_amd$X)[2], Logistic_TS_util_fun, 
 #                   x_sp = glc_amd$X, y_sp = ifelse(glc_amd$y == 'GLC', 1, 0), 
  #                  betaNull = 0, multiTest = FALSE, norm = 'L1'))


#index = which(TS > 0)

#save(TS, file = 'TS_glc_amd.RData')


x = glc_amd$X
y = ifelse(glc_amd$y == 'GLC', 1, 0)

### calculate p-value for genes screened in
obj = list()
for (i in 1:dim(x)[2]){
  print(i)
  try({
    obj[[i]] = Net.Resample.Logistic(X = x, Y = y, which.covariate = i, parallel = TRUE,
                                     betaNull = 0, multiTest = FALSE, B = 500, beta.init = 'adaptive',
                                     beta.true = 0)
  })
}

save(obj, file = 'pval_glc_amd.RData')




proj_glc_amd = lasso.proj(x, y, standardize = TRUE, parallel = TRUE, family = 'binomial', ncores = 40)



save(proj_glc_amd, file = 'lasso_glc_amd.RData')



