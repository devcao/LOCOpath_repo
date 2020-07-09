rm(list = ls())
.libPaths('~/R')
setwd("~/hdi_simu")
require(hdi)

source("path_related.R")
source("path_resample.R")
source("path_power.R")

require(LOCOpath)

data(riboflavin)

load("screened_index_lars.RData")

util_fun = function(...){
	return(loco_resample(...)$pval)
}
p = dim(riboflavin$x)[2]
obj = list()

x = riboflavin$x[, index]

y = riboflavin$y

for (i in 1:length(index)){
	obj[[i]]=loco_resample(path_type = "lars", x = x, y = y, 
              whichCov = i, betaNULL = 0, 
              s = 1, t = 1,
		enet.control = list(),
             # enet.control = list(nlambda = 10, alpha = 1, family = "gaussian"),
              B = 10000, n_threads = -1,
              plot_null = FALSE)

}
#print(pval)
save(obj, file = 'LL_lars_pval.RData')
