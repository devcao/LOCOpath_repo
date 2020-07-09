.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_simu")
source("path_power.R")
require(mvnfast)


source('~/hdi_simu/NetResampleLogisticTS.R')
source('~/hdi_simu/NetTS.R')
source('~/hdi_simu/NetResampleTS.R')

require(LOCOpath)
# source("~/hdi_path/bin/pathwise_power.R")

print(depenDesign)
#######################################################################################################################################
# Set simulation parameters (to be done with command-line arguments)

# Execute this from within the directory containing this R script:
############################################################################


options(echo=TRUE)

args <- commandArgs(trailingOnly = TRUE)

print(args)

# args <- c("1000","2","1","20","3","5","3",".92",".96",".95",".98","1")

n <- as.numeric(args[1])

p <- as.numeric(args[2])

iter <- as.numeric(args[3])  # goes from 1 to 12

B <- as.numeric(args[4])

save.name <- args[5]

beta_i = as.numeric(args[6])
###################################################
###################################################
###################################################

bb = 0 # beta_i / 2
#for (rho in list(0, 0.5, 0.9, 'equl', 'weak_equl')){


if(beta_i == 0){
  norm = 'L2.squared'
}else if(beta_i == 1){
  norm = 'L1'
}else if (beta_i == 2){
  norm = 'L_inf'
}
print(norm)

rho = 0
results = Net.Resample.Logistic.Power(n = n, p = p, beta=c(0+bb,rep(1,9),rep(0,990)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = norm, beta.init = 'adaptive', beta.init.alpha=0.7, alpha=0.2)

print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", norm, "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

##############
rho = 0.5
results = Net.Resample.Logistic.Power(n = n, p = p, beta=c(0+bb,rep(1,9),rep(0,990)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = norm, beta.init = 'adaptive', beta.init.alpha=0.6, alpha=0.8)
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", norm, "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

##############
rho = 0.9
results = Net.Resample.Logistic.Power(n = n, p = p, beta=c(0+bb,rep(1,9),rep(0,990)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = norm, beta.init = 'adaptive', beta.init.alpha=0.4)
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", norm, "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

#}


