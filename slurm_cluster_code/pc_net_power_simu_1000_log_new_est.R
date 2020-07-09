.libPaths(new="~/R")

rm(list=ls())
require(mvnfast)
setwd("~/hdi_simu")
source("path_power.R")

source('~/hdi_simu/NetTS.R')
source('~/hdi_simu/NetResampleTS.R')

source('~/hdi_simu/NetResampleLogisticTS.R')

require(LOCOpath)
# source("~/hdi_path/bin/pathwise_power.R")


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

bb = beta_i / 2
for (rho in list(0, 0.5, 0.9, 'equl', 'weak_equl')){

results = Net.Resample.Logistic.Power(n = n, p = p, beta=c(0+bb,rep(1,9),rep(0,990)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = 'L2.squared', beta.init = 'adaptive', beta.null.estimate=TRUE)

print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L2', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

##############

results = Net.Resample.Logistic.Power(n = n, p = p, beta=c(0+bb,rep(1,9),rep(0,990)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = 'L1', beta.init = 'adaptive', beta.null.estimate=TRUE)
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L1', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

##############

results = Net.Resample.Logistic.Power(n = n, p = p, beta=c(0+bb,rep(1,9),rep(0,990)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = 'L_inf', beta.init = 'adaptive', beta.null.estimate=TRUE)
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L_inf', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

}


