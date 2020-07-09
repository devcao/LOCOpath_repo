.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_simu")

require(mvnfast)
source("path_related.R")
source("path_power.R")
source("path_resample.R")

source("NetTS.R")
source("NetResampleTS.R")
source("NetResampleLogisticTS.R")
source("Logistic_Enet.R")

require(LOCOpath)
require(gglasso)
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
#for (rho in list(0, 0.5, 0.9, 'equl', 'weak_equl')){
rho = 0.9

results = Net.Resample.Logistic.Multi.Power(n = n, p = p, beta=c(0+bb,rep(4,2),bb,bb,rep(0,75)), 
                              rho=rho, iter = iter, B = B,
                              which.covariate = list(c(1,4,5)), 
                              betaNull = list(c(0,0,0)), 
                              parallel = TRUE, 
                              norm = 'L2.squared',
                              beta.init = 'adaptive')


print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L2', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

##############

results = Net.Resample.Logistic.Multi.Power(n = n, p = p, beta=c(0+bb,rep(4,2),bb,bb,rep(0,75)), 
                              rho=rho, iter = iter, B = B, 
                              which.covariate = list(c(1,4,5)), 
                              betaNull = list(c(0,0,0)), 
                              parallel = TRUE, 
                              norm = 'L1',  
                              beta.init = 'adaptive')
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L1', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

##############

results = Net.Resample.Logistic.Multi.Power(n = n, p = p, beta=c(0+bb,rep(4,2),bb,bb,rep(0,75)), 
                              rho=rho, iter = iter, B = B, 
                              which.covariate = list(c(1,4,5)), 
                              betaNull = list(c(0,0,0)), 
                              parallel = TRUE, 
                              norm = 'L_inf', 
                              beta.init = 'adaptive')
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L_inf', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

#}


