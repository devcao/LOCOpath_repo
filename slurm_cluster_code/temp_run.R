.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_simu")

require(mvnfast)
source("~/hdi_simu/compare_power.R")


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

beta_i = as.numeric(args[4])
####

##### if not running on a cluster
#n = 100
#p = 12
#iter = 500
#B = 500
#####


###################################################
###################################################
###################################################



bb = 5
rho=0.5

results = desparse.Logistic.Power(n = n, p = p, beta=c(bb,rep(1,9),rep(0, 990)), rho=rho, iter = iter, which.covariate = 1, betaNull = 0)

print(mem_used())
f1 = paste0("~/hdi_simu/results/log_pc_new_proj_",p,"_", 'rho',rho,'beta_',bb,".RData")
save(results,file = f1)





  

