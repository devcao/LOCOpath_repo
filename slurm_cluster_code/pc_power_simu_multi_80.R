.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_simu")
source("path_related.R")
source("path_power.R")
source("path_resample.R")
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

bb = beta_i / 10
for (rho in list(0, 0.5, 0.9, 'equl', 'weak_equl')){

results = Path.Resample.Power(n = n, p = p, beta=c(1+bb,rep(1,2),rep(0,77)), 
                              rho=rho, iter = iter, B = B, setting = 'dep',
                              which.covariate = list(c(1,4,5)), 
                              betaNull = list(c(1,0,0)), 
                              multiTest = TRUE, parallel = TRUE, 
                              norm = 'L2.squared', path.method ='lars', 
                              beta.init = 'adaptive')


print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_",'p_80', 'L2', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

##############

results = Path.Resample.Power(n = n, p = p, beta=c(1+bb,rep(1,2),rep(0,77)), 
                              rho=rho, iter = iter, B = B, setting = 'dep', 
                              which.covariate = list(c(1,4,5)), 
                              betaNull = list(c(1,0,0)), 
                              multiTest = TRUE, parallel = TRUE, 
                              norm = 'L1', path.method ='lars', 
                              beta.init = 'adaptive')
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_",'p_80' ,'L1', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

##############

results = Path.Resample.Power(n = n, p = p, beta=c(1+bb,rep(1,2),rep(0,77)), 
                              rho=rho, iter = iter, B = B, setting = 'dep', 
                              which.covariate = list(c(1,4,5)), 
                              betaNull = list(c(1,0,0)), 
                              multiTest = TRUE, parallel = TRUE, 
                              norm = 'L_inf', path.method ='lars', 
                              beta.init = 'adaptive')
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_",'p_80', 'L_inf', "_rho_", rho, 'bb_', bb, ".RData")
print(f1)
save(results,file = f1)

}


