.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_simu")
require(LOCOpath)

source("path_related.R"); source("path_resample.R"); source("path_power.R")

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
####


###################################################
###################################################
###################################################



bb = beta_i/20


results = Path.Resample.Power(n = n, p = p, beta=beta, rho="weak_equl", iter = iter, B = B, setting = 'dep', which.covariate = list(1:10), betaNull = list(rep(1,10)), multiTest = TRUE, parallel = TRUE, norm = norm, path.method = path.method, beta.init = beta.init)

print(mem_used())
f1 = paste0("~/hdi_path/results/quick/",norm,"_",p,"_","M1_WkEq_",save.name,bb,".RData")
save(results,file = f1)



results = Path.Resample.Power(n = n, p = p, beta=c(rep(1,10),rep(0,p-10)), rho="weak_equl", iter = iter, B = B, setting = 'dep', which.covariate = list(11:p), betaNull = list(rep(0,p-10)), multiTest = TRUE, parallel = TRUE, norm = norm, path.method = path.method, beta.init = beta.init)

print(mem_used())
f1 = paste0("~/hdi_path/results/quick/",norm,"_",p,"_","M0_WkEq_",save.name,bb,".RData")
save(results,file = f1)






#p0=runif(10,0,2)
#results = Path.Resample.Power(n = n, p = p, beta=c(rep(1,10),bb,bb,rep(0,988)), rho=0, multiTest = TRUE, iter = iter, B = B, setting = 'dep', which.covariate = list(c(1,2,11,12)), betaNull = list(c(1,1,0,0)), parallel = TRUE, norm = norm, path.method = path.method, beta.init = beta.init)


#print(mem_used())
#f3 = paste0("~/hdi_path/results/L2.sq/Multiple_Exp_",save.name,bb,".RData")
#save(results,file = f3)


#p0 = runif(10,0,2)
#results = Path.Resample.Power(n = n, p = p, beta=c(rep(1,10),bb,bb,rep(0,988)), rho="equl", iter = iter, B = B,multiTest = TRUE,  setting = 'dep', which.covariate = list(c(1,2,11,12)), betaNull = list(c(1,1,0,0)), parallel = TRUE, norm = norm, path.method = path.method, beta.init = beta.init)


#print(mem_used())
#f4 = paste0("~/hdi_path/results/L2.sq/Multiple_Equal_",save.name,bb,".RData")
#save(results,file = f4)




  

