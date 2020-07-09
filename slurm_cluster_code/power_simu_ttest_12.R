.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_simu")

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



bb = beta_i/10

for (rho in list(0, 0.5, 0.9, 'weak_equl','equl')){
results = TTest.Power(n = n, p = p, beta=c(bb,rep(1,2),rep(0, 9)), rho=rho, iter = iter,  setting = 'dep', which.covariate = 1, betaNull = 0)

print(mem_used())
f1 = paste0("~/hdi_simu/results/pc_ttest_",p,"_", 'rho',rho,'beta_',bb,".RData")
save(results,file = f1)

}
#results = desparse.Power(n = n, p = p, beta=c(bb,rep(1,9),rep(0,p-10)), rho=0.9, iter = iter,  setting = 'dep', which.covariate = 1, betaNull = 0)

#print(mem_used())
#f1 = paste0("~/hdi_path/results/SI/proj_AR09_p_",p,"_",bb,".RData")
#save(results,file = f1)


#results = desparse.Power(n = n, p = p, beta=c(bb,rep(1,9),rep(0,p-10)), rho="equl", iter = iter,  setting = 'dep', which.covariate = 1, betaNull = 0)

#print(mem_used())
#f1 = paste0("~/hdi_path/results/SI/proj_Eq_p_",p,"_",bb,".RData")
#save(results,file = f1)



#results = desparse.Power(n = n, p = p, beta=c(bb,rep(1,9),rep(0,p-10)), rho="weak_equl", iter = iter,  setting = 'dep', which.covariate = 1, betaNull = 0)

#print(mem_used())
#f1 = paste0("~/hdi_path/results/SI/proj_WkEq_p_",p,"_",bb,".RData")
#save(results,file = f1)







 


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




  

