.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_path")

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


# path.method <- args[6]

# beta.init <- args[7]

save.name <- args[5]

beta_i = as.numeric(args[6])
####


###################################################
###################################################
###################################################



if (beta_i ==0){
	rho = 0
}else if (beta_i ==1){
	rho = 0.5	
}else if (beta_i == 2){
	rho = 0.9
}else if (beta_i == 3){
	rho = 'equl'
}else if (beta_i == 4){
	rho = 'weak_equl'
}else{
	stop("wrong rho")
}

results = Path.Resample.Power(n = n, p = p, beta=c(0,rep(1,9),rep(0,90)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = 'L2.squared', path.method ='lars', beta.init = 'adaptive')

print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L2', "_", rho ,".RData")
print(f1)
save(results,file = f1)




results = Path.Resample.Power(n = n, p = p, beta=c(0,rep(1,9),rep(0,90)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = 'L1', path.method ='lars', beta.init = 'adaptive')
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L1',"_", rho, ".RData")
print(f1)
save(results,file = f1)



results = Path.Resample.Power(n = n, p = p, beta=c(0,rep(1,9),rep(0,90)), rho=rho, iter = iter, B = B, setting = 'dep', which.covariate = 1, betaNull = 0, multiTest = FALSE, parallel = TRUE, norm = 'L_inf', path.method ='lars', beta.init = 'adaptive')
print(results)
f1 = paste0("~/hdi_simu/results/",save.name,"_", 'L_inf',"_", rho, ".RData")
print(f1)
save(results,file = f1)




