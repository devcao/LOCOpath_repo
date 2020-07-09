.libPaths("~/R")
require(LOCOpath)
setwd("~/hdi_simu")
source("path_related.R")
source("path_resample.R")
source("path_power.R")



beta_true = c(c(0, rep(1, 9)), rep(0, 90))
n = 100; p = 100;

all_result_l1_0 = power.loco.test(n=n, p=p, beta=beta_true, rho=0, s=1,t=1,
                           iter = 500, B = 500,
                           whichCov=1, betaNULL=0, norm = 'L1')

save(all_result_l1_0, file = "all_result_l1_0.RData")

all_result_l2_0 = power.loco.test(n=n, p=p, beta=beta_true, rho=0, s=2,t=2,
                           iter = 500, B = 500,
                           whichCov=1, betaNULL=0, norm = 'L2.squared')
save(all_result_l2_0, file = "all_result_l2_0.RData")


all_result_l2_5 = power.loco.test(n=n, p=p, beta=beta_true, rho=0.5, s=2,t=2,
                           iter = 500, B = 500,
                           whichCov=1, betaNULL=0, norm = 'L2.squared')

save(all_result_l2_5, file = "all_result_l2_5.RData")


all_result_l1_5 = power.loco.test(n=n, p=p, beta=beta_true, rho=0.5, s=1,t=1,
                           iter = 500, B = 500,
                           whichCov=1, betaNULL=0, norm = 'L1')

save(all_result_l1_5, file = "all_result_l1_5.RData")


all_result_l2_9 = power.loco.test(n=n, p=p, beta=beta_true, rho=0.5, s=2,t=2,
                           iter = 500, B = 500,
                           whichCov=1, betaNULL=0, norm = 'L2.squared')

save(all_result_l2_9, file = "all_result_l2_9.RData")


all_result_l1_9 = power.loco.test(n=n, p=p, beta=beta_true, rho=0.5, s=1,t=1,
                           iter = 500, B = 500,
                           whichCov=1, betaNULL=0, norm = 'L1')
save(all_result_l1_9, file = "all_result_l1_9.RData")


