.libPaths(new="~/R")

rm(list=ls())

setwd("~/hdi_simu")

#require(LOCOpath)
source('graphLASSO.R')


n = 100
p = 50

results = simu_graph_screen(n = n, p = p, type = 'A', Iter = 250)
save(results, file = 'type_A_n_100_p_50.RData')




