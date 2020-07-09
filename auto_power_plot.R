
require(latex2exp)

## generate power curve automatically
new_net_power_plot = function(path_1, path_2, file_name, beta_sp = 0.1, cex_label=1.4, 
                              recalculate = FALSE, comp_method=NULL,
                              path_2_alter = FALSE, beta_imp, poi_80=FALSE){
  # Args:
  # path_1: the main part of the name of your output R data
  # path_2: can be ignored. If specified, path_2 is the method you want to compare the power.
  # file_name: name of the output pdf
  # beta_sp: the space of the beta sequence
  # cex_label: cex of the legend
  # recalculate: if true, will recalculate the power and omit NA
  # comp_method: the name of the method you want to campare
  # path_2_alter: if true, extract data from path_2 as ..$power, otherwise extract the whole data
  # beta_imp: i forgot, use default please
  # poi_80: fix some small issue for poisson regression, you can ignore
  # Return: NULL, will output a pdf power curve to your current directory
  if(beta_sp == 0.5){
    beta = seq(0,5,0.5)
    xlim_r = 5
  }else if(beta_sp == 0.1){
    beta = seq(0,1,0.1)
    xlim_r = 1
  }else if(beta_sp == 0.02){
    beta = seq(0,0.2,0.02)
    xlim_r = 0.2
  }else{
    stop('wrong beta spacing')
  }
  # 
  
  
  rho = 0
  power_rho0_1 = matrix(0,11,4)
  power_rho0_2 = matrix(0,11,4)
  power_rho0_3 = matrix(0,11,4)
  power_rho0_4 = matrix(0,11,4)
  
  for (i in 1:11){
    path1 = paste0(path_1,'_L1_rho_',rho, 'bb_',beta[i],".RData")
    path2 = paste0(path_1,'_L2_rho_',rho, 'bb_',beta[i],".RData")
    path3 = paste0(path_1,'_L_inf_rho_',rho, 'bb_',beta[i],".RData")
    
    try({
      load(file = path1)
      power_rho0_1[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho0_1[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })
    try({
      load(file = path2)
      power_rho0_2[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho0_2[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })   
    try({
      load(file = path3)
      power_rho0_3[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho0_3[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })
    if (!missing('path_2')){
      
      path4 = paste0(path_2,rho, 'beta_',beta[i],".RData")
      print(path4)
      
      try({
        if(path_2_alter){
          load(file = path4)
          power_rho0_4[i,] = results$power
        }else{
        load(file = path4)
        power_rho0_4[i,] = results
        print(results)
        }
      })
    }
  }
  
  rho = 0.5
  power_rho5_1 = matrix(0,11,4)
  power_rho5_2 = matrix(0,11,4)
  power_rho5_3 = matrix(0,11,4)
  power_rho5_4 = matrix(0,11,4)
  
  for (i in 1:11){
    path1 = paste0(path_1,'_L1_rho_',rho, 'bb_',beta[i],".RData")
    path2 = paste0(path_1,'_L2_rho_',rho, 'bb_',beta[i],".RData")
    path3 = paste0(path_1,'_L_inf_rho_',rho, 'bb_',beta[i],".RData")
    
    try({
      load(file = path1)
      power_rho5_1[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho5_1[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })
    try({
      load(file = path2)
      power_rho5_2[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho5_2[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })   
    try({
      load(file = path3)
      power_rho5_3[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho5_3[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })
    if (!missing('path_2')){
      
      path4 = paste0(path_2,rho, 'beta_',beta[i],".RData")
      print(path4)
      
      try({
        if(path_2_alter){
          load(file = path4)
          power_rho5_4[i,] = results$power
        }else{
          load(file = path4)
          power_rho5_4[i,] = results
          print(results)
        }
      })
    }
  }
  
  rho = 0.9
  power_rho9_1 = matrix(0,11,4)
  power_rho9_2 = matrix(0,11,4)
  power_rho9_3 = matrix(0,11,4)
  power_rho9_4 = matrix(0,11,4)
  
  for (i in 1:11){
    path1 = paste0(path_1,'_L1_rho_',rho, 'bb_',beta[i],".RData")
    path2 = paste0(path_1,'_L2_rho_',rho, 'bb_',beta[i],".RData")
    path3 = paste0(path_1,'_L_inf_rho_',rho, 'bb_',beta[i],".RData")
    
    try({
      load(file = path1)
      power_rho9_1[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho9_1[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })
    try({
      load(file = path2)
      power_rho9_2[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho9_2[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })   
    try({
      load(file = path3)
      power_rho9_3[i,] = results$power
      if(i==1) print(results$power)
      if(recalculate){
        power_rho9_3[i,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })
    if (!missing('path_2')){
      
      path4 = paste0(path_2,rho, 'beta_',beta[i],".RData")
      print(path4)
      
      try({
        if(path_2_alter){
          load(file = path4)
          power_rho9_4[i,] = results$power
        }else{
          load(file = path4)
          power_rho9_4[i,] = results
          print(results)
        }
      })
    }
  }
  
  xlim_l = 0
  
  
  
  cex = cex_label
  file_names = paste0(file_name,".pdf") 
  pdf(file_names, width = 5, height = 9.61)
  
  par(mfrow=c(3,1),mar=c(0,0,0,0),oma = c(5.1, 4.1, 4.1, 2.1))
  
  plot(NA,xlim=c(xlim_l,xlim_r),ylim=c(0,1),xaxt="n")
  lines(power_rho0_1[,3]~beta, pch = 19, type = 'b', lwd = 1.5, col = 1)
  lines(power_rho0_2[,3]~beta, pch = 5, type = 'b', lwd = 1.5, col = 1)
  lines(power_rho0_3[,3]~beta, pch = 3, type = 'b', lwd = 1.5, col = 1)
  if (!missing('path_2')){
    lines(power_rho0_4[,3]~beta, pch = 1, type = 'b', lwd = 1.5, col = 2)
  }
  abline(h = 0.05, lty = 2)
  
  if(beta_sp==0.5){
  legend( x = -0.2,#grconvertX(10.5, from = 'nfc', to = 'user'),
          y = 1,#grconvertX(0.6, from = 'nfc', to = 'user'),
          legend=TeX('$\\Sigma = \\mathbf{I}_p$ '),
          bty="n", cex = cex)
  }else if(beta_sp==0.1){
    if (poi_80==TRUE){
      legend( x = 0.7,#$grconvertX(0, from = 'nfc', to = 'user'),
              y = 0.9,#grconvertX(1, from = 'nfc', to = 'user'),
              legend=TeX('$\\Sigma = \\mathbf{I}_p$ '),
              bty="n", cex = cex)
    }else{
    legend( x = -0.05,#grconvertX(10.5, from = 'nfc', to = 'user'),
            y = 1,#grconvertX(0.6, from = 'nfc', to = 'user'),
            legend=TeX('$\\Sigma = \\mathbf{I}_p$ '),
            bty="n", cex = cex)  
    }
  }else if(beta_sp==0.02){
    legend( x = -0.01,#grconvertX(10.5, from = 'nfc', to = 'user'),
            y = 1,#grconvertX(0.6, from = 'nfc', to = 'user'),
            legend=TeX('$\\Sigma = \\mathbf{I}_p$ '),
            bty="n", cex = cex)  
  }
  
  plot(NA,xlim=c(xlim_l,xlim_r),ylim=c(0,1),xaxt="n")
  lines(power_rho5_1[,3]~beta, pch = 19, type = 'b', lwd = 1.5, col = 1)
  lines(power_rho5_2[,3]~beta, pch = 5, type = 'b', lwd = 1.5, col = 1)
  lines(power_rho5_3[,3]~beta, pch = 3, type = 'b', lwd = 1.5, col = 1)
  if (!missing('path_2')){
    lines(power_rho5_4[,3]~beta, pch = 1, type = 'b', lwd = 1.5, col = 2)
  }
  abline(h = 0.05, lty = 2)
  
  if(beta_sp==0.5){
  legend( x = -0.2,#$grconvertX(0, from = 'nfc', to = 'user'),
          y = 1,#grconvertX(1, from = 'nfc', to = 'user'),
          legend=TeX('$\\Sigma = (0.5^{|i-j|})_{1\\leq i , j \\leq p}$ '),
          bty="n", cex = cex)
  }else if(beta_sp==0.1){
    if (poi_80==TRUE){
      legend( x = 0.6,#$grconvertX(0, from = 'nfc', to = 'user'),
              y = 0.9,#grconvertX(1, from = 'nfc', to = 'user'),
              legend=TeX('$\\Sigma = (0.5^{|i-j|})_{1\\leq i , j \\leq p}$ '),
              bty="n", cex = cex)
    }else{
    legend( x = -0.05,#$grconvertX(0, from = 'nfc', to = 'user'),
            y = 1,#grconvertX(1, from = 'nfc', to = 'user'),
            legend=TeX('$\\Sigma = (0.5^{|i-j|})_{1\\leq i , j \\leq p}$ '),
            bty="n", cex = cex)
    }
  }else if(beta_sp==0.02){
    legend( x = -0.01,#grconvertX(10.5, from = 'nfc', to = 'user'),
            y = 1,#grconvertX(0.6, from = 'nfc', to = 'user'),
            legend=TeX('$\\Sigma = (0.5^{|i-j|})_{1\\leq i , j \\leq p}$ '),
            bty="n", cex = cex)  
  }
  
  
  
  plot(NA,xlim=c(xlim_l,xlim_r),ylim=c(0,1))
  lines(power_rho9_1[,3]~beta, pch = 19, type = 'b', lwd = 1.5, col = 1)
  lines(power_rho9_2[,3]~beta, pch = 5, type = 'b', lwd = 1.5, col = 1)
  lines(power_rho9_3[,3]~beta, pch = 3, type = 'b', lwd = 1.5, col = 1)
  if (!missing('path_2')){
    lines(power_rho9_4[,3]~beta, pch = 1, type = 'b', lwd = 1.5, col = 2)
  }
  abline(h = 0.05, lty = 2)
  
  if(beta_sp==0.5){
    legend( x = -0.2,#$grconvertX(0, from = 'nfc', to = 'user'),
            y = 1,#grconvertX(1, from = 'nfc', to = 'user'),
            legend=TeX('$\\Sigma = (0.9^{|i-j|})_{1\\leq i , j \\leq p}$ '),
            bty="n", cex = cex)
  }else if(beta_sp==0.1){
    if (poi_80==TRUE){
      legend( x = 0.6,#$grconvertX(0, from = 'nfc', to = 'user'),
              y = 0.9,#grconvertX(1, from = 'nfc', to = 'user'),
              legend=TeX('$\\Sigma = (0.9^{|i-j|})_{1\\leq i , j \\leq p}$ '),
              bty="n", cex = cex)
    }else{
    legend( x = -0.05,#$grconvertX(0, from = 'nfc', to = 'user'),
            y = 1,#grconvertX(1, from = 'nfc', to = 'user'),
            legend=TeX('$\\Sigma = (0.9^{|i-j|})_{1\\leq i , j \\leq p}$ '),
            bty="n", cex = cex)
    }
  }else if(beta_sp==0.02){
    legend( x = -0.01,#grconvertX(10.5, from = 'nfc', to = 'user'),
            y = 1,#grconvertX(0.6, from = 'nfc', to = 'user'),
            legend=TeX('$\\Sigma = (0.9^{|i-j|})_{1\\leq i , j \\leq p}$ '),
            bty="n", cex = cex)  
  }
  
  
  ### may need to change this 
  if (!missing('beta_imp')){
    mtext(side=1, outer=TRUE, expression(paste(beta[1],'-1')),line=3)
    mtext(side=2, outer=TRUE, "Empirical power",line=2.5)
  }else{
    mtext(side=1, outer=TRUE, expression(beta[1]),line=3)
    mtext(side=2, outer=TRUE, "Empirical power",line=2.5)
  }
  
  
  if(beta_sp==0.1){
    legend_x = grconvertX(0.1, from = 'nfc', to = 'user')#0.15#-4.4
    legend_y = grconvertX(3.15, from = 'nfc', to = 'user')
  }else if(beta_sp == 0.5){
    legend_x = grconvertX(0.1, from = 'nfc', to = 'user')#1.1#-4.4
    legend_y = grconvertX(0.66, from = 'nfc', to = 'user')
  }else if(beta_sp == 0.02){
    legend_x = grconvertX(0.1, from = 'nfc', to = 'user')#1.1#-4.4
    legend_y = grconvertX(15.5, from = 'nfc', to = 'user')
  }else {
    stop('wrong beta spacing')
  }
  
  if(!missing('path_2')){
    legend( x = legend_x,
            y = legend_y,
            legend=c( expression(T[1](1,1)),
                      expression(T[1](2,2)),
                      expression(T[1](infinity,infinity)),
                      comp_method),
            pch=c(19,5,3,1),
            col = c(1,1,1,2),
            lwd = 1,
            bty="n", xpd = NA, horiz = TRUE, cex=1.2)
  }else{
    legend( x = legend_x,#grconvertX(-0.1, from = 'nfc', to = 'user'),
            y = legend_y,#grconvertX(0.5, from = 'nfc', to = 'user'),
            legend=c( expression(T[1](1,1)),
                      expression(T[1](2,2)),
                      expression(T[1](infinity,infinity))
            ),
            pch=c(19,5,3),
            col = c(1,1,1),
            lwd = 1,
            bty="n", xpd = NA, horiz = TRUE, cex=1.4)
  }
  
  
  
  dev.off()
}
#####################################################
    
##############################
### some examples
##############################    
    
new_net_power_plot(
  path_1 = 'poi_poi_3signal_new_pc_net_hdi_simu_1000_dep',
  file_name = 'poi_neww_1000',
  beta_sp=0.02
)


new_net_power_plot(
  path_1 = 'poi_poi_3signal_multi_pc_hdi_simu_dep',
  file_name = 'poi_new_multi_1000',
  beta_sp=0.02
)


new_net_power_plot(
  path_1 = 'log_mmulti_pc_hdi_simu_dep',
  file_name = 'log_multi_1000',
  beta_sp=0.5
)
#log_mmulti_pc_hdi_simu_dep

new_net_power_plot(
  path_1 = 'log_test_multi_80_pc_hdi_simu_dep',
  file_name = 'log_multi_80',
  beta_sp=0.5
)

new_net_power_plot(
  path_1 = 'non0_log_pc_net_hdi_simu_80_dep',
  
  file_name = 'non0_log_80',

  beta_sp=0.5,
  recalculate = TRUE,
  beta_imp=seq(1,6,0.5)
)



new_net_power_plot(
  path_1 = 'log_new_pc_net_hdi_simu_1000_dep',
  path_2 = 'log_pc_new_proj_1000_rho',
  file_name = 'log_1000',
  comp_method = 'De-sparsified',
  beta_sp=0.5
)



new_net_power_plot(
  path_1 = 'log_new_pc_net_hdi_simu_1000_dep',
  path_2 = 'log_pc_new_proj_1000_rho',
  file_name = 'log_1000',
  comp_method = 'De-sparsified',
  beta_sp=0.5
)



new_net_power_plot(
  path_1 = 'log_new_new_pc_net_hdi_simu_80_dep',
  path_2 = 'log_pc_new_proj_80_rho',
  file_name = 'log_80',
  comp_method = 'De-sparsified',
  beta_sp=0.5,
  path_2_alter = TRUE
)


log_mmulti_pc_hdi_simu_dep

new_net_power_plot(
  path_1 = 'poi_3signal_new_pc_net_hdi_simu_1000_dep',
  file_name = 'poi_1000',
  beta_sp=0.1
)

new_net_power_plot(
  path_1 = 'poi_3signal_multi_pc_hdi_simu_dep',
  file_name = 'poi_multi_1000',
  beta_sp=0.1
)




new_net_power_plot(
  path_1 = 'pc_right_eq_hdi_simu_1000_dep',
  file_name = 'eq_1000',
  beta_sp=0.1,
  cex_label = 1.4,
  beta_imp=TRUE
)



new_net_power_plot(
  path_1 = 'pc_right_eq_hdi_simu_80_dep',
  path_2 = 'eq_pc_ttest_80_rho', # 
  file_name = 'eq_80',
  comp_method = 'T-test',
  beta_sp=0.1,
  cex_label = 1.4,
  beta_imp=TRUE
)



new_net_power_plot(
  path_1 = 'poi_new_pc_net_hdi_simu_80_dep',
  path_2 = 'poi_pc_80ttest_80_rho', 
  file_name = 'poi_80_1',
  comp_method = 'Wald',
  beta_sp=0.1,
  cex_label = 1.4,
  poi_80=TRUE
)



new_net_power_plot(
  path_1 = 'poi_80_multi_pc_hdi_simu_dep',
  file_name = 'poi_multi_80',
  beta_sp=0.1,
  cex_label = 1.4,
  poi_80=TRUE
)




#path_1 = 'poi_new_pc_net_hdi_simu_80_dep'
#path_2 = 'poi_pc_80ttest_80_rho'
#beta_sp=0.1

