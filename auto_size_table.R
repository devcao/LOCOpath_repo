
"
\begin{table}[htbp]
\centering
\begin{tabular}{llllll}
\hline
\hline
Design & Method & $\alpha=0.20$ & $\alpha=0.10$ & $\alpha=0.05$& $\alpha=0.01$\\
\hline
\hline
$\Sigma = \mathbf{I}_p$ & $T_1(1,1)$ & 0.208 & 0.130 & 0.070 & 0.018 \\
& $T_1(2,2)$ & 0.242 & 0.124 & 0.066&  0.018 \\
& $T_1(\infty,\infty)$ & 0.268 & 0.126 & 0.068 & 0.024  \\
\hline
$\Sigma = (0.5^{|i-j|})_{1\leq i , j \leq p}$&  $T_1(1,1)$ & 0.230 & 0.128 & 0.078 & 0.018 \\
& $T_1(2,2)$ & 0.184 & 0.100 & 0.056 & 0.018 \\
& $T_1(\infty,\infty)$ & 0.186 & 0.100 & 0.050 & 0.014  \\
\hline
$\Sigma = (0.9^{|i-j|})_{1\leq i , j \leq p}$ & $T_1(1,1)$ & 0.210 & 0.134 & 0.076 & 0.020  \\
& $T_1(2,2)$ & 0.230 & 0.146 & 0.090 & 0.024 \\
& $T_1(\infty,\infty)$ & 0.212 & 0.102 & 0.046 & 0.000  \\
\hline
$\Sigma = (0.5^{\mathbf{1}(i \neq j)})_{1\leq i , j \leq p}$ & $T_1(1,1)$ & 0.184 & 0.102 & 0.058 &  0.012  \\
& $T_1(2,2)$ & 0.206 & 0.084 & 0.044 & 0.008 \\
& $T_1(\infty,\infty)$ & 0.176 & 0.096 & 0.058 & 0.016  \\
\hline
$\Sigma = (0.8^{\mathbf{1}(i \neq j)})_{1\leq i , j \leq p}$ & $T_1(1,1)$ & 0.214 & 0.142 & 0.098 & 0.042 \\
& $T_1(2,2)$ & 0.188 & 0.086 & 0.046 & 0.016 \\
& $T_1(\infty,\infty)$ & 0.148 & 0.054 & 0.026 & 0.004   \\
\hline
\hline
\end{tabular}
\caption{Multiple testing empirical size under different $\Sigma$ with $n=100$, $p=1000$.}
\label{tab:m_size_all_1000}
\end{table}
"

first_str = "\\begin{table}[htbp]
\\centering
\\begin{tabular}{llllll}
\\hline
\\hline
Design & Method & $\\alpha=0.20$ & $\\alpha=0.10$ & $\\alpha=0.05$& $\\alpha=0.01$\\\\
\\hline
\\hline
$\\Sigma = \\mathbf{I}_p$ & $T_1(1,1)$ &"

last_str = "\\hline
\\end{tabular}
\\caption{Multiple testing empirical size under different $\\Sigma$ with $n=100$, $p=1000$.}
\\label{tab:m_size_all_1000}
\\end{table}
"


auto_size_table = function(path_1, path_2, caption_name, label_name, path_2_change=FALSE,
                           rho, beta_sp = 0.1, recalculate = FALSE, comp_method=NULL
  
){
  
  
  if(beta_sp == 0.5){
    beta = seq(0,5,0.5)
    xlim_r = 5
  }else if(beta_sp == 0.1){
    beta = seq(0,1,0.1)
    xlim_r = 1
  }else{
    stop('wrong beta spacing')
  }
  # 
  
  size_table = list()
  
  i = 1
  for(rho in list(0, 0.5 ,0.9, 'weak_equl', 'equl')){
    
    power_00 = matrix(0,4,4)
    
    path1 = paste0(path_1,'_L1_rho_',rho, 'bb_',0,".RData")
    path2 = paste0(path_1,'_L2_rho_',rho, 'bb_',0,".RData")
    path3 = paste0(path_1,'_L_inf_rho_',rho, 'bb_',0,".RData")
    
    try({
      load(file = path1)
      power_00[1,] = results$power
     
      if(recalculate){
        power_00[1,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })
    try({
      load(file = path2)
      power_00[2,] = results$power
  
      if(recalculate){
        power_00[2,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })   
    try({
      load(file = path3)
      power_00[3,] = results$power
  
      if(recalculate){
        power_00[3,] =  apply(results$path.power, c(2,3), mean, na.rm=TRUE) 
      }
    })
    if (!missing('path_2')){
      
      path4 = paste0(path_2,rho, 'beta_',0,".RData")
      print(path4)
      
      try({
        load(file = path4)
        if(path_2_change){
          power_00[4,] = results$power  
        }else{
        power_00[4,] = results
        }
      })
    }
    size_table[[i]] = power_00
    i = i + 1
  }
 
  #print(size_table)
  str1 = paste0(first_str, 
    paste(size_table[[1]][1,], collapse = ' & '), '\\\\' ,'\n',
    '& $T_1(2,2)$ &', 
    paste(size_table[[1]][2,], collapse = ' & '), '\\\\' ,'\n',
    '& $T_1(\\infty,\\infty)$ &', 
    paste(size_table[[1]][3,], collapse = ' & '), '\\\\' ,'\n',
    '\\hline', '\n'
    )
  
  str_1 = paste0(first_str, 
                paste(size_table[[1]][1,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(2,2)$ &', 
                paste(size_table[[1]][2,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(\\infty,\\infty)$ &', 
                paste(size_table[[1]][3,], collapse = ' & '), '\\\\' ,'\n',
                '& ', comp_method, ' &', 
                paste(size_table[[1]][4,], collapse = ' & '), '\\\\' ,'\n',
                '\\hline', '\n'
  )
    
  str2 = paste0("$\\Sigma = (0.5^{|i-j|})_{1\\leq i , j \\leq p}$&  $T_1(1,1)$ &", 
                     paste(size_table[[2]][1,], collapse = ' & '), '\\\\' ,'\n',
                     '& $T_1(2,2)$ &', 
                     paste(size_table[[2]][2,], collapse = ' & '), '\\\\' ,'\n',
                     '& $T_1(\\infty,\\infty)$ &', 
                     paste(size_table[[2]][3,], collapse = ' & '), '\\\\' ,'\n',
                     '\\hline','\n'
  )
  
  str_2 = paste0("$\\Sigma = (0.5^{|i-j|})_{1\\leq i , j \\leq p}$&  $T_1(1,1)$ &", 
                paste(size_table[[2]][1,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(2,2)$ &', 
                paste(size_table[[2]][2,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(\\infty,\\infty)$ &', 
                paste(size_table[[2]][3,], collapse = ' & '), '\\\\' ,'\n',
                '& ', comp_method, ' &', 
                paste(size_table[[2]][4,], collapse = ' & '), '\\\\' ,'\n',
                '\\hline','\n'
  )
  
  str3 = paste0("$\\Sigma = (0.9^{|i-j|})_{1\\leq i , j \\leq p}$&  $T_1(1,1)$ &", 
                      paste(size_table[[3]][1,], collapse = ' & '), '\\\\' ,'\n',
                      '& $T_1(2,2)$ &', 
                      paste(size_table[[3]][2,], collapse = ' & '), '\\\\' ,'\n',
                      '& $T_1(\\infty,\\infty)$ &', 
                      paste(size_table[[3]][3,], collapse = ' & '), '\\\\' ,'\n',
                      '\\hline','\n'
  )
  
  str_3 = paste0("$\\Sigma = (0.9^{|i-j|})_{1\\leq i , j \\leq p}$&  $T_1(1,1)$ &", 
                paste(size_table[[3]][1,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(2,2)$ &', 
                paste(size_table[[3]][2,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(\\infty,\\infty)$ &', 
                paste(size_table[[3]][3,], collapse = ' & '), '\\\\' ,'\n',
                '& ', comp_method, ' &', 
                paste(size_table[[3]][4,], collapse = ' & '), '\\\\' ,'\n',
                '\\hline','\n'
  )
  
  str4 = paste0("$\\Sigma = (0.5^{\\mathbf{1}(i \\neq j)})_{1\\leq i , j \\leq p}$ & $T_1(1,1)$ &", 
                      paste(size_table[[4]][1,], collapse = ' & '), '\\\\' ,'\n',
                      '& $T_1(2,2)$ &', 
                      paste(size_table[[4]][2,], collapse = ' & '), '\\\\' ,'\n',
                      '& $T_1(\\infty,\\infty)$ &', 
                      paste(size_table[[4]][3,], collapse = ' & '), '\\\\' ,'\n',
                      '\\hline','\n'
  )
  
  str_4 = paste0("$\\Sigma = (0.5^{\\mathbf{1}(i \\neq j)})_{1\\leq i , j \\leq p}$ & $T_1(1,1)$ &", 
                paste(size_table[[4]][1,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(2,2)$ &', 
                paste(size_table[[4]][2,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(\\infty,\\infty)$ &', 
                paste(size_table[[4]][3,], collapse = ' & '), '\\\\' ,'\n',
                '& ', comp_method, ' &', 
                paste(size_table[[4]][4,], collapse = ' & '), '\\\\' ,'\n',
                '\\hline','\n'
  )
  
  str5 = paste0("$\\Sigma = (0.8^{\\mathbf{1}(i \\neq j)})_{1\\leq i , j \\leq p}$ & $T_1(1,1)$ &", 
                      paste(size_table[[5]][1,], collapse = ' & '), '\\\\' ,'\n',
                      '& $T_1(2,2)$ &', 
                      paste(size_table[[5]][2,], collapse = ' & '), '\\\\' ,'\n',
                      '& $T_1(\\infty,\\infty)$ &', 
                      paste(size_table[[5]][3,], collapse = ' & '), '\\\\' ,'\n',
                      '\\hline','\n'
  )
  
  str_5 = paste0("$\\Sigma = (0.8^{\\mathbf{1}(i \\neq j)})_{1\\leq i , j \\leq p}$ & $T_1(1,1)$ &", 
                paste(size_table[[5]][1,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(2,2)$ &', 
                paste(size_table[[5]][2,], collapse = ' & '), '\\\\' ,'\n',
                '& $T_1(\\infty,\\infty)$ &', 
                paste(size_table[[5]][3,], collapse = ' & '), '\\\\' ,'\n',
                '& ', comp_method, ' &', 
                paste(size_table[[5]][4,], collapse = ' & '), '\\\\' ,'\n',
                '\\hline','\n'
  )
  
  
  last_str = paste0("\\hline
  \\end{tabular}
  \\caption{", caption_name,
  "}
  \\label{tab:", label_name, 
  "}
  \\end{table}
  ")
  
  
  
  # $\Sigma = \mathbf{I}_p$ & $T_1(1,1)$ & 0.208 & 0.130 & 0.070 & 0.018 \\
  # & $T_1(2,2)$ & 0.242 & 0.124 & 0.066&  0.018 \\
  # & $T_1(\infty,\infty)$ & 0.268 & 0.126 & 0.068 & 0.024  \\
  # \hline
  if(missing('path_2')){
    return(cat(paste0(str1,str2,str3,str4,str5,last_str)))  
  }else{
    return(cat(paste0(str_1,str_2,str_3,str_4,str_5,last_str)))  
  }
  
   
}




auto_size_table(
  path_1 = 'log_new_new_pc_net_hdi_simu_80_dep',
  caption_name = 'Logistic new method p=80',  
  label_name = 'tab1',
  rho = rho,
  beta_sp = 0.5)

auto_size_table(
  path_1 = 'log_new_est_pc_net_hdi_simu_1000_dep',
  caption_name = 'Logistic new method p=1000',  
  label_name = 'tab1',
  rho = rho,
  beta_sp = 0.5)


auto_size_table(
  path_1 = 'log_pc_net_hdi_simu_1000_dep',
  caption_name = 'Logistic old method p=1000',  
  label_name = 'tab1',
  rho = rho,
  beta_sp = 0.5)



  auto_size_table(
    path_1 = 'log_new_est_pc_net_hdi_simu_80_dep',
    caption_name = 'Logistic new method p=80',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.5)

  
  
  auto_size_table(
    path_1 = 'log_new_new_pc_net_hdi_simu_80_dep',
    path_2 = 'log_pc_new_proj_80_rho',
    caption_name = 'Logistic old method p=80',  
    label_name = 'tab1',
    comp_method = 'De-sparsified',
    rho = rho,
    beta_sp = 0.5)
  
  auto_size_table(
    path_1 = 'poi_new_pc_net_hdi_simu_80_dep',
    path_2 = 'poi_pc_80ttest_80_rho', 
    caption_name = 'poi old method p=80',  
    label_name = 'tab1',
    comp_method = 'Wald',
    rho = rho,
    beta_sp = 0.1)
  
  

  
  
  ###
  new_net_power_plot(
    path_1 = 'poi_new_pc_net_hdi_simu_80_dep',
    path_2 = 'poi_pc_80ttest_80_rho', 
    file_name = 'poi_80_1',
    comp_method = 'Wald',
    beta_sp=0.1,
    cex_label = 1.4,
    poi_80=TRUE
  )
  
  ###
  
  
  
  auto_size_table(
    path_1 = 'log_new_est_pc_net_hdi_simu_12_dep',
    path_2 = 'log_pc_12ttest_12_rho', # 
    caption_name = 'Logistic new method p=12',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1,
    comp_method = "Wald")
  
  
  
  auto_size_table(
    path_1 = 'log_NNnew_pc_net_hdi_simu_12_dep',
    path_2 = 'log_pc_12ttest_12_rho', # 
    caption_name = 'Logistic old method p=12',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1,
    comp_method = "Wald")
  
  
  auto_size_table(
    path_1 = 'log_NNnew_pc_net_hdi_simu_12_dep',
    path_2 = 'log_pc_12ttest_12_rho', # 
    caption_name = 'Logistic old method p=12',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1,
    comp_method = "Wald")
  
  
  
  
  
  auto_size_table(
    path_1 = 'log_test_multi_80_pc_hdi_simu_dep',
    caption_name = 'loistic multi p=80',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  auto_size_table(
    path_1 = 'log_glasso_multi_beta05_80_pc_hdi_simu_dep',
    caption_name = 'loistic glasso multi p=80',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  auto_size_table(
    path_1 = 'log_mmulti_pc_hdi_simu_dep',
    caption_name = 'loistic multi p=1000',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  auto_size_table(
    path_1 = 'log_glasso_multi_beta05_1k_pc_hdi_simu_dep',
    caption_name = 'loistic glasso multi p=1000',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  
  auto_size_table(
    path_1 = 'non0_log_pc_net_hdi_simu_80_dep',
    caption_name = 'loistic beta_1=1 p=80',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1, recalculate = TRUE)
  
  
  
  auto_size_table(
    path_1 = 'non0_log_pc_net_hdi_simu_12_dep',
    caption_name = 'loistic beta_1=1 p=12',
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  
  
  
  
  auto_size_table(
    path_1 = 'poi_3signal_new_pc_net_hdi_simu_1000_dep',
    caption_name = 'Poisson p=1000',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  
  
  auto_size_table(
    path_1 = 'pc_right_eq_hdi_simu_80_dep',
    path_2 = 'eq_pc_ttest_80_rho', # 
    caption_name = 'linear beta1=beta2 p=80',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1,
    comp_method = "T-test")
  
  
  
  auto_size_table(
    path_1 = 'pc_right_eq_hdi_simu_12_dep',
    path_2 = 'eq_pc_ttest_12_rho', # 
    caption_name = 'linear beta1=beta2 p=12',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1,
    comp_method = "T-test")
  
  
  
  
  
  auto_size_table(
    path_1 = 'pc_right_eq_hdi_simu_1000_dep',
    caption_name = 'linear beta1=beta2 p=1000',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  auto_size_table(
    path_1 = 'log_size_new_pc_net_hdi_simu_1000_dep',
    caption_name = 'Logistic p=1000 B=10k',
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  auto_size_table(
    path_1 = 'size_log_net_1000',
    caption_name = 'logistic p=1000 new alpha',
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  auto_size_table(
    path_1 = 'log_new_pc_net_hdi_simu_1000_dep',
    path_2 = 'log_pc_new_proj_1000_rho',
    caption_name = 'log_1000',
    label_name = 'tab_log_1000',
    rho = rho,
    beta_sp = 0.1,
    path_2_change = TRUE)
  
  
  
  #### ??? #### p=1000 poisson ??
  #### p = 80 poisson
  #### poi_pc_80ttest need to label how many time failed
  
  
  auto_size_table(
    path_1 = 'log_new_new_pc_net_hdi_simu_80_dep',
    caption_name = 'log_80_new',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  auto_size_table(
    path_1 = 'log_size_3inter_new_pc_net_hdi_simu_80_dep',
    caption_name = 'log_80_B_10k',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  
  auto_size_table(
    path_1 = 'poi_3signal_multi_pc_hdi_simu_dep',
    caption_name = 'poi_1k_multi',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  
  
  auto_size_table(
    path_1 = 'poi_80_multi_pc_hdi_simu_dep',
    caption_name = 'poi_80_multi',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  
  
  
  auto_size_table(
    path_1 = 'poi_new_pc_net_hdi_simu_80_dep',
    path_2 = 'poi_pc_80ttest_80_rho', # 
    caption_name = 'poisson p=80',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1,
    comp_method = "Wald")
  
  
  auto_size_table(
    path_1 = 'poi_new_pc_net_hdi_simu_12_dep',
    path_2 = 'poi_pc_12ttest_12_rho', # 
    caption_name = 'poisson p=12',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1,
    comp_method = "Wald")
  
  
  
  auto_size_table(
    path_1 = 'poi_80_multi_pc_hdi_simu_dep',
    caption_name = 'poisson multi p=80',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1)
  
  
  auto_size_table(
    path_1 = 'log_new_alter_multi_80_pc_hdi_simu_dep', 
    path_2 = 'log_new_alter_glasso_multi_beta05_80_pc_hdi_simu_dep', 
    caption_name = 'log comp p=80',  
    label_name = 'tab1',
    rho = rho,
    beta_sp = 0.1,
    comp_method = "G-lasso")
  
  
  


