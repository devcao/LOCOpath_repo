
# This file contains basic code for graphical LASSO 

###### Loading some useful package, not all required ########
require(glasso)
require(CVglasso)
require(LOCOpath)
require(huge)
#source('graphLASSO.R')
source('NetTS.R')
require(igraph)
#######################################

########## A mofifiled version of the huge.plot function in huge package ##############
huge.plot_v1 = function (G, epsflag = FALSE, graph.name = "default", cur.num = 1, 
                         location = NULL, ...) 
{
  # A mofifiled version of the huge.plot function in huge package
  # All arguments the same as huge.plot, check huge.plot for more details
  gcinfo(FALSE)
  if (missing(location)) 
    location = tempdir()
  oldlocation = getwd()
  a
  setwd(location)
  g = graph.adjacency(as.matrix(G != 0), mode = "undirected", 
                      diag = FALSE)
  layout.grid = layout.fruchterman.reingold(g)
  if (epsflag == TRUE) 
    postscript(paste(paste(graph.name, cur.num, sep = ""), 
                     "eps", sep = "."), width = 8, height = 8)
  
  plot(g, layout = layout.grid, edge.color = "gray50", vertex.color = "red", 
       vertex.size = 2, vertex.label = NA, ...)
  rm(g, location)
  gc()
  if (epsflag == TRUE) 
    dev.off()
  setwd(oldlocation)
}



### 2 functoins to make the imgage function in R generate correct matrix plot ###
rotate <- function(x) t(apply(x, 2, rev))  # working function#\
image_v1 = function(mat,...){image(rotate(mat),...)} # new version of the image_funcion



###### graphical model variable screen #######
simu_graph_screen = function(n = 100, p = 50, type = 'A', Iter = 250){
  # graphical model variable screen, output will be used to generete ROC curve
  # Args:
  # n, p; int, sample size and dimension of matrix
  # type: allow 2 options ('A', 'C'), 2 types of matices
  # Iter: int, number if iterations, use defualt 50
  # Returns:
  #  A List of: Theta_list: a list of all randomly generated precision matrix
  #             glasso_results: a list of glasso estimates
  #             LOCO_results: a list of calulated LOCO variable importance
  
  Theta_list = list()
  glasso_results = list()
  LOCO_results = list()
  for (i in 1:Iter){
    
    #aa = chol(Sigma)
    Theta = generate_precision(p=p, type = type)
    Theta_list[[i]] = Theta
    Sigma=solve(Theta)
    
    Mu=rep(0,p)
    X=rmvn(n=n, mu = Mu,sigma = Sigma)
    
    S <- var(X)
    
    a = CVglasso(X=X, S=S)
    glasso_results[[i]] = a$Omega
    
    TS_sigma = graph_TS(S = S, n_rho = 50)
    LOCO_results[[i]] = TS_sigma
    #TS_sd = TS_sigma/sum(TS_sigma, na.rm=TRUE)
    #diag(TS_sd) = 0
    
    #results[i, 1] = mean( (Theta != 0) == (a$Omega != 0) )
    #results[i, 2] = mean( (Theta != 0) == (ifelse(TS_sd>thresh, 1,0) != 0) )
    cat('Now:',i, '\n')
  }
  return(list(Theta_list = Theta_list, glasso_results = glasso_results, LOCO_results = LOCO_results))
}


###### generate 2 types of precision matrix as described in GRASS paper #######
generate_precision = function(p=50, type = 'A'){
  # Generate 2 types of precision matrix as described in GRASS paper
  # Args: 
  # p: int, dimension of matrix
  # type: allow 2 options ('A', 'C'), 2 types of matices
  # Returns: matrix, precision matrix
  # 
  
  if(type == 'A'){
    A = matrix(ifelse(runif(p*p, min=0,max=1)<=0.01, 1, 0), p, p)
    non_0_edge = which(A == 1, arr.ind=TRUE)
    upper_edge =  non_0_edge[which((non_0_edge[,2] - non_0_edge[,1]) > 0), ] 
    if(is.matrix(upper_edge)){
      lower_edge = upper_edge
      lower_edge[,1] = upper_edge[,2]; lower_edge[,2] = upper_edge[,1]
      random_vector = runif(dim(upper_edge)[1], min=-0.3,max=0.7)
      A_prime = matrix(0, p, p)
      A_prime[upper_edge] = random_vector
      A_prime[lower_edge] = random_vector
      diag(A_prime) = 1
      
    }else{
      lower_edge[1] = upper_edge[2]; lower_edge[2] = upper_edge[1]
      random_vector = runif(1, min=-0.3,max=0.7)
      A_prime = matrix(0, p, p)
      A_prime[upper_edge[1],upper_edge[2]] = random_vector
      A_prime[upper_edge[2],upper_edge[1]] = random_vector
      diag(A_prime) = 1
      
    }
    lam_min = min(eigen(A_prime)$values)
    Theta = A_prime + (0.1-lam_min)*diag(rep(1,p))
    Theta = Theta / Theta[1,1]
    
  }else if (type == 'C'){
    A_prime = matrix(0, p, p)
    for(i in 1:(p-2)){
      random_value = runif(2, min=-0.3,max=0.7)
      A_prime[i,i+1] = random_value[1]
      A_prime[i,i+2] = random_value[2]
      A_prime[i+1,i] = A_prime[i,i+1]
      A_prime[i+2,i] = A_prime[i,i+2]
    }
    
    lam_min = min(eigen(A_prime)$values)
    Theta = A_prime + (0.1-lam_min)*diag(rep(1,p))
    Theta = Theta / Theta[1,1]
    
  }else{
    stop("Only support type = 'A' or 'C' for now" )
  }
  return(Theta)
  
}



###### graphical model variable screen (my own version, not for the GRASS paper) #######
simu_graph_screen_v2 = function(Theta, Iter = 200){
  # graphical model variable screen, output will be used to generete ROC curve
  # Args:
  # THeta: matrix, the true precision matrix
  # Iter: int, number if iterations, use defualt 50
  # Returns:
  #  A List of: Theta_list: a list of all randomly generated precision matrix
  #             glasso_results: a list of glasso estimates
  #             LOCO_results: a list of calulated LOCO variable importance
  Sigma=solve(Theta)
  
  Theta_list = list()
  glasso_results = list()
  LOCO_results = list()
  
  for (i in 1:Iter){
    Theta_list[[i]] = Theta
    
    Mu=rep(0,p)
    X=rmvn(n=100, mu = Mu,sigma = Sigma)
    
    S<- var(X)
    
    a = CVglasso(X=X, S=S)
    
    glasso_results[[i]] = a$Omega
    
    TS_sigma = graph_TS(S = S, n_rho = 50)
    LOCO_results[[i]] = TS_sigma
    #TS_sd = TS_sigma/sum(TS_sigma, na.rm=TRUE)
    #diag(TS_sd) = 0
    
    #results[i, 1] = mean( (Theta != 0) == (a$Omega != 0) )
    #results[i, 2] = mean( (Theta != 0) == (ifelse(TS_sd>thresh, 1,0) != 0) )
    cat('Now:',i, '\n') 
   
  }
  
  return(list(Theta_list = Theta_list, glasso_results = glasso_results, LOCO_results = LOCO_results))
}



#### A wrapper function to calculate all precision matrix LOCO statistic for all entries of a precision matrix ###
graph_TS = function(S, ...){
  # A wrapper function to calculate all LOCO statistic for all entries of a matrix
  # Args: 
  # S: matrix, emprical covariance matrix
  # ..., arguments for ExactNet.TS.Graph
  # Returns: matrix of LOCO statistic
  p = dim(S)[1]
  TS = matrix(NA, p, p)
  for (i in 1:(p-1)){
    for (j in seq(i+1, p)){
      #for (i in 1:(p)){
      #  for (j in seq(1, p)){
      cat(i, j, '\n')
      try({
        ts = ExactNet.TS.Graph(S, c1 = i, c2 = j, ...)
      })
      TS[i,j] = ts; TS[j,i] = ts; 
    }
  }
  return(TS)
}


#### A wrapper function to calculate all covariance matrix LOCO statistic for all entries of of a matrix ###
graph_TS_sigma = function(S, ...){
  # A wrapper function to calculate all LOCO statistic for all entries of a  covariance matrix
  # Args: 
  # S: matrix, emprical covariance matrix
  # ..., arguments for ExactNet.TS.Graph.Sigma
  # Returns: matrix of LOCO statistic
  
  p = dim(S)[1]
  TS = matrix(NA, p, p)
  for (i in 1:(p-1)){
    for (j in seq(i+1, p)){
      #for (i in 1:(p)){
      #  for (j in seq(1, p)){
      cat(i, j, '\n')
      try({
        ts = ExactNet.TS.Graph.Sigma(S, c1 = i, c2 = j, ...)
      })
      TS[i,j] = ts; TS[j,i] = ts; 
    }
  }
  return(TS)
}


#### A function to transform the outputs of graphical LASSO code
#### to a matrix so we can compute the LOCO statistic easier
trans_lars_matrix = function(obj, use.glasso = FALSE){
  # Args:
  # obj: cound be glasso output or graphLASSO.path output
  # use.glasso: if obj is glasso output, use FALSE. Otherwise, set to be TRUE
  # Return: a matrix,, each row contains the upper triangle value for a fixed lambda
  if (use.glasso){
    
    p = dim(obj$wi)[1]
    n_lam = dim(obj$wi)[3]
    
    beta_hat = matrix(0, n_lam, p*(p-1)/2)
    for(i in 1:n_lam){
      mat_temp = obj$wi[,,i]
      beta_hat[i, ] = mat_temp[upper.tri(mat_temp, diag = FALSE)]
    }
    
  }else{
    
    p = dim(obj$wi[[1]])[1]
    n_lam = length(obj$rholist)
    
    beta_hat = matrix(0, n_lam, p*(p-1)/2)
    for(i in 1:n_lam){
      mat_temp = obj$wi[[i]]
      beta_hat[i, ] = mat_temp[upper.tri(mat_temp, diag = FALSE)]
    }
    
  }
  
  return(beta_hat)
}



######### Our own version of graphical LASSO ########
graphLASSO = function(S, lambda, maxIt = 100, tol = 1e-6, use.lars=FALSE){
  # Our own version of graphical LASSO
  # Args:
  # S: int. sample covariance matrix
  # lambda: positive float, the regularization parameter
  # maxIt: int, default 100, maximum itersion allowed
  # tol: float, tolerance criterion, if results between 2 itersion are < tol, quit from intersion
  # use.lars: use lars or glment as backend. TRUE for lars, FALSE for glmnet
  #           Defulat FALSE (recommened), glmnet runs faster than lars
  # Returns:
  #  A list of: lambda: lambda value, 
  #             W: estimated covariance matrix
  #             Theta: estimated precision matrix, invert of W
  
  p = dim(S)[1]
  W = S + diag(rep(lambda, p))
  W_old = W
  i = 0
  
  while (i < maxIt){
    i = i+1
    for (j in p:1){
      sub_index = (1:p)[-j]
      eigen_obj = eigen(W[sub_index, sub_index], symmetric = TRUE)
      V = eigen_obj$vectors
      d = eigen_obj$values
      X = V %*% diag(sqrt(d)) %*% t(V) #  W_11^(-1/2) * s_12
      Y = V %*% diag( 1/sqrt(d) ) %*% t(V) %*% S[sub_index,j]
      
      if(use.lars){
        lars_obj = lars(X, Y, type='lasso', intercept=FALSE, normalize=FALSE, use.Gram = FALSE) 
        b = predict(lars_obj, s=lambda, type='coefficients', mode='lambda')$coefficients
      }else{
        b = coef( glmnet(X,Y,lambda = lambda/(p-1), intercept = FALSE, standardize = FALSE) )[-1]  
      }
      
      
      W[sub_index, j] = W[sub_index, sub_index] %*% b
      W[j, sub_index] = t(W[sub_index, j])
    }
    #print(W)
    if( norm(W-W_old, type = '1') < tol ){
      break
    }
    
    W_old = W
  }
  if (i == maxIt){
    warning('Maximum number of iteration reached, glasso may not converge.')
  }
  return(list(lambda=lambda, W = W, Theta = solve(W)))  
}


###### Our own version of graphical LASSO, with constraint Sigma_{c1,c2} = 0 #######
graphLASSO.c = function(S, lambda, c1, c2, maxIt = 100, tol = 1e-6, use.lars=FALSE){
  # Our own version of graphical LASSO, with constraint Sigma_{c1,c2} = 0
  # Args:
  # S: int. sample covariance matrix
  # lambda: positive float, the regularization parameter
  # c1, c2: int, specify the correspoing entry of under constraint Sigma_{c1,c2} = 0 
  # maxIt: int, default 100, maximum itersion allowed
  # tol: float, tolerance criterion, if results between 2 itersion are < tol, quit from intersion
  # use.lars: use lars or glment as backend. TRUE for lars, FALSE for glmnet
  #           Defulat FALSE (recommened), glmnet runs faster than lars
  # Returns:
  #  A list of: lambda: lambda value, 
  #             W: estimated covariance matrix
  #             Theta: estimated precision matrix, invert of W
  
  if(c1==c2){
    stop('wrong input of c1,c2, must be different')
  }else{
    c1=min(c1,c2)
    c2=max(c1,c2)
  }
  
  p = dim(S)[1]
  W = S + diag(rep(lambda, p))
  W_old = W
  i = 0
  
  while (i < maxIt){
    i = i+1
    for (j in p:1){
      if(j == c2){
        sub_index = (1:p)[-c(c1,c2)]
        eigen_obj = eigen(W[sub_index, sub_index], symmetric = TRUE)
        V = eigen_obj$vectors
        d = eigen_obj$values
        X = V %*% diag(sqrt(d)) %*% t(V) #  W_11^(-1/2) * s_12
        Y = V %*% diag( 1/sqrt(d) ) %*% t(V) %*% S[sub_index,j]
        if(use.lars){
          lars_obj = lars(X, Y, type='lasso', intercept=FALSE, normalize=FALSE, use.Gram = FALSE) 
          b = predict(lars_obj, s=lambda, type='coefficients', mode='lambda')$coefficients
        }else{
          b = coef( glmnet(X,Y,lambda = lambda/(p-2), intercept = FALSE, standardize = FALSE) )[-1]  
        }
        #w_12_constraint = W[sub_index, sub_index] %*% b
        W[sub_index, j] = W[sub_index, sub_index] %*% b
        W[j, sub_index] = t(W[sub_index, j])
        W[c1, j] = 0
        W[j, c1] = 0
        
      }else if (j == c1){
        sub_index = (1:p)[-c(c1,c2)]
        eigen_obj = eigen(W[sub_index, sub_index], symmetric = TRUE)
        V = eigen_obj$vectors
        d = eigen_obj$values
        X = V %*% diag(sqrt(d)) %*% t(V) #  W_11^(-1/2) * s_12
        Y = V %*% diag( 1/sqrt(d) ) %*% t(V) %*% S[sub_index,j]
        if(use.lars){
          lars_obj = lars(X, Y, type='lasso', intercept=FALSE, normalize=FALSE, use.Gram = FALSE) 
          b = predict(lars_obj, s=lambda, type='coefficients', mode='lambda')$coefficients
        }else{
          b = coef( glmnet(X,Y,lambda = lambda/(p-2), intercept = FALSE, standardize = FALSE) )[-1]  
        }
        W[sub_index, j] = W[sub_index, sub_index] %*% b
        W[j, sub_index] = t(W[sub_index, j])
        W[c2, j] = 0
        W[j, c2] = 0
        
      }else{
        sub_index = (1:p)[-j]
        eigen_obj = eigen(W[sub_index, sub_index], symmetric = TRUE)
        V = eigen_obj$vectors
        d = eigen_obj$values
        X = V %*% diag(sqrt(d)) %*% t(V) #  W_11^(-1/2) * s_12
        Y = V %*% diag( 1/sqrt(d) ) %*% t(V) %*% S[sub_index,j]
        
        if(use.lars){
          lars_obj = lars(X, Y, type='lasso', intercept=FALSE, normalize=FALSE, use.Gram = FALSE) 
          b = predict(lars_obj, s=lambda, type='coefficients', mode='lambda')$coefficients
        }else{
          b = coef( glmnet(X,Y,lambda = lambda/(p-1), intercept = FALSE, standardize = FALSE) )[-1]  
        }
        
        
        W[sub_index, j] = W[sub_index, sub_index] %*% b
        W[j, sub_index] = t(W[sub_index, j])
      }
      
    }
    #print(W)
    if( norm(W-W_old, type = '1') < tol ){
      break
    }
    
    W_old = W
  }
  if (i == maxIt){
    warning('Maximum number of iteration reached, glasso may not converge.')
  }
  return(list(lambda=lambda, W = W, Theta = solve(W)))  
}

####################################################################################

##### Wrapper function the graphLASSO function, genereate graph LASSO solution path ##########
graphLASSO.path = function(S, rholist = NULL, n_rho = 10, ...){
  # Wrapper function the graphLASSO function, genereate graph LASSO solution path
  # Args:
  # S: matrix, sample covariance matrix
  # rholist: vector, sequence of rho, default NULL, if specified, need to be a vector
  # n_rho: int, length of the sequnce of rho, default 10
  # Returns: a list of:
  #             w: list of estimated covariance matrix
  #             wi: list of estimated precision matrix
  #             rholist: a vector of sequence of rho         
  #
  if (is.null(rholist)) {
    rholist = seq(max(abs(S))/n_rho, max(abs(S)), length = n_rho)
  }
  w = list(); wi = list()
  i = 0
  for (rho in rholist){
    i = i + 1
    glasso_obj = graphLASSO(S = S, lambda = rho, ...)
    w[[i]] = glasso_obj$W
    wi[[i]] = glasso_obj$Theta
  }
  return(list(w = w, wi = wi, rholist = rholist))
}


# Wrapper function the graphLASSO.c function, genereate graph LASSO solution path under constraint Sigma_{c1,c2} = 0 #######
graphLASSO.path.c = function(S, rholist = NULL, n_rho = 10, ...){
  # Args:
  # S: matrix, sample covariance matrix
  # rholist: vector, sequence of rho, default NULL, if specified, need to be a vector
  # n_rho: int, length of the sequnce of rho, default 10
  # ..., other arguments for function graphLASSO.c
  # Returns: a list of:
  #             w: list of estimated covariance matrix
  #             wi: list of estimated precision matrix
  #             rholist: a vector of sequence of rho         
  #
  if (is.null(rholist)) {
    rholist = seq(max(abs(S))/n_rho, max(abs(S)), length = n_rho)
  }
  w = list(); wi = list()
  i = 0
  for (rho in rholist){
    i = i + 1
    glasso_obj = graphLASSO.c(S = S, lambda = rho, ...)
    w[[i]] = glasso_obj$W
    wi[[i]] = glasso_obj$Theta
  }
  return(list(w = w, wi = wi, rholist = rholist))
}

# Wrapper function the glasso, genereate graph LASSO solution path #######
glassopath.g = function(S, rholist = NULL, n_rho = 10, ...){
  # Args:
  # S: matrix, sample covariance matrix
  # rholist: vector, sequence of rho, default NULL, if specified, need to be a vector
  # n_rho: int, length of the sequnce of rho, default 10
  # ..., other arguments for function glasso
  # Returns: a list of:
  #             w: list of estimated covariance matrix
  #             wi: list of estimated precision matrix
  #             rholist: a vector of sequence of rho         
  #
  
  
  if (is.null(rholist)) {
    rholist = seq(max(abs(S))/n_rho, max(abs(S)), length = n_rho)
  }
  p = dim(S)[1]
  w = list(); wi = list()
  i = 0
  for (rho in rholist){
    i = i + 1
    
    glasso_obj = glasso(s = S, rho = rho, ...)
    w[[i]] = glasso_obj$w
    wi[[i]] = glasso_obj$wi
  }
  return(list(w = w, wi = wi, rholist = rholist))
}



# Constriant version: wrapper function the glasso function, genereate graph LASSO solution path under constraint Sigma_{c1,c2} = 0 #######
glassopath.c = function(S, c1 = 1, c2 = 2, large_pen = 10e4, rholist = NULL, n_rho = 10, ...){
  # Args:
  # S: matrix, sample covariance matrix
  # c1, c2: int, specify the correspoing entry of under constraint Sigma_{c1,c2} = 0 
  # large_pen: the penalty on the (c1,c2) entry of S, should be large number
  #            just use the default 10^4 unless you are experimenting 
  # rholist: vector, sequence of rho, default NULL, if specified, need to be a vector
  # n_rho: int, length of the sequnce of rho, default 10
  # ..., other arguments for function glasso
  # Returns: a list of:
  #             w: list of estimated covariance matrix
  #             wi: list of estimated precision matrix
  #             rholist: a vector of sequence of rho         
  #
  
  
  if (is.null(rholist)) {
    rholist = seq(max(abs(S))/n_rho, max(abs(S)), length = n_rho)
  }
  p = dim(S)[1]
  w = list(); wi = list()
  i = 0
  for (rho in rholist){
    i = i + 1
    rho_matrix = matrix(rho, p, p)
    rho_matrix[c1,c2] = large_pen
    rho_matrix[c2,c1] = large_pen
    
    glasso_obj = glasso(s = S, rho = rho_matrix, ...)
    w[[i]] = glasso_obj$w
    wi[[i]] = glasso_obj$wi
  }
  return(list(w = w, wi = wi, rholist = rholist))
}



