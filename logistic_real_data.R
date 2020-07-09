




## Leukemia dataset ##

require(LOCOpath)
require(hdi)

leukemia = readRDS('Golub1999.rds')
str(leukemia)

## a utiltily function for parrallel computing 
Logistic_TS_util_fun = function(x_sp, y_sp, which.covariate = 1, betaNull = 0, multiTest = FALSE, norm = 'L1', ...){
  return(
    ExactNet.TS.Logistic(X = x_sp, Y = y_sp, which.covariate = which.covariate, betaNull = betaNull,
                 multiTest = multiTest, ...)
  )
}

## normalization 
leukemia$X = scale(leukemia$X)
leukemia$y = ifelse(leukemia$y == 'ALL', 1, 0)

## start parrallel computing 
n_threads = detectCores()
cl = makeCluster(n_threads, type = "FORK")

TS=unlist(parLapply(cl, X=1:dim(leukemia$X)[2], Logistic_TS_util_fun, 
                    x_sp = leukemia$X, y_sp = leukemia$y, 
                    betaNull = 0, multiTest = FALSE, norm = 'L1'))

index = which(TS > 0)


## variable screening 
x = leukemia$X[, index]
y = ifelse(leukemia$y == 'ALL', 1, 0)

### calculate p-value after variable screening

obj = list()
for (i in 1:length(index)){
  print(i)
  try({
    obj[[i]] = Net.Resample.Logistic(X = leukemia$X[, index], Y = leukemia$y, which.covariate = i,
                           betaNull = 0, multiTest = FALSE, B = 500, 
                           beta.init = 'adaptive', beta.null.estimate = TRUE,
                           beta.true = 0)
    print(obj[[i]]$pval)
  })
}
save(obj, file = 'pval_leuk.RData')

pval_leukemia = lapply(obj, FUN = function(x){x$pval}) # extract p-value

## desparsified LASSO
lasso.proj(leukemia$X, ifelse(leukemia$y == 'ALL', 1, 0), standardize = TRUE, parallel = TRUE, family = 'binomial', ncores = 28)


## prostate dataset ##
prostate = readRDS('Singh2002.rds')

## normalization
prostate$X = scale(prostate$X)
prostate$y = ifelse(prostate$y == 'normal', 1, 0)

## start parrallel computing 
n_threads = detectCores()
cl = makeCluster(n_threads, type = "FORK")

TS=unlist(parLapply(cl, X=1:dim(prostate$X)[2], Logistic_TS_util_fun, 
                    x_sp = prostate$X, y_sp = prostate$y,
                    betaNull = 0, multiTest = FALSE, norm = 'L1'))

index = which(TS > 0)

save(TS, file = 'TS_prostate.RData')


### calculate p-value for genes screened in
index = which(TS > 0)

pp = list()
for (i in c(19,37,49)){
  print(i)
  try({
    pp[[i]] = Net.Resample.Logistic(X = prostate$X[, index], Y = prostate$y, which.covariate = i,
                                     betaNull = 0, multiTest = FALSE, B = 50000, 
                                     beta.init = 'adaptive', beta.null.estimate = FALSE,
                                     beta.true = 0)$pval
    print(pp[[i]]$pval)
  })
}


save(obj, file = 'pval_prostate.RData')

pval_prostate = unlist(lapply(obj, FUN = function(x){x$pval})) # extract p-value

colnames(prostate$X)[ index[which(pval_prostate<0.05)] ] 

# yield 2 significant genes
# "X83543" "X07732"
# 0.0014 < 0.0001 


## desparsified LASSO
proj_prostate = lasso.proj(prostate$X, ifelse(prostate$y == 'normal', 1, 0), standardize = TRUE, parallel = TRUE, family = 'binomial', ncores = 28)

##
#0.003053213
#X07732
#6185
#



###### colon dataset ######

# load the package to get the dataset
require(HiDimDA)
## load data
data(AlonDS)

# preprocessing 
Y = ifelse(AlonDS$grouping == 'healthy', 0, 1)
X = matrix(0, 62, 2000)
for(i in 2:2001){
 X[,i-1] = AlonDS[[i]] 
}  
X = log10(X)
X = scale(X)

Logistic_TS_util_fun = function(x_sp, y_sp, which.covariate = 1, betaNull = 0, multiTest = FALSE, norm = 'L1', ...){
  return(
    ExactNet.TS.Logistic(X = x_sp, Y = y_sp, which.covariate = which.covariate, betaNull = betaNull,
                         multiTest = multiTest, ...)
  )
}

## parallel computing starts
n_threads = detectCores()
cl = makeCluster(n_threads, type = "FORK")

TS=unlist(parLapply(cl, X=1:dim(X)[2], Logistic_TS_util_fun, 
                    x_sp = X, y_sp = Y, 
                    betaNull = 0, multiTest = FALSE, norm = 'L1'))

save(TS, file = 'TS_colon.RData')

index = which(TS > 0)

### calculate p-value for genes screened in
obj = list()
for (i in 1:length(index)){
  print(i)
  try({
    obj[[i]] = Net.Resample.Logistic(X = X[,index], Y = Y, which.covariate = i,
                                     betaNull = 0, multiTest = FALSE, B = 500, 
                                     beta.init = 'adaptive', beta.null.estimate = FALSE,
                                     beta.true = 0)
    print(obj[[i]]$pval)
  })
}
save(obj, file = 'pval_colon.RData')

pval_view = unlist(lapply(obj, FUN = function(x){x$pval})) # extract p-value

## depsarsified LASSO 
lasso.proj(X, Y, standardize = TRUE, parallel = TRUE, family = 'binomial', ncores = 28)




#################################
## Lymphoma dataset 
## need to install R packgage KODAMA manually
## The tar ball file can be downloaded from 
# https://cran.r-project.org/src/contrib/Archive/KODAMA/

install.packages('~/R/KODAMA_1.5.tar.gz', repos = NULL, type="source")

## load dataset
require(KODAMA)
data(lymphoma)

## preprocessing
lymphoma$X = scale(lymphoma$data)
lymphoma$y = ifelse(lymphoma$class == 'DLBCL', 1, 0)


## parallel computing
n_threads = detectCores()
cl = makeCluster(n_threads, type = "FORK")

TS=unlist(parLapply(cl, X=1:dim(lymphoma$X)[2], Logistic_TS_util_fun, 
                    x_sp = lymphoma$X, y_sp = lymphoma$y,
                    betaNull = 0, multiTest = FALSE, norm = 'L1'))

save(TS, file = 'TS_lym.RData')


index = which(TS > 0)

### calculate p-value for genes screened in
obj = list()
for (i in 1:length(index)){
  print(i)
  try({
    obj[[i]] = Net.Resample.Logistic(X = lymphoma$X[, index], Y = lymphoma$y, which.covariate = i,
                                     betaNull = 0, multiTest = FALSE, B = 500, 
                                     beta.init = 'adaptive', beta.null.estimate = FALSE,
                                     beta.true = 0)
    print(obj[[i]]$pval)
  })
}
save(obj, file = 'pval_lym.RData')

pval_view = lapply(obj, FUN = function(x){x$pval}) #extract p-values

## desparsified LASSO
lasso.proj(lymphoma$X, lymphoma$y, standardize = TRUE, parallel = TRUE, family = 'binomial', ncores = 28)


############## variable importance bar plot #######################


par(mfrow=c(1,1))

load('./results/TS_prostate.RData')
## TS_temp is for adding the label on the bar plot
TS_temp = TS
names(TS_temp) = colnames(prostate$X)

pdf('TS_prostate.pdf')
barplot(sort(TS/sum(TS), decreasing = TRUE)[1:30], col="grey50", 
        main="",
        ylab="variable importance",
        xlab = " ",
        space=1)
text(seq(1.5, 60, by=2), par("usr")[3], 
     srt = 60, adj= 1, xpd = TRUE,
     labels = (names(sort(TS_temp/sum(TS_temp), decreasing = TRUE)[1:30])), cex=0.5)
dev.off()


load('./results/TS_leukemia.RData')

pdf('TS_leukemia.pdf')
barplot(sort(TS/sum(TS), decreasing = TRUE)[1:50], col="grey50", 
        main="",
        ylab="variable importance",
        xlab = "",
        space=1)
dev.off()

load('./results/TS_colon.RData')

pdf('TS_colon.pdf')
barplot(sort(TS/sum(TS), decreasing = TRUE)[1:50], col="grey50", 
        main="",
        ylab="variable importance",
        xlab = "",
        space=1)
dev.off()


load('./results/TS_lym.RData')

pdf('TS_lym.pdf')
barplot(sort(TS/sum(TS), decreasing = TRUE)[1:50], col="grey50", 
        main="",
        ylab="variable importance",
        xlab = "",
        space=1)

dev.off()



