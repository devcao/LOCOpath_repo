

require(CVglasso)
require(huge)

# This package contains the riboflavin dataset
require(hdi)
data("riboflavin")

## from riboflavin dataset to riboflavinV100
var_ribo = diag(var(riboflavin$x))
order_ribo = order(var_ribo, decreasing = TRUE)
riboflavinV100 = riboflavin$x[,order_ribo[1:100]]

## preporcessing normalization
X.npn = huge.npn(riboflavinV100)

## glasso results
glasso_ribo = CVglasso(X=X.npn, S=var(X.npn))

## LOCO path results
TS_ribo = graph_TS(S = var(X.npn), n_rho = 50)

a1 = glasso_ribo$Omega
a2 = TS_ribo/sum(TS_ribo,na.rm = TRUE)  # normalize

## generate the network structure for 4 different thereshold
par(mfrow=c(1,4))
huge.plot_v1(Matrix(ifelse(a2>quantile(a2, prob = 0.9, na.rm=TRUE),1,0)), main = 'q=0.9')

huge.plot_v1(Matrix(ifelse(a2>quantile(a2, prob = 0.95, na.rm=TRUE)
                              ,1,0)), main = 'q=0.95')
huge.plot_v1(Matrix(ifelse(a2>quantile(a2, prob = 0.97, na.rm=TRUE)
                              ,1,0)), main = 'q=0.97')
huge.plot_v1(Matrix(ifelse(a2>quantile(a2, prob = 0.99, na.rm=TRUE)
                              ,1,0)), main = 'q=0.99')



# This package contains the flow cytometry dataset
require(sparsebn)

# load the flow cytometry dataset
data(cytometryContinuous)
# preprocessing
X.npn = huge.npn(cytometryContinuous$data)

# glasso results
glasso_cyto = CVglasso(X=X.npn, S=var(X.npn))
# LOCO path results
TS_cyto = graph_TS(S = var(X.npn), n_rho = 50)

# normalize

a2 = TS_cyto/sum(TS_cyto,na.rm = TRUE) 

## generate the network structure for 4 different thereshold
par(mfrow=c(1,4))
huge.plot_v1(Matrix(ifelse(a2>quantile(a2, prob = 0.75, na.rm=TRUE),1,0)), main = 'q=0.75')

huge.plot_v1(Matrix(ifelse(a2>quantile(a2, prob = 0.80, na.rm=TRUE)
                           ,1,0)), main = 'q=0.80')
huge.plot_v1(Matrix(ifelse(a2>quantile(a2, prob = 0.85, na.rm=TRUE)
                           ,1,0)), main = 'q=0.85')
huge.plot_v1(Matrix(ifelse(a2>quantile(a2, prob = 0.90, na.rm=TRUE)
                           ,1,0)), main = 'q=0.90')



par(mfrow=c(1,1))

TS = a2[upper.tri(a2)]

### generate the variable importance bar plot
pdf('theta_imp_cyto.pdf')
barplot(sort(TS/sum(TS), decreasing = TRUE), col="grey50", 
        main="",
        ylab="variable importance",
        xlab = "",
        space=1)
dev.off()


pdf('theta_imp_ribo_200.pdf')
a2 = TS_ribo/sum(TS_ribo,na.rm = TRUE)
TS = a2[upper.tri(a2)]

barplot(sort(TS/sum(TS), decreasing = TRUE)[1:200], col="grey50", 
        main="",
        ylab="variable importance",
        xlab = "",
        space=1)

dev.off()


plot(sort(TS, decreasing = TRUE)[1:100], type='h', xlab = 'index', ylab = 'variable importance')
points(sort(TS, decreasing = TRUE)[1:100], pch = 19, cex = 0.5)
abline(v=26, col = 'red', lty = 3, lwd = 2)




#######
