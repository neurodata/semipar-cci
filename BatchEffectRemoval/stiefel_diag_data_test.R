setwd("~/git/semipar-cci/BatchEffectRemoval/")

require("rstiefel")
require("msm")

load("./RDa/5dataset_list.RDa")
load("./RDa/5dataset_sex.RDa")

source("stiefel_diagonalize.r")
source("classification.r")


y_list = c(A_list[[1]], A_list[[3]],A_list[[5]])
g = c(label_list[[1]], label_list[[3]],label_list[[5]])

r=30
stiefelDecomp = stiefel_diagonalize(y_list, r=30)

D_list = stiefelDecomp$D_list

ts.plot(stiefelDecomp$trace_D[,1:5])

m = length(y_list)
D<- matrix(unlist(D_list),r)
plot(c(1:r), seq(min(D),max(D),length.out = r),type = "n")
for(j in 1:m){
  lines(c(1:r),D[,j],col=g[j])
}

# plot(colMeans(trace_D))

knn_error = knnLOOVaryingK(M = t(D[2:10,]), label = g)


plot(rowMeans(D[,g==1]))
lines(rowMeans(D[,g==2]))

plot(knn_error)

A_lowtri = sapply(y_list,function(x){x[lower.tri(x)]})
knn_raw_error = knnLOOVaryingK(M = t(A_lowtri), label = g)

plot(knn_error,ylim=c(0,1))
lines(knn_raw_error)
