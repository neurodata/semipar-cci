setwd("c://work/git/semipar-cci/BatchEffectRemoval/")

require("rstiefel")
require("msm")

loadData<- function(path){
  load(path)
  Data1
}

DataList = list()
DataList[[1]] = loadData("../Data/processed/BNU1.RDa")
DataList[[2]] = loadData("../Data/processed/BNU3.RDa")
DataList[[3]] = loadData("../Data/processed/KKI2009.RDa")
DataList[[4]] = loadData("../Data/processed/MRN114.RDa")
DataList[[5]] = loadData("../Data/processed/SWU4.RDa")

source("stiefel_diagonalize.r")

y_list = list()
g=numeric()
age=numeric()
for(i in 1:5){
  y_list= c(y_list,DataList[[i]]$A)
  g = c(g, DataList[[i]]$SEX)
  age = c(age, DataList[[i]]$age)
}

pick = age==20

y_list_pick = list()
g_pick = g[pick]
age_pick = age[pick]

for(i in 1:length(y_list)){
  if(pick[i]){
    y_list_pick[[length(y_list_pick)+1]] = y_list[[i]]
  }
}


r=30
stiefelDecomp = stiefel_diagonalize(y_list_pick, r=r,diagD = T)

estP = getPEst(stiefelDecomp$U0,  stiefelDecomp$D_list)

image(estP[[1]])


D_list = stiefelDecomp$D_list
ts.plot(stiefelDecomp$trace_D[,1:5])

D = colMeans(stiefelDecomp$trace_D[101:200,])



m = length(y_list_pick)
D<- matrix(D,r*r)
plot(c(1:r*r), seq(min(D),max(D),length.out = r),type = "n")
for(j in 1:m){
  lines(c(1:(r*r)),D[,j],col=g_pick[j])
}

# plot(colMeans(trace_D))

source("classification.r")

#knn error based on adjacency matrix
A_lowtri = sapply(y_list_pick,function(x){x[lower.tri(x)]})
knn_raw_error = knnLOOVaryingK(M = t(A_lowtri), label = g_pick,max_k = 20)

#knn error based on core
knn_error = knnLOOVaryingK(M = t(D), label = g_pick,max_k = 20)

save(knn_raw_error,file= "knn_raw_error.RDa")
save(knn_error,file= "knn_error.RDa")

load("knn_raw_error.RDa")
load("knn_error.RDa")

pdf("knnErrorEmbeddingVsRaw.pdf",8,6)
plot(knn_error,ylim=c(0,1),type="l",col="red")
lines(knn_raw_error,col="blue")
dev.off()

min(knn_error)
min(knn_raw_error)

###test
# 
# U0 = stiefelDecomp$U0
# 
# D = stiefelDecomp$D_list[[1]]
#   
# fit1<- U0%*%diag(diag(D))%*%t(U0)
# fit2<- U0%*%D%*%t(U0)

# fit1-fit2


require("BayesLogit")
require("pROC")

X = cbind(1,t(D))
logit_fit<- logit(g_pick-1,X,P0 = diag(1E-3,(r^2+1)))
fit_p<- 1/(1+exp(-X%*%colMeans(logit_fit$beta)))
roc_core<- roc(g_pick-1,c(fit_p))


X = cbind(1,t(A_lowtri))
logit_fit<- logit(g_pick-1,X,P0 = diag(1E-3,ncol(X)))
fit_p2<- 1/(1+exp(-X%*%colMeans(logit_fit$beta)))
roc_core2<- roc(g_pick-1,c(fit_p2))



