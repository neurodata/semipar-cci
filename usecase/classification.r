require("TensorEmbedding")
require("base")
setwd("~/git/semipar-cci/")

load("Data/processed/BNU3.RDa")

m<- length(BNU3)
n<- nrow(BNU3[[1]]$A)

A<- array(0,dim = c(n,n,m))
for(i in 1:m){
  A[,,i]<- BNU3[[i]]$A
}


id<- unlist(lapply(BNU3, function(x){x$id}))

k<- 30
fit <- TensorEmbedding::symm_group_tensor_decomp(A, id, n, m, k, 500, 1E-5, 1E-5)


SEX<- unlist(lapply(BNU3, function(x){x$SEX}))


y<- SEX-1


plot(c(1,k),range(fit$diagC),type="n")
for(i in 1:m){
  lines(fit$diagC[i,],col=y[i]+1)
}

require("ggplot2")

# 
# overlayHist<- function(x,y){
#   dat<- data.frame(
#     "class"=y,
#     "feature"=x
#   )
#   ggplot(dat, aes(x=feature, fill=as.factor(class))) + geom_histogram(alpha=0.2, position="identity", bins = 20)
#   
# }
# 
# overlayHist(fit$diagC[,30],y)




looCV <- function(i){
  
  
  computeProb<-function(x){
    d1<-(x-m1)
    log_prob1<- -(d1%*%solve(cov1,d1) + log(det(cov1)))/2
    
    d0<-(x-m0)
    log_prob0<- -(d0%*%solve(cov0,d0) + log(det(cov0)))/2
    
    prob<-  1/(1+exp(log_prob0-log_prob1))
    prob 
  }
  
  x<- fit$diagC[-i,][y[-i]==1,]
  m1<- colMeans(x)
  cov1<- cov(x)
  diag(cov1)<- diag(cov1)+ 1E-8
  
  x<- fit$diagC[-i,][y[-i]==0,]
  m0<- colMeans(x)
  cov0<- cov(x)
  diag(cov0)<- diag(cov0)+ 1E-8
  
  
  pred<- computeProb(fit$diagC[i,])
  
  1*(as.numeric(pred)>0.5 )== y[i]
}


estCorrect<- unlist(lapply(1:m, looCV))

1-sum(estCorrect)/m
