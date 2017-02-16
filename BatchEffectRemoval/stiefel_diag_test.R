setwd("~/git/semipar-cci/BatchEffectRemoval/")

require("rstiefel")
require("msm")

source("stiefel_diagonalize.r")

n<- 70
r<- 10
m<- 30

zeroB<- matrix(0,r,r)
zeroA<- matrix(0,n,n)

U<- matrix(rnorm(n*r),n)
svdU = svd(U)
U= svdU$u

symmetrize<- function(X){
  X= t(X)
  X[lower.tri(X)]<- t(X)[lower.tri(X)]
  X
}

Y_list = list()
Z_list = list()
for(j in 1:m){
  D<- c(1:r)*6+ rnorm(r)
  if(j> (m/2))
    D<- 10/c(1:r)+ rnorm(r)
  UDU<- U%*% diag(D)%*% t(U)
  Z_list[[j]]<- UDU + symmetrize(matrix(rnorm(n*n),n))
  Y_list[[j]]<- (Z_list[[j]]>0)*1
  diag(Y_list[[j]])=2
}
rm(U)
rm(Z_list)

stiefelDecomp = stiefel_diagonalize(Y_list, r)

D_list = stiefelDecomp$D_list

ts.plot(stiefelDecomp$trace_D[,1:5])

D<- matrix(unlist(D_list),r)
g<- (c(1:m)> (m/2))*1+1
plot(c(1:r), seq(min(D),max(D),length.out = r),type = "n")
for(j in 1:m){
  lines(c(1:r),D[,j],col=g[j])
}

# plot(colMeans(trace_D))
