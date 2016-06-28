require("TensorEmbedding")


#dimensions

n<- 100

#latent reduced dimension 
k<- 3


#number of groups
mG <- 10
weight<- rep(4, mG)

g<- numeric()

for(i in 1:mG){
  g <- c( g, rep(i, weight[i]))
}

#total number of subjects
m <- length(g)

#data tensor
A <- array(0, dim = c(n,n,m))




# A_mean<- (A_mean+t(A_mean))/2


g0<- g[1]-1

for(i in 1:m){

  if(g[i]!=g0){
    g0<- g[i]
    P<- matrix(rnorm(n*k,sd = 0.5),n,k)
    A_mean1<- P%*%t(P)
  }
    
    A[,,i]<-  A_mean1
}

propMat<- 1/(1+exp(-A))

simA <- array(runif(n*n*m), dim = c(n,n,m))

for(i in 1:m){
  simA[,,i]<-  (simA[,,i] + t(simA[,,i]))/2
}

simA <- (simA < propMat)*1


testObj<- TensorEmbedding::symm_group_tensor_decomp(simA, g, n, m, k, 10, 1E-4, 1E-4, tol = 1E-12,restrictCoreToDiag = FALSE)

L <-  testObj$L
C <- testObj$C


m
g