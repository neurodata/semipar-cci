require("TensorEmbedding")

n<- 10
p<- 3

k<- 3 # reduced dimension

P<- matrix(rnorm(n*k,sd = 0.5),n,k)
A_mean1<- P%*%t(P)

P<- matrix(rnorm(n*k,sd = 0.5),n,k)
A_mean2<- P%*%t(P)

# A_mean<- (A_mean+t(A_mean))/2

A <- array(0, dim = c(n,n,p))
for(i in 1:p){
  noise<- matrix(rnorm(n*n),n,n)
  noise<- (noise+t(noise))/2
  if(i<= (p/2))
  A[,,i]<-  A_mean1 + noise
  else
    A[,,i]<-  A_mean2  + noise
}

propMat<- 1/(1+exp(-A))

simA <- (array(runif(n*n*p), dim = c(n,n,p)))

for(i in 1:p){
  simA[,,i]<-  (simA[,,i] + t(simA[,,i]))/2
}

simA <- (simA < propMat)*1

# testObj<- TensorEmbedding::symmetric_tensor_decomp(simA,n,p, k, 10, 1E-4, 1E-4, logistic = T, tol = 1E-12,restrictCoreToDiag = FALSE)

 testObj<- TensorEmbedding::symmetric_decomp_shared_core(simA,n,p, k, 10, 1E-4, 1E-4, logistic = T, tol = 1E-12,restrictCoreToDiag = FALSE)

L <-  testObj$L
C <- testObj$C


loglogit<- function(x, theta){
   sum( -x*theta + log(1 + exp(theta)))
}

loss<- function(L,C){
    l<-0
    for(i in 1:p){
        # l<- l+ loglogit(A[,,i],   L %*%C[,,i]%*% t(L))
                l<- l+ loglogit(simA[,,i],   L[,,i] %*%C%*% t(L[,,i]))

    }
    return (l + sum(L*L/1E4)/2 +sum(C*C/1)/2)
}


numGradL<- function(i,j,k){

    L1<- L
    L1[i,j,k]<- L[i,j,k]+1E-8

    (loss(L1,C)- loss(L,C))/1E-8

}



numGradC<- function(i,j){

    C1<- C
    C1[i,j]<- C[i,j]+1E-8
    C1[j,i]<- C[i,j]+1E-8

    (loss(L,C1)- loss(L,C))/1E-8

}

testObj$gradC
numGradC(2,3)

