require("TensorEmbedding")



n<- 20
p<- 100

k<- 10 # reduced dimension

A <- array(rnorm(n*n*p), dim = c(n,n,p))
for(i in 1:p){
  A0<- A[,,i]
  A[,,i]<-  (A0+t(A0))/2
}


testObj<- TensorEmbedding::rcpparma_bothproducts(A,n,p, k)

# A-testObj$A

L <-  testObj$L
C <- testObj$C


grad_L<-function(L,C){
    grad<- L * 0
    for(i in 1:p){
      LC <- L%*%C[,,i]
      diff<- LC%*%t(L)-A[,,i]
      grad<- grad + 4*diff %*% LC
    }
   grad
}


grad_C<-function(L,C){
  grad<- array(0,dim =c(k,k,p))
  n<- dim(C)[1]
  
  for(m in 1:p){
    LC <- L%*%C[,,m]
    diff<- LC%*%t(L)-A[,,m]
    for(i in 1:k){
      for(j in 1:k){
        K<- matrix(0, k,k)
        K[i,j]<- 1
        K[j,i]<- 1
        grad[i,j,m]<- sum(2*diff* (L%*%K%*%t(L)))
      }
    }
  }
  grad
}


max(abs(grad_C(L,C) - testObj$gradC))