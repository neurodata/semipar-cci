require("TensorEmbedding")

n<- 20
p<- 5

k<- 5 # reduced dimension

P<- matrix(rnorm(n*k),n,k)
A_mean<- P%*%t(P)
# A_mean<- (A_mean+t(A_mean))/2

A <- array(0, dim = c(n,n,p))
for(i in 1:p){
  noise<- matrix(rnorm(n*n)*0.1,n,n)
  noise<- (noise+t(noise))/2
  A[,,i]<-  A_mean + noise
}


testObj<- TensorEmbedding::symmetric_tensor_decomp(A,n,p, k, 1000, 1E-4, 1E-4)

# A-testObj$A

L <-  testObj$L
C <- testObj$C

image(A_mean, zlim=c(-20,20))
image(L%*% C[,,3]%*%t(L), zlim=c(-20,20))


# grad_L<-function(L,C){
#     grad<- L * 0
#     for(i in 1:p){
#       LC <- L%*%C[,,i]
#       diff<- LC%*%t(L)-A[,,i]
#       grad<- grad + 4*diff %*% LC
#     }
#    grad
# }


# grad_C<-function(L,C){
#   grad<- array(0,dim =c(k,k,p))
#   n<- dim(C)[1]
  
#   for(m in 1:p){
#     LC <- L%*%C[,,m]
#     diff<- LC%*%t(L)-A[,,m]
#     for(i in 1:k){
#       for(j in 1:k){
#         K<- matrix(0, k,k)
#         K[i,j]<- 1
#         K[j,i]<- 1
#         grad[i,j,m]<- sum(2*diff* (L%*%K%*%t(L)))
#       }
#     }
#   }
#   grad
# }


# max(abs(grad_C(L,C) - testObj$gradC))