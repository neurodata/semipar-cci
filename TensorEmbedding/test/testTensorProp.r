require("TensorEmbedding")

n<- 30
p<- 10

k<- 5 # reduced dimension

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


testObj<- TensorEmbedding::symmetric_tensor_decompEM(simA,n,p, k, steps = 5000, 1E-4, 1E-4, tol = 1E-12,restrictCoreToDiag = TRUE,loss_type = 1)

# testObj<- TensorEmbedding::symmetric_tensor_decomp(simA,n,p, k, steps = 1000, 1E-4, 1E-3, tol = 1E-12,restrictCoreToDiag = TRUE)

L <-  testObj$L
C <- testObj$C

require(ggplot2)
require(reshape)

# C[C>100]<-0

cMelted<- numeric()

for(i in 1:p){
  cGroup1 <- melt(C[, , i])
  cMelted<- rbind( cMelted, cbind(cGroup1, "group"=i))
}


p1 <- ggplot(cMelted, aes(X1, X2)) + geom_tile(aes(fill = value),
      colour = "white") + scale_fill_gradient(low = "white",high = "red") + facet_wrap(~group, ncol=5)

p1
# pdf("compare_2_groups.pdf",8,6)
# p1
# dev.off()

out<- seq(range(C)[1], range(C)[2], length.out = k)


plot(c(1,k), range(C), type = "n")
for(i in 1:(p/2)){
  lines(c(1:k), C[i,], col = "blue")
}
for(i in (p/2+1):p){
  lines(c(1:k), C[i,], col = "red")
}

testObj<- TensorEmbedding::symmetric_tensor_decomp(simA,n,p, k, steps = 1000, 1E-4, 1E-3, tol = 1E-12,restrictCoreToDiag = TRUE)


image(simA[,,1])
image(1/(1+exp(-L%*% diag(C[1,])%*%t(L))),zlim=c(0,1))


plot( (L%*% diag(C[1,]) %*%t(L)), A[,,1], xlim=c(-5,5),ylim = c(-5,5))
plot( testObj$L%*% testObj$C[,,1]%*%t(L), A[,,1], xlim=c(-5,5),ylim = c(-5,5))


plot( 1/(1+exp(-L%*% diag(C[1,]) %*%t(L))),  1/(1+exp(-A[,,1])),xlim =  range((A[,,1])))
