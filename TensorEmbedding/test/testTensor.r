require("TensorEmbedding")

n<- 300
p<- 10

k<- 5 # reduced dimension

P<- matrix(rnorm(n*k,sd = 0.5),n,k)
A_mean1<- P%*%t(P)

P<- matrix(rnorm(n*k,sd = 0.5),n,k)
A_mean2<- P%*%t(P)

# A_mean<- (A_mean+t(A_mean))/2

A <- array(0, dim = c(n,n,p))
for(i in 1:p){
  noise<- matrix(rnorm(n*n),n,n)*5
  noise<- (noise+t(noise))/2
  if(i<= (p/2))
  A[,,i]<-  A_mean1 + noise
  else
    A[,,i]<-  A_mean2 + noise
  
}

testObj<- TensorEmbedding::symmetric_tensor_decomp(A,n,p, k, 1000, 1E-4, 1E-4)

L <-  testObj$L
C <- testObj$C

# image(A[,,1], zlim=c(-3,3))
# image(A[,,2], zlim=c(-3,3))
# 
# image(A_mean1, zlim=c(-3,3))
# image(L%*% C[,,1]%*%t(L), zlim=c(-3,3))
# 
# 
# hist(A_mean1 - L%*% C[,,1]%*%t(L))
# hist(A_mean2 - L%*% C[,,6]%*%t(L))
# 
# image(C[,,1])
# image(C[,,6])

require(ggplot2)
require(reshape)

cGroup2 <- melt(C[, , 2])

cMelted<- numeric()

for(i in 1:p){
  cGroup1 <- melt(C[, , i])
  cMelted<- rbind( cMelted, cbind(cGroup1, "group"=i))
}

cMelted<- rbind(  cbind(cGroup1, "group"=1)
                  , cbind(cGroup2, "group"= 2))


p1 <- ggplot(cMelted, aes(X1, X2)) + geom_tile(aes(fill = value),
      colour = "white") + scale_fill_gradient(low = "white",high = "red") + facet_wrap(~group, ncol=5)

pdf("compare_2_groups.pdf",8,6)
p1
dev.off()
