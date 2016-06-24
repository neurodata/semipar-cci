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
  A[,,i]<-  A_mean1 #+ noise
  else
    A[,,i]<-  A_mean2 # + noise
}

propMat<- 1/(1+exp(-A))

simA <- (array(runif(n*n*p), dim = c(n,n,p)))

for(i in 1:p){
  simA[,,i]<-  (simA[,,i] + t(simA[,,i]))/2
}

simA <- (simA < propMat)*1

k=5
testObj<- TensorEmbedding::symmetric_tensor_decomp(simA,n,p, k, 1000, 1E-5, 1E-5, logistic = T, tol = 1E-12,restrictCoreToDiag = TRUE)

L <-  testObj$L
C <- testObj$C


logit<- function(x){
  1/(1+exp(-x))
}

plot(propMat[,,1], logit(L%*%C[,,1]%*%t(L)),xlim=c(0,1),ylim=c(0,1))


# 
k=5
testObj2<- TensorEmbedding::symmetric_decomp_shared_core(simA,n,p, k, 100, 1E-4, 1E-4, logistic = T, tol = 1E-3,restrictCoreToDiag = TRUE)

L2 <-  testObj2$L
C2 <- testObj2$C
 

plot(propMat[,,1], logit(L2[,,1]%*%C2%*%t(L2[,,1])),xlim=c(0,1),ylim=c(0,1))


library(ROCR)


pred <- prediction(c(logit(L%*%C[,,1]%*%t(L))), c(simA[,,1]))
perf <- performance(pred, "tpr", "fpr")
performance(pred,"auc")
plot(perf)


pred <- prediction(c(logit(L2[,,1]%*%C2%*%t(L2[,,1]))), c(simA[,,1]))
perf <- performance(pred, "tpr", "fpr")
performance(pred,"auc")
plot(perf)



pred <- prediction(c(logit(A[,,1])), c(simA[,,1]))
perf <- performance(pred, "tpr", "fpr")
performance(pred,"auc")
plot(perf)

require(ggplot2)
require(reshape)

cMelted<- numeric()

for(i in 1:p){
  cGroup1 <- melt((C[,,i]))
  cMelted<- rbind( cMelted, cbind(cGroup1, "group"=i))
}


p1 <- ggplot(cMelted, aes(X1, X2)) + geom_tile(aes(fill = value),
                                               colour = "white") + scale_fill_gradient(low = "white",high = "red") + facet_wrap(~group, ncol=5)


p1

out<- seq(range(C)[1], range(C)[2], length.out = k)
# 
# 
plot(c(1:k), out, type = "n")
for(i in 1:(p/2)){
  lines(c(1:k), diag(C[,,i]), col = "blue")
}
for(i in (p/2+1):p){
  lines(c(1:k), diag(C[,,i]), col = "red")
}