require("TensorEmbedding")
require("gtools")

n<- 100
p<- 30

k<- 5 # reduced dimension
B<- 2 # blocks
weight<- rdirichlet(1, rep(1,B))

# P<- matrix(rnorm(n*k,sd = 0.5),n,k)

pInBlocks<-  matrix(rnorm(B*k,sd = 1),B,k)

Z<- apply(t(rmultinom(n,1,weight)), 1, function(x){
  c(1:B)[as.logical(x)]
})


P <- pInBlocks[Z,]

vecGroup1<- runif(k,min = -1, max = 1)
vecGroup2<- runif(k,min = -1, max = 1)

# plot(vecGroup1, type="l", col="red")
# lines(vecGroup2)

A_mean1<- P%*%diag(vecGroup1)%*%t(P)
A_mean2<- P%*%diag(vecGroup2)%*%t(P)

##
hist(A_mean1)
hist(A_mean2)

# A_mean<- (A_mean+t(A_mean))/2

A <- array(0, dim = c(n,n,p))
for(i in 1:p){
  noise<- matrix(rnorm(n*n),n,n)
  noise<- (noise+t(noise))/2
  if(i<= (p/2))
  A[,,i]<-  A_mean1 #+ noise
  else
    A[,,i]<-  A_mean2#  + noise
}

propMat<- 1/(1+exp(-A))

simA <- (array(runif(n*n*p), dim = c(n,n,p)))

for(i in 1:p){
  l_pos<- lower.tri( simA[,,i],diag = F)
  simA[,,i][l_pos]<-  t(simA[,,i])[l_pos]
}

simA <- (simA < propMat)*1

k_latent = k

# testObj<- TensorEmbedding::symmetric_tensor_decomp(A,n,p, k_latent, 1000, 1E-4, 1E-3, loss_type = 2, tol = 1E-12,restrictCoreToDiag = TRUE)

testObj<- TensorEmbedding::symmetric_tensor_decomp(simA,n,p, k_latent, 1000, 1E-3, 1E-3, logistic = TRUE, tol = 1E-12,restrictCoreToDiag = TRUE)
# testObj<- TensorEmbedding::symmetric_tensor_decomp(simA,n,p, k_latent, 1000, 1E-4, 1E-4, loss_type = 2, tol = 1E-12,restrictCoreToDiag = TRUE)


L <- testObj$L
C <- testObj$C

#plot 


logit<-function(x){
  1/(1+exp(-x))
}

plot( A[,,1], L%*%C[,,1]%*%t(L) ,xlim=c(-3,3),ylim=c(-3,3))
plot( A[,,6], L%*%C[,,6]%*%t(L) ,xlim=c(-3,3),ylim=c(-3,3))



image(simA[,,(1)], zlim=c(0,1))
image(logit(A[,,(1)]),  zlim=c(0,1))
image(logit(L%*%C[,,(1)]%*%t(L)),  zlim=c(0,1))

km_fit_CL<- kmeans(logit(c(L%*%C[,,(1)]%*%t(L))),3)
km_fit_CL$centers
unique(c(logit(A[,,(1)])))






image(simA[,,(p/2+1)], zlim=c(0,1))
image(logit(A[,,(p/2+1)]),  zlim=c(0,1))
image(logit(L%*%C[,,(p/2+1)]%*%t(L)),  zlim=c(0,1))

km_fit_CL<- kmeans(logit(c(L%*%C[,,(p/2+1)]%*%t(L))),3)
km_fit_CL$centers
unique(c(logit(A[,,(p/2+1)])))


km_L<- kmeans(L, centers = B)

P_fit<- km_L$centers[km_L$cluster,]

image(logit(P_fit%*%C[,,(p/2+1)]%*%t(P_fit)))
image(logit(A[,,(p/2+1)]))


table(logit(P_fit%*%C[,,(p/2+1)]%*%t(P_fit)))
table(logit(A[,,(p/2+1)]))

plot(logit(L%*%C[,,(p/2+1)]%*%t(L)), logit(P_fit%*%C[,,(p/2+1)]%*%t(P_fit)))


######
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

# pdf("compare_2_groups.pdf",8,6)
# p1
# dev.off()

out<- seq(range(C)[1], range(C)[2], length.out = k_latent)


plot(c(1:k_latent), out, type = "n" )#, ylim = c(0,50))
for(i in 1:(p/2)){
  lines(c(1:k_latent), diag(C[,,i]), col = "blue")
}
for(i in (p/2+1):p){
  lines(c(1:k_latent), diag(C[,,i]), col = "red")
}



# km_L <- kmeans(L,centers = B)

# block_c<- km_L$centers

# P_fit<- block_c[km_L$cluster,]

# image(P_fit%*%C[,,6]%*%t(P_fit))
