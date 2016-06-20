require("TensorEmbedding")

n<- 300
p<- 30

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


#avg 1

sum1<- matrix(0, n,n)
sum2<- matrix(0, n,n)
for(i in 1:(p/2)){
  sum1<- sum1+simA[,,i]
  sum2<- sum2+simA[,,p/2+i]
  }

mean_g1<- sum1/(p/2)
mean_g2<- sum2/(p/2)

svd(mean_g1)
