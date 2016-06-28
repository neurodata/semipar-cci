require("TensorEmbedding")

n<- 100

weight<- c(2,3)
g <- c( rep(0, weight[1]), rep(1, weight[2]))

p <- length(g)

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
  if(g[i]==0){
    A[,,i]<-  A_mean1 + noise
  }
  if(g[i]==1){
    A[,,i]<-  A_mean2 + noise
  }
}

propMat<- 1/(1+exp(-A))

simA <- (array(runif(n*n*p), dim = c(n,n,p)))

for(i in 1:p){
  simA[,,i]<-  (simA[,,i] + t(simA[,,i]))/2
}

simA <- (simA < propMat)*1


Aavg = array(0, dim = c(n,n,2))


for(i in 1:p){
  if(g[i]==0){
    Aavg[,,1] <- Aavg[,,1]+  simA[,,i] / weight[1]
  }
  if(g[i]==1){
    A[,,i]<-  A_mean2 + noise
    Aavg[,,2] <- Aavg[,,2]+ simA[,,i] / weight[2]
  }
}


testObj<- TensorEmbedding::symm_group_tensor_decomp(simA, g, n, p, k, 10, 1E-4, 1E-4, tol = 1E-12,restrictCoreToDiag = FALSE)

L <-  testObj$L
C <- testObj$C




loglogit<- function(x, theta){
   sum( -x*theta + log(1 + exp(theta)))
}

loss<- function(L,C){
    l<-0
    for(i in 1:2){
        l<- l+ weight[i]* loglogit(Aavg[,,i],   L %*%C[,,i]%*% t(L))
    }
    return (l + sum(L*L/1E4)/2 +sum(C*C/1E4)/2)
}

loss(L,C)

numGradL<- function(i,j){

    L1<- L
    L1[i,j]<- L[i,j]+1E-8

    (loss(L1,C)- loss(L,C))/1E-8

}



numGradC<- function(i,j, k){

    C1<- C
    C1[i,j,k]<- C[i,j,k]+1E-8
    C1[j,i,k]<- C[i,j,k]+1E-8

    (loss(L,C1)- loss(L,C))/1E-8

}

testObj$gradL[55,3]
numGradL(55,3)

