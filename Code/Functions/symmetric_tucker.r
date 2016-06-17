require("tensr")

n<- 3
p<- 3

k<- 2 # reduced dimension

A <- array(rnorm(n*n*p), dim = c(n,n,p))

for(i in 1:p){
  A0<- A[,,i]
  A[,,i]<-  (A0+t(A0))/2
}

L<- matrix(rnorm(n*k),n,k)

#core
C <- array(rnorm(k*k*p), dim = c(k,k,p))

for(i in 1:p){
  C0<- C[,,i]
  C[,,i]<-  (C0+t(C0))/2
}

sqr_loss<- function(L,C){
  l<-0
  for(i in 1:p){
    LCL<-L%*%C[,,i]%*%t(L)
    l<- l+sum((A[,,i]-LCL)^2)
  }
  l
}

sqr_loss(L,C)

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
        grad[i,j,m]<- sum(2*diff* (L%*%K%*%t(L)))
      }
    }
  }
  grad
}


grad_C(L,C)


grad_L_num<-function(L,C, i, j){
    l1<- sqr_loss(L,C)
    L0<- L
    L0[i,j] <-L0[i,j]+1E-8
    l2<- sqr_loss(L0,C)
    (l2-l1)/1E-8
}

grad_C_num<-function(L,C, i, j,m){
  l1<- sqr_loss(L,C)
  C0<- C
  C0[i,j,m] <-C0[i,j,m]+1E-8
  l2<- sqr_loss(L,C0)
  (l2-l1)/1E-8
}

grad_L_num(L,C,2,2)
grad_L(L,C)[2,2]


grad_C_num(L,C,2,2,2)
(grad_C(L,C))[2,2,2]
