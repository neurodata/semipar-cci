require("rstiefel")


plot(rmovMF(100,c(10,10)))


n<- 70
r<- 5
U<- matrix(rnorm(n*r),n)
D<- rnorm(r)
V<- matrix(rnorm(n*r),n)

Y<- U%*% diag(D)%*% t(V) + rnorm(n^2)

svdY<- svd(Y,nu = r,nv = 5)

U<- svdY$u
D<- svdY$d[1:r]
V<- svdY$v


B<- matrix(0,r,r)
A<- matrix(0,n,n)


trace_u<-numeric()
for(i in 1:100){
  U = rbmf.matrix.gibbs(A, B, Y%*%V %*% diag(D), U)
  V = rbmf.matrix.gibbs(A, B, t(Y)%*%U %*% diag(D), V)
  trace_u<- c(trace_u, U[1,1])
}


plot(c(U%*%diag(D)%*%t(V)),c(Y))


C = U*10

trace_u<-numeric()
U_list=list()
for(i in 1:1000){
  U<- rbmf.matrix.gibbs(A, B, C, U)
  U_list[[i]] = U
  trace_u<- c(trace_u, U[1,1])
}

acf(trace_u)

plot(U_list[[1]], C)
plot(U_list[[1]], C)


U_sum = matrix(rowSums(matrix(unlist(U_list),n*r)),n)

Cp<-U
for(j in 1:200){
  Cp<- rbmf.matrix.gibbs(A, B, U_sum, Cp)
}


plot(U_sum/1000, Cp)
plot(U_sum, C)
plot(Cp,C)
