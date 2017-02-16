require("rstiefel")
require("msm")

n<- 50
r<- 10
m<- 10

zeroB<- matrix(0,r,r)
zeroA<- matrix(0,n,n)

U<- matrix(rnorm(n*r),n)
svdU = svd(U)
U= svdU$u

symmetrize<- function(X){
  X= t(X)
  X[lower.tri(X)]<- t(X)[lower.tri(X)]
  X
}

Y_list = list()
Z_list = list()
for(j in 1:m){
  D<- c(1:r)*6+ rnorm(r)
  if(j> (m/2))
    D<- 10/c(1:r)+ rnorm(r)
  UDU<- U%*% diag(D)%*% t(U)
  Z_list[[j]]<- UDU + symmetrize(matrix(rnorm(n*n),n))
  Y_list[[j]]<- (Z_list[[j]]>0)*1
  diag(Y_list[[j]])=2
}
rm(U)

ub_list = list()
lb_list = list()
for(j in 1:m){
  mat_ub = matrix(Inf,n,n) 
  mat_lb = matrix(-Inf,n,n) 
  
  mat_lb[Y_list[[j]] ==1] =0
  mat_ub[Y_list[[j]] ==0] =0
  ub_list[[j]]=mat_ub
  lb_list[[j]]=mat_lb
}

k= 10
tau = 100


avgZ = matrix(rowMeans(matrix(unlist(Y_list),n*n)),n)

# Z_stacked = do.call("cbind",Z_list)

diag(avgZ)=0.99
avgZ[avgZ==1]=0.99
avgZ[avgZ==0]=0.01
r=20
svd1<- svd(qnorm(avgZ),nu = r,nv = r)
U0<- svd1$u

U_list = list()
for(j in 1:m){
  U_list[[j]] = rbmf.matrix.gibbs(zeroA, zeroB, k*U0, U0)
}


trace_D<- numeric()
for(i in 1:200){
  
  #update D
  D_list = list()
  for(j in 1:m){
    # vinvm = apply(U_list[[j]], 2, function(x){  sum(x* (Z_list[[j]]%*%x))/2})
    vinvm = apply(U0, 2, function(x){  sum(x* (Z_list[[j]]%*%x))/2})
    v = 2* tau/(2+tau)
    mean_D = v* vinvm
    D_list[[j]] = rnorm(r, mean_D, sqrt(v))
  }
  
  #update Z
  
  Z_list = list()
  for(j in 1:m){
    matZ = matrix(0,n,n)
    mean_Z = U0%*% diag(D_list[[j]])%*% t(U0)
    
    lowertri_idx = lower.tri(mean_Z)
    mean_Z_lowertri = mean_Z[lowertri_idx]
    lb_lowertri = lb_list[[j]][lowertri_idx]
    ub_lowertri = ub_list[[j]][lowertri_idx]
    z = rtnorm(length(mean_Z_lowertri),mean = mean_Z_lowertri,sd = 1,lb_lowertri,ub_lowertri)
    
    matZ[lowertri_idx] = z
    matZ = symmetrize(matZ)
    diag(matZ) = rnorm(n, diag(mean_Z),1)
    Z_list[[j]] = matZ
  }
  
  if(FALSE){
    #update U0
    U_sum = matrix(rowSums(matrix(unlist(U_list),n*r)),n)
    
    C= k * U_sum
    svdC = svd(C)
    U0 = rbmf.matrix.gibbs(zeroA, zeroB, svdC$u%*%diag(svdC$d), U0%*%svdC$v) %*% t(svdC$v)
    # U0 = rbmf.matrix.gibbs(zeroA, zeroB, k * U_sum, U0)
    
    #update U
    for(j in 1:m){
      U_list[[j]] = rbmf.matrix.gibbs(Z_list[[j]]/2, diag(D_list[[j]]), k*U0, U_list[[j]])
    }
  }
  
  trace_D = rbind(trace_D, unlist(D_list))
  print(i)
}


acf(trace_D[,1],lag.max = 40)


# cor(c(U0),c( U_list[[1]]))
# plot(U0, c(U_list[[1]]))
# ts.plot(trace_D)


D<- matrix(unlist(D_list),r)
g<- (c(1:m)> (m/2))*1+1
plot(c(1:r), seq(min(D),max(D),length.out = r),type = "n")
for(j in 1:m){
  lines(c(1:r),D[,j],col=g[j])
}

# plot(colMeans(trace_D))
