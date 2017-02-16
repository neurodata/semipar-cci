require("rstiefel")
require("msm")

symmetrize<- function(X){
  X= t(X)
  X[lower.tri(X)]<- t(X)[lower.tri(X)]
  X
}

stiefel_diagonalize<- function(Y_list, r){
  
  n = nrow(Y_list[[1]])
  m = length(Y_list)
  
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
  
  # k= 10
  tau = 100
  
  
  avgZ = matrix(rowMeans(matrix(unlist(Y_list),n*n)),n)
  
  # Z_stacked = do.call("cbind",Z_list)
  
  diag(avgZ)=0.99
  avgZ[avgZ==1]=0.99
  avgZ[avgZ==0]=0.01
  svd1<- svd(qnorm(avgZ),nu = r,nv = r)
  U0<- svd1$u
  
  # U_list = list()
  # for(j in 1:m){
  #   U_list[[j]] = rbmf.matrix.gibbs(zeroA, zeroB, k*U0, U0)
  # }
  
  Z_list = list()
  for(j in 1:m){
    Z_list[[j]] = Y_list[[j]] - 0.5
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
    
    # if(FALSE){
    #   #update U0
    #   U_sum = matrix(rowSums(matrix(unlist(U_list),n*r)),n)
    #   
    #   C= k * U_sum
    #   svdC = svd(C)
    #   U0 = rbmf.matrix.gibbs(zeroA, zeroB, svdC$u%*%diag(svdC$d), U0%*%svdC$v) %*% t(svdC$v)
    #   # U0 = rbmf.matrix.gibbs(zeroA, zeroB, k * U_sum, U0)
    #   
    #   #update U
    #   for(j in 1:m){
    #     U_list[[j]] = rbmf.matrix.gibbs(Z_list[[j]]/2, diag(D_list[[j]]), k*U0, U_list[[j]])
    #   }
    # }
    
    trace_D = rbind(trace_D, unlist(D_list))
    print(i)
  }
 
  return(list("D_list"=D_list, "trace_D"=trace_D)) 
}
