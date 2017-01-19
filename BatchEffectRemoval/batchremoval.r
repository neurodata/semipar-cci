runBatchRemoval<- function(A_list, r, tot_iter){
  
  trace_sigma2<- numeric()
  
  samplePsi= function(X, w, K, Binv, b){
    r<- ncol(X)
    V<- solve(t(X)%*% (X*w) + Binv)
    m<- V%*% ( t(X)%*%K+ Binv%*%b)
    t(chol(V))%*% rnorm(r)+m
  }
  
  
  m_j = length(A_list)
  n = nrow(A_list[[1]][[1]])
  
  #initialize F and C with SVD
  
  avgA = matrix(rowMeans(matrix(unlist(A_list),n^2)),n)
  avgA[avgA==1]<- 1- 1E-6
  avgA[avgA==0]<- 1E-6
  avgP <- log(avgA)/(log(1-avgA))
  eigenLogitAvgA= svd(avgP,nu = r,nv = r)
  avgC = eigenLogitAvgA$d[1:r]
  avgF = eigenLogitAvgA$u
  
  F = avgF
  
  # randomly initialize F
  F_list<- list()
  for(j in 1:m_j){
    F_list[[j]]<-  avgF + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)
  }
  # randomly initialize C
  C_list<- list()
  for(j in 1:m_j){
    C_list[[j]] <- avgC + matrix( rnorm(r*length(A_list[[j]]), sd = 0.1), r)
  }
  
  for(iter in 1:tot_iter){
    #update w_list
    w_list<- lapply(c(1:m_j),function(j){
      lapply( c(1:ncol(C_list[[j]])), function(x){ 
        
        F_local = F_list[[j]]
        psi<- (F_local%*% (t(F_local)*C_list[[j]][,x]))
        w<- rpg(n*(n-1)/2,1,psi[upper.tri(psi)])
        W<- matrix(0,n,n)
        W[upper.tri(W)]<- w
        W[lower.tri(W)]<- t(W)[lower.tri(W)]
        W
      } )
    })
    
    #update F_list
    for(j in 1:m_j){
      for(k in 1:n){
        
        n_j<- ncol(C_list[[j]])
        #update F
        F_wo_k<-  F_list[[j]][-k,]
        
        y<- numeric()
        w<- numeric()
        X<- numeric()
        for(i in 1:n_j){
          y_i<- A_list[[j]][[i]][k,-k]
          w_i<- w_list[[j]][[i]][k,-k]
          X_i<- F_wo_k %*% diag(C_list[[j]][,i]) 
          y<- c(y,y_i)
          w<- c(w,w_i)
          X<- rbind(X,X_i)
        }
        
        K<- y- 1/2
        Binv<- diag(1/sigma2, r)
        b<- F[k,]
        F_list[[j]][k,]<- samplePsi(X, w, K, Binv,b)
      }
    }
    
    #update C_list
    
    for(j in 1:m_j){
      n_j<- ncol(C_list[[j]])
      
      for(i in 1:n_j){
        
        
        y<- numeric()
        w<- numeric()
        X<- numeric()
        
        for(k in 1: (n-1)){
          X_i <- F[-(1:k),] %*%diag( F[k,]) 
          y_i <- A_list[[j]][[i]][k,-(1:k)]
          w_i<- w_list[[j]][[i]][k,-(1:k)]
          y<- c(y,y_i)
          w<- c(w,w_i)
          X<- rbind(X,X_i)
        }
        
        K<- y- 1/2
        Binv<- diag(1, r)
        b<- rep(0,r)
        C_list[[j]][,i]<- samplePsi(X, w, K, Binv,b)
      }
    }
    #update F
    
    F_sum<- F_list[[1]]
    for(j in 2:length(F_list)){
      F_sum = F_sum + F_list[[j]]
    }
    
    F<- F_sum/length(F_list) + rnorm( n*r, sd=sqrt(sigma2/length(F_list)))
    
    #update sigma
    
    F_diff2_sum<- (F_list[[1]]-F)^2
    for(j in 2:length(F_list)){
      F_diff2_sum =   F_diff2_sum + (F_list[[j]]-F)^2
    }
    
    sigma2 = 1/rgamma(1, n*r*m_j/2 , sum(F_diff2_sum)/2)
    
    trace_sigma2<- c(trace_sigma2, sigma2)
    
    print(iter)
  }
  
  return(list("F_list"=F_list,
              "C_list"=C_list,
              "F"=F,
              "sigma2"= sigma2,
              "trace_sigma2"= trace_sigma2
              ))
  
}
