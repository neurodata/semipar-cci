require("BayesLogit")

runBatchRemoval<- function(A_list, r, tot_iter, alpha =1, q=1.5){
  
  m_j = length(A_list)
  
  trace_sigma2<- numeric()
  
  sigma2<- rep(0.1,m_j)
  
  
  #gamma for var of C
  

  gamma_a= alpha * q^(3* ((1:r)-1))
  gamma_b=  q^(2* ((1:r)-1))
  var_c = 1/rgamma(r, gamma_a, rate=gamma_b)
  
  # func_logit<-  function(x){1/(1+exp(-x))}
  func_loglogit<-  function(x){
    l = -log(1+exp(-x))
    l[l< -1E6]<- -1E6
    l[l> 1E6]<- 1E6
    l
  }
  
  samplePsi= function(X, w, K, Binv, b){
    r<- ncol(X)
    V<- solve(t(X)%*% (X*w) + Binv)
    m<- V%*% ( t(X)%*%K+ Binv%*%b)
    t(chol(V))%*% rnorm(r)+m
  }
  
  symmetrize<- function(X){
    X= t(X)
    X[lower.tri(X)]<- t(X)[lower.tri(X)]
    X
  }
  
  
  computeLoss = function(F, F_list,C_list, Z_list,psi_list,sigma2,nu2){
    
    loss = 0
    
    loss = loss - sum(F^2)/2/nu2
    
    for(j in 1:m_j){
      
      F_local = F_list[[j]]
      
      loss =  loss - sum((F_local- F)^2)/2/sigma2[j] - n*r*log(sigma2[j])/2
      
      Z_local = Z_list[[j]]
      
      for(x in 1:ncol(C_list[[j]])){ 
        
        A_local = A_list[[j]][[x]]
        psi<- psi_list[[j]][[x]]
        
        
        loss = loss + sum(A_local * func_loglogit(psi + Z_local) + (1-A_local)* func_loglogit(- psi - Z_local))
        
      }
    }
    
    loss
    
  }
  
  
  n = nrow(A_list[[1]][[1]])
  
  #initialize F and C with SVD
  
  avgA = matrix(rowMeans(matrix(unlist(A_list),n^2)),n)
  avgA[avgA==1]<- 1- 1E-6
  avgA[avgA==0]<- 1E-6
  avgP <- log(avgA)/(log(1-avgA))
  eigenLogitAvgA= svd(avgP,nu = r,nv = r)
  avgC = eigenLogitAvgA$d[1:r]
  avgF = eigenLogitAvgA$u
  
  avgZ = avgP - avgF%*% diag(avgC) %*%t(avgF)
  
  F = avgF
  
  # randomly initialize F
  F_list<- list()
  for(j in 1:m_j){
    F_list[[j]]<-  avgF + matrix(rnorm(n*r,sd=sqrt(sigma2[j])),n)
  }
  # randomly initialize C
  C_list<- list()
  for(j in 1:m_j){
    C_list[[j]] <- avgC + matrix( rnorm(r*length(A_list[[j]]), sd = 0.1), r)
  }
  #randomly initialize Z
  Z_list<- list()
  for(j in 1:m_j){
    Z_list[[j]]<-  matrix(0,n,n) #0 + symmetrize(matrix(rnorm(n*n,sd=1),n))
  }
  
  
  best_loss = -Inf
  MAP_est = list()
  
  #Compute psi = FCF
  
  psi_list<- lapply(c(1:m_j),function(j){
    lapply( c(1:ncol(C_list[[j]])), function(x){ 
      
      F_local = F_list[[j]]
      psi<- (F_local%*% (t(F_local)*C_list[[j]][,x]))
      psi
    } )
  })
  
  for(iter in 1:tot_iter){
    
    
    
    #update w_list
    w_list<- lapply(c(1:m_j),function(j){
      lapply( c(1:ncol(C_list[[j]])), function(x){ 
        
        psi0<- psi_list[[j]][[x]] + Z_list[[j]]
        b = psi0[lower.tri(psi0)]
        w<- rpg(n*(n-1)/2,1,b)
        
        W<- matrix(0,n,n)
        W[lower.tri(W)]<- w
        W= symmetrize(W)
        W
      } )
    })
    
    #update Z_list
    for(j in 1:m_j){
      
      lowerWs = sapply(w_list[[j]], function(x) x[lower.tri(x)])
      lowerPsi = sapply(psi_list[[j]], function(x) x[lower.tri(x)])
      lowerYs = sapply(A_list[[j]], function(x) x[lower.tri(x)])
      
      v = 1/(rowSums(lowerWs) + 1E-2)
      m = v * rowSums( - lowerWs* lowerPsi + lowerYs -0.5)
      
      z = rnorm(length(m),m, sqrt(v))
      
      Z_list[[j]] [lower.tri(Z_list[[j]])] = z
      Z_list[[j]] = symmetrize(Z_list[[j]])
    }
    
    #update F_list
    for(j in 1:m_j){
      for(k in 1:n){
        
        p<- ncol(C_list[[j]])
        #update F
        F_wo_k<-  F_list[[j]][-k,]
        
        # y<- numeric()
        # w<- numeric()
        # X<- numeric()
        # for(i in 1:n_j){
        #   y_i<- A_list[[j]][[i]][k,-k]
        #   w_i<- w_list[[j]][[i]][k,-k]
        #   X_i<- F_wo_k %*% diag(C_list[[j]][,i]) 
        #   y<- c(y,y_i)
        #   w<- c(w,w_i)
        #   X<- rbind(X,X_i)
        # }
        
        y<- c(sapply(c(1:p), function(i) A_list[[j]][[i]][k,-k]))
        z<- c(sapply(c(1:p), function(i) Z_list[[j]][k,-k]))
        w<- c(sapply(c(1:p), function(i) w_list[[j]][[i]][k,-k]))
        X<- matrix(sapply(c(1:p), function(i) c(C_list[[j]][,i]) * t(F_wo_k)),ncol=r,byrow = T)
        # y<- numeric()
        # w<- numeric()
        # X<- numeric()
        # for(i in 1:p){
        #   y_i<- A_list[[j]][[i]][k,-k]
        #   w_i<- w_list[[j]][[i]][k,-k]
        #   X_i<- F_wo_k %*% diag(C_list[[j]][,i]) 
        #   y<- c(y,y_i)
        #   w<- c(w,w_i)
        #   X<- rbind(X,X_i)
        # }

        K<- y- 1/2 - w*z
        Binv<- diag(1/sigma2[j], r)
        b<- F[k,]
        F_list[[j]][k,]<- samplePsi(X, w, K, Binv,b)
      }
    }
    
    #update C_list
    
    for(j in 1:m_j){
      p<- ncol(C_list[[j]])
      
      if(iter <10){
        F_local=F
      }      else{
        F_local = F_list[[j]]
      }      
      
      X<- matrix( unlist(sapply(c(1:(n-1)), function(k) (t(F_local[-(1:k),])* ( F_local[k,])))), ncol=r,byrow = T)
      
      for(i in 1:p){
        
        y<- unlist(sapply(c(1:(n-1)), function(k) A_list[[j]][[i]][k,-(1:k)]))
        w<- unlist(sapply(c(1:(n-1)), function(k) w_list[[j]][[i]][k,-(1:k)]))
        
        z<- unlist(sapply(c(1:(n-1)), function(k) Z_list[[j]][k,-(1:k)]))
        

        K<- y- 1/2 - w *z
        Binv<- diag(1/var_c)#diag(1E-5, r)
        b<- rep(0,r)
        C_list[[j]][,i]<- samplePsi(X, w, K, Binv,b)
      }
    }
    #update F with prior (0, nu2)
    
    nu2 = 1
    
    F_sum<- F_list[[1]]/sigma2[1]
    for(j in 2:length(F_list)){
      if(m_j>1)
        F_sum = F_sum + F_list[[j]]/sigma2[j]
    }
    
    var = 1/( sum(1/sigma2) + 1/nu2)
    m =  var* (F_sum)
    

    F<- m + rnorm( n*r, sd=sqrt(var))
    
    #Compute psi = FCF
    
    psi_list<- lapply(c(1:m_j),function(j){
      lapply( c(1:ncol(C_list[[j]])), function(x){ 
        
        F_local = F_list[[j]]
        psi<- (F_local%*% (t(F_local)*C_list[[j]][,x]))
        psi
      } )
    })
    
    #update sigma
    
    
    for(j in 1:m_j){
      F_diff2_sum<- (F_list[[j]]-F)^2
      # if(EM){
      #   sigma2[j] =  (sum(F_diff2_sum)/2+ 0.01 )/((n*r)/2+0.01)
      # }      else{
        sigma2[j] = 1/rgamma(1, n*r/2-1 , sum(F_diff2_sum)/2)
        # }
    }
    sigma2[sigma2>2] = 2
    
    
    # for(j in 2:length(F_list)){
    #   if(m_j>1)
    #     F_diff2_sum =   F_diff2_sum + (F_list[[j]]-F)^2
    # }
    
    #update var_c
    
    C_mat = do.call("cbind",C_list)
    
    var_c = 1/rgamma(r, gamma_a + ncol(C_mat)/2, rate=gamma_b + rowSums(C_mat^2)/2)
    
    
    # trace_sigma2<- c(trace_sigma2, sigma2)
    
    print(iter)
    print(var_c)
    
    #code to capture and store MAP
    
    est = list("F_list"=F_list,
               "C_list"=C_list,
               "Z_list"=Z_list,
               "F"=F,
               "sigma2"= sigma2)
    
    cur_loss = computeLoss(F, F_list,C_list,Z_list,psi_list,sigma2,nu2)
    
    if(!is.na(cur_loss)){
      if(cur_loss> best_loss){
        MAP_est = est
        best_loss = cur_loss
      }
    }
    
    # if(is.na(cur_loss)){
    # break
    # }
    
    print(c(cur_loss,best_loss))
    print(sigma2)
  }
  
  return(list("F_list"=F_list,
              "C_list"=C_list,
              "Z_list"=Z_list,
              "F"=F,
              "sigma2"= sigma2,
              "trace_sigma2"= trace_sigma2,
              "MAP" = MAP_est
  ))
  
}
