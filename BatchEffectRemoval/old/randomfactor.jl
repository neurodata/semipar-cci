function runBatchRemoval(A_list, r, tot_iter, EM= true)
  
  function func_loglogit(x)
    l = -log(1+exp(-x))
    l[l .< -1E6]= -1E6
    l[l .> 1E6]= 1E6
    l
  end
  
  n_subjects = map(length, A_list)
  
  g = length(A_list)
  n = size(A_list[1][1],1)
  
  sigma2 = 0.5^2
  
  F = zeros(n,r)
  F_list = []
  C_list = []
  
  psi_list= deepcopy(A_list)
  w_list= deepcopy(A_list)
  
  for j = 1:g
    push!( F_list,  zeros(n,r) )
    push!( C_list,  zeros(r, length(A_list[j]))    )
  end
  
  
  function samplePsi(X, w, K, Binv, b)
    r= size(X,2);
    V = inv( X'* (X.* w) + Binv)
    m= V * ( X' * K+ Binv * b)
    if(EM)
      m
    else
      chol(V)' * rand(r) + m
    end
  end
  
  function computeLoss(F, F_list,C_list,psi_list,sigma2,nu2)
    
    loss = 0
    loss  -= sum(F.^2)./2./nu2
    for j=1:g
      local F_local = F_list[j]
      loss =  loss - sum((F_local- F).^2)/2/sigma2 - n*r*log(sigma2)/2
      for x= 1:size(C_list[j],2)
        A_local = A_list[j][x]
        psi= psi_list[j][x]
        loss = loss + sum(A_local .* func_loglogit(psi) + (1-A_local).* func_loglogit(-psi))
      end
    end
    loss
  end
  
  
  #initialize F and C with SVD
  
  
  A_flat=[]
  for  A1 = A_list
    for A2 = A1
      A_flat = vcat(A_flat, vec(A2))
    end
  end
  
  total_m = sum(n_subjects)
  avgA =  reshape( sum(reshape(A_flat, n^2,total_m),2)/total_m, n,n)
  avgA[avgA.==1]= 1- 1E-6
  avgA[avgA.==0]= 1E-6
  avgP = log(avgA)./(log(1-avgA))
  eigenLogitAvgA= svd(avgP)
  
  
  avgC = eigenLogitAvgA[2][1:r]
  avgF = ((eigenLogitAvgA[1] .+ eigenLogitAvgA[3])/2)[:,1:r]
  
  F = deepcopy(avgF)
  
  # randomly initialize F and C
  for j=1:g
    F_list[j]=  avgF + randn(n,r) .* sqrt(sigma2)
    C_list[j] = randn( r, length(A_list[j])) * 0.1 .+ avgC
  end
  
  
  best_loss = -Inf
  MAP_est = []
  
  #Compute psi = FCF
  
  for j = 1:g
    for x = 1: size(C_list[j],2)
      local F_local = F_list[j]
      psi_list[j][x] = F_local *  diagm(C_list[j][:,x]) * F_local'
    end
  end
  
  trace_sigma2 = []
  
  for iter= 1:tot_iter
    
    #update w_list
    
    @parallel for j = 1:g
      for x = 1: size(C_list[j],2)
        psi= psi_list[j][x]
        b = psi
        if(EM)
          w= 1/2./b .* tanh(b./2)
          # else
          # w= rpg(n*(n-1)/2,1,b)
          # W= matrix(0,n,n)
          # W[upper.tri(W)]= w
          # W[lower.tri(W)]= t(W)[lower.tri(W)]
        end
        w_list[j][x] = w
      end
    end
    
    #update F_list
    @parallel for j= 1:g
      for k= 1:n
        
        n_j= size(C_list[j],2)
        #update F
        F_wo_k=  F_list[j][1:end .!= k,:]
        
        y= []
        w= []
        X= zeros(0,size(F_wo_k,2))
        for i in 1:n_j
          y_i= A_list[j][i][k,1:end .!= k]
          w_i= w_list[j][i][k,1:end .!= k]
          X_i= F_wo_k * diagm(C_list[j][:,i])
          y= vcat(y,y_i)
          w= vcat(w,w_i)
          X= vcat(X,X_i)
        end
        
        K= y .- 1/2
        Binv= eye(r) .* 1/sigma2
        b= F[k,:]
        F_list[j][k,:]= samplePsi(X, w, K, Binv,b)
      end
    end
    
    #update C_list
    
    @parallel for j = 1:g
      n_j= size(C_list[j],2)
      
      for i = 1:n_j
        
        # if iter>10
        local F_local = F_list[j]
        # else
        # local F_local = F
        # end
        y= []
        w= []
        X= zeros(0,size(F_local,2))
        
        for k = 1: (n-1)
          X_i = F_local[(k+1):end,:] *diagm( F_local[k,:])
          y_i = A_list[j][i][k,(k+1):end]
          w_i = w_list[j][i][k,(k+1):end]
          y= vcat(y,y_i)
          w= vcat(w,w_i)
          X= vcat(X,X_i)
        end
        
        K= y .- 1/2
        Binv= eye(r) .* 1E-5
        b= zeros(r)
        C_list[j][:,i]= samplePsi(X, w, K, Binv,b)
      end
    end
    #update F with prior (0, nu2)
    
    nu2 = 1
    
    F_sum = deepcopy(F_list[1])
    for j = 2:length(F_list)
      if(g>1)
        F_sum = F_sum + F_list[j]
      end
    end
    
    var = 1/(g /sigma2 + 1/nu2)
    m =  var* (F_sum ./sigma2)
    
    if(EM)
      F= m
      # else
      # F= m + rnorm( n*r, sd=sqrt(var))
    end
    
    #Compute psi = FCF
    @parallel for j = 1:g
      for x = 1: size(C_list[j],2)
        local F_local = F_list[j]
        psi_list[j][x] = F_local *  diagm(C_list[j][:,x]) * F_local'
      end
    end
    
    #update sigma
    
    F_diff2_sum= (F_list[1]-F).^2
    for j =2:length(F_list)
      if(g>1)
        F_diff2_sum =   F_diff2_sum + (F_list[j]-F).^2
      end
    end
    
    if(EM)
      sigma2 =  sum(F_diff2_sum)/(n*r*g)
      # else
      # sigma2 = 1/rgamma(1, n*r*g/2 , sum(F_diff2_sum)/2)
    end
    
    cur_loss = computeLoss(F, F_list,C_list,psi_list,sigma2,nu2)
    
    println(cur_loss)
    
    trace_sigma2 = vcat(trace_sigma2, sigma2)
  end
  
  Dict("F_list"=>F_list,
  "C_list"=>C_list,
  "F"=>F,
  "sigma2"=> sigma2,
  "trace_sigma2"=> trace_sigma2
  )
  
  # "MAP" = MAP_est)
  
end
