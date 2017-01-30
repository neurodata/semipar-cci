setwd("~/git/semipar-cci/BatchEffectRemoval/")

runSim <- function(A_list, r, iter =200){
  
  source("batchremoval.r")
  
  batchRemoval<- tryCatch(runBatchRemoval(A_list, r=r, iter, EM=T),error = function(e) NULL)
  
  A_list_flat<- c(A_list[[1]],A_list[[2]])
  A_list_new<- list()
  A_list_new[[1]]<- A_list_flat
  
  batchNotRemoval<- tryCatch(runBatchRemoval(A_list_new, r=r, iter, EM=T),error = function(e) NULL)
  
  list("A_list" = A_list, "batchRemoval" = batchRemoval , "batchNotRemoval"= batchNotRemoval)
}


A_list0 <- lapply( c(1:10), function(x){
  n=50
  m=50
  r=5
  sigma2=0.5^2
  
  F<- matrix(rnorm(n*r),n)
  
  F_list<- list()
  F_list[[1]]<-  F + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)
  F_list[[2]]<-  F + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)
  
  m_j = 2
  
  C_list<- list()
  C<- matrix( c(rep( 0.1* c(1:r), m/2), rep( 1/ c(1:r)^3 ,m/2)), nrow = r)
  
  C_list[[1]] <- C + matrix( rnorm(r*m, sd =0.01), r)
  C_list[[2]] <- C + matrix( rnorm(r*m, sd =0.01), r)
  C_list0 <- C_list
  
  func_logit<-  function(x){1/(1+exp(-x))}
  func_ilogit<-  function(x){log(x/(1-log(x)))}
  
  p_list<- lapply(c(1:m_j),function(j){
    lapply( c(1:ncol(C_list[[j]])), function(x){ 
      
      F_local = F_list[[j]]
      p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
      p
    } )
  })
  
  A_list<- lapply(c(1:m_j),function(j){
    lapply( c(1:ncol(C_list[[j]])), function(x){ 
      
      p<- p_list[[j]][[x]]
      A<- (p> runif(length(p)))*1
      A[lower.tri(A)]<- t(A)[lower.tri(A)]
      A
    } )
  })
  
  rm(C_list)
  rm(F_list)
  A_list
})



require("snow")


varyingRank<- list()

for(i in c(1:10)){
  A = A_list0[[i]]
  clus <- makeCluster(32)
  clusterExport(clus,"runSim")
  clusterExport(clus,"A")
  varyingRank[[i]]<- parLapply(clus,seq(2,50,length.out = 25), function(r){runSim(A,r,200)})
  stopCluster(clus)
}

save(varyingRank,file="varyingRank.RDa")



