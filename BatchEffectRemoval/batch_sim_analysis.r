setwd("~/git/semipar-cci/BatchEffectRemoval/")


load(file="varyingSigma.RDa")
load(file="varyingM.RDa")
load(file="varyingN.RDa")
load(file="varyingRank.RDa")


# length(varyingSigma)
# 
# image(varyingSigma[[3]]$batchRemoval$C_list[[2]])
# image(varyingSigma[[3]]$batchNotRemoval$C_list[[1]])
# 
# 
# varyingSigma[[3]]$batchRemoval$C_list[[2]]
# varyingSigma[[3]]$batchNotRemoval$C_list[[1]][1,]
# 
# varyingSigma[[3]]$batchRemoval$C_list[[1]][1,]
# varyingSigma[[3]]$batchRemoval$C_list[[2]][1,]

require("transport")

computeDist <- function(data){
  
  # m<- length(data[[1]]$A_list[[1]])
  
  dist<- sapply(data, function(x){
    
    m<-length( x$A_list[[1]])
    
    C1<- x$batchRemoval$MAP$C_list[[1]][,1:(m/2)]
    C2<- x$batchRemoval$MAP$C_list[[2]][,1:(m/2)]
    F<- x$batchRemoval$F
    if(!is.null(F)){
      C1<- C1*apply(F,2,sd)
      C2<- C2*apply(F,2,sd)
    }
    
    d1<- tryCatch( wasserstein(pp(t(C1)),pp(t(C2)),prob=F,p=1),error = function(e) NULL)
    # d1 <-  sqrt(mean((rowMeans(C1) - rowMeans(C2))^2))
    # d1<- tryCatch( sqrt(mean((rowMeans(C1) - rowMeans(C2))^2)),error = function(e) NULL)
    
    
    C1<- x$batchNotRemoval$MAP$C_list[[1]][, 1:(m/2)]
    C2<- x$batchNotRemoval$MAP$C_list[[1]][, (m+1):(m+ m/2)]
    F<- x$batchNotRemoval$F
    if(!is.null(F)){
      C1<- C1*apply(F,2,sd)
      C2<- C2*apply(F,2,sd)
    }
    
    # d2<- tryCatch( sqrt(mean((rowMeans(C1) - rowMeans(C2))^2)),error = function(e) NULL)
    d2<- tryCatch( wasserstein(pp(t(C1)),pp(t(C2)),prob=F,p=1),error = function(e) NULL)
    
    c(d1,d2)
  })
  if(is.list(dist))
    dist = do.call("cbind",dist)
  dist
}


require("ggplot2")

distSigma<- computeDist(varyingSigma)
distM<- computeDist(varyingM)
distN<- computeDist(varyingN)
distRank<- computeDist(varyingRank)

# save(distSigma, file="distSigma.RDa")
# save(distM, file="distM.RDa")
# save(distN, file="distN.RDa")

load( file="distSigma.RDa")
load( file="distM.RDa")

pdf("sigma_dist.pdf",6,4)
df<-data.frame("sigma"=seq(0.05,1,length.out = 20),"Distance" = c(distSigma[1,],distSigma[2,]),"method"=rep(c("random factor model","shared factor model"),each=20))
ggplot(data=df, aes(x=sigma,y=Distance))+geom_line(aes(group=method,col=method))
dev.off()




pdf("N_dist.pdf",6,4)
df<-data.frame("n"=seq(10,200,by = 10),"Distance" = c(distN[1,],distN[2,]),"method"=rep(c("random factor model","shared factor model"),each=20))
ggplot(data=df, aes(x=n,y=Distance))+geom_line(aes(group=method,col=method))
dev.off()

pdf("M_dist.pdf",6,4)
df<-data.frame("m"=seq(10,200,by = 10),"Distance" = c(distM[1,],distM[2,]),"method"=rep(c("random factor model","shared factor model"),each=20))
ggplot(data=df, aes(x=m,y=Distance))+geom_line(aes(group=method,col=method))
dev.off()





require("ggplot2")

plotCore<- function(x){
  C = x$MAP$C_list
  F = x$MAP$F
  if(length(C)==2)
    m<- ncol(C[[1]])
  else
    m<- ncol(C[[1]])/2
  r<- nrow(C[[1]])
  batchID = rep(1:2, each=m)
  label = rep(1:2, each= m/2)
  
  df = data.frame("value"=unlist(C)*apply(F,MARGIN = 2,sd),"index"=c(1:r),"subject"=rep(c(1:m),each=r),"batch"=rep(batchID,each=r),"label"=as.factor(rep(label,each=r)))
  ggplot(data=df, aes(x=index,y=value))+geom_line(aes(group=subject,col=label))+ facet_grid(~ batch)
}

plotCore(varyingSigma[[10]]$batchRemoval)
plotCore(varyingSigma[[10]]$batchNotRemoval)

require("class")



C_combined1<- t(matrix(unlist(varyingSigma[[10]]$batchNotRemoval$C_list),5))
C_combined2<- t(matrix(unlist(varyingSigma[[10]]$batchRemoval$C_list),5))

KNNtest<-function(k, C_combined){
  m<- nrow(C_combined)
  
  label<- rep(c(1,2,1,2),each= m/4)
  
  cl<- sapply(1:length(label), function(i){
    true_cl = label[i]
    true_cl==knn(C_combined[-i,],test = C_combined[i,],cl=label[-i],k = k)
  })
  1-sum(cl)/m
}

KNNtest(1, C_combined1)
KNNtest(1, C_combined2)
