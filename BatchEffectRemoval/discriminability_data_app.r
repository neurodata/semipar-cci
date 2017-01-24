setwd("~/git/semipar-cci/BatchEffectRemoval/")

loadND<- function(path){
  load(path)
  dataList
}

BNU1<- loadND("../Data/processed/BNU1.RDa")
KKI2009<- loadND("../Data/processed/KKI2009.RDa")
MRN114<- loadND("../Data/processed/MRN114.RDa")

length(BNU1)
length(MRN114)


A_list<-list()

A_list[[1]]<- lapply( BNU1, function(x){x$A})
A_list[[2]]<- lapply( KKI2009, function(x){x$A})
A_list[[3]]<- lapply( MRN114, function(x){x$A})

label_list<- list()
label_list[[1]]<- unlist(lapply( BNU1, function(x){x$SEX}))
label_list[[2]]<- unlist(lapply( KKI2009, function(x){x$SEX}))
label_list[[3]]<- unlist(lapply( MRN114, function(x){x$SEX}))

source("batchremoval.r")

r=10
testRun<- runBatchRemoval(A_list, r=r, 200)

F_list =  testRun$MAP$F_list
C_list =  testRun$MAP$C_list
F = testRun$MAP$F
sigma2 = testRun$MAP$sigma2
trace_sigma2 = testRun$trace_sigma2

# ts.plot(trace_sigma2)

sqrt(sigma2)/sd(F)


image(cbind(C_list[[1]],C_list[[2]],C_list[[3]]))



image(F)
image(F_list[[1]])
image(F_list[[2]])
image(F_list[[3]])
m_j<-3

func_logit<- function(x){
  p<- 1/(1+exp(-x))
  p[p==1]<- 1-1E-6
  p[p==0]<- 1E-6
  p
  }



p_est_list<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list[[j]])), function(x){ 
    F_local = F_list[[j]]
    p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
    p
  } )
})


p_list_wo_batch<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list[[j]])), function(x){ 
    F_local = F
    p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
    p
  } )
})


image(A_list[[2]][[2]],zlim = c(0,1))
image(p_est_list[[2]][[2]],zlim = c(0,1))
image(p_list_wo_batch[[2]][[2]],zlim = c(0,1))


B1<- (F_list[[1]] -F)
FCB = F%*%diag(C_list[[1]][,1])%*%t(B1)
BCB = B1%*%diag(C_list[[1]][,1])%*%t(B1)

overallB = FCB+t(FCB)+BCB
hist(overallB)
image(overallB)

require("mclust")




M1<- kmeans(t(C_list[[1]]*apply(F,MARGIN = 2,sd)),2)$cluster
M2<- kmeans(t(C_list[[2]]*apply(F,MARGIN = 2,sd)),2)$cluster
M3<- kmeans(t(C_list[[3]]*apply(F,MARGIN = 2,sd)),2)$cluster

adjustedRandIndex(M1,label_list[[1]])
adjustedRandIndex(M2,label_list[[2]])
adjustedRandIndex(M3,label_list[[3]])

require("ggplot2")

ns<-  (unlist(lapply(A_list, length)))

batchID<- c(rep(1,ns[[1]]),rep(2,ns[[2]]),rep(3,ns[[3]]))
label<- unlist(label_list)


df = data.frame("value"=unlist(C_list)*apply(F,MARGIN = 2,sd),"index"=c(1:r),"subject"=rep(c(1:sum(ns)),each=r),"batch"=rep(batchID,each=r),"label"=as.factor(rep(label,each=r)))

# df = data.frame("value"=unlist(C_list),"index"=c(1:r),"subject"=rep(c(1:sum(ns)),each=r),"batch"=rep(batchID,each=r),"label"=as.factor(rep(label,each=r)))


pdf("batch_removal_data.pdf",10,6)
ggplot(data=df, aes(x=index,y=value))+geom_line(aes(group=subject,col=label))+ facet_grid(~ batch)
dev.off()

save(testRun,file="resultDataA.RDa")


A_list_flat<- c(A_list[[1]],A_list[[2]],A_list[[3]])
A_list_new<- list()
A_list_new[[1]]<- A_list_flat
testRun2<- runBatchRemoval(A_list_new, r=r, 200)

save(testRun2,file="resultDataANoBE.RDa")




C_list =  testRun2$MAP$C_list
F =  testRun2$MAP$F

df = data.frame("value"=unlist(C_list)*apply(F,MARGIN = 2,sd),"index"=c(1:r),"subject"=rep(c(1:sum(ns)),each=r),"batch"=rep(batchID,each=r),"label"=as.factor(rep(label,each=r)))

pdf("BE_not_removed_data.pdf",10,6)
ggplot(data=df, aes(x=index,y=value))+geom_line(aes(group=subject,col=label))+ facet_grid(~ batch)
dev.off()

X<-cbind(1,t(matrix(unlist(testRun2$C_list),r)))

fit2<- logit(label-1, X,burn = 1000,samp = 1000)
library(pROC)
pred_p2<- c(func_logit(X%*%colMeans(fit2$beta)))

roc2<-roc(label,pred_p2)


X<- cbind(1,t(matrix(unlist(testRun$C_list),r)))

fit1<- logit(label-1, X,burn = 1000,samp = 1000)

pred_p1<- c(func_logit(X%*%colMeans(fit1$beta)))

roc1<-roc(label,pred_p1)

plot(roc2)
lines(roc1,col="red")

auc(roc2)
auc(roc1)
