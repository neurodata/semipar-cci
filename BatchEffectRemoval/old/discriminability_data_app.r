setwd("~/git/semipar-cci/BatchEffectRemoval/")

loadND<- function(path){
  load(path)
  dataList
}


func_logit<- function(x){
  p<- 1/(1+exp(-x))
  p[p==1]<- 1-1E-6
  p[p==0]<- 1E-6
  p
}


BNU1<- loadND("../Data/processed/BNU1.RDa")
BNU3<- loadND("../Data/processed/BNU3.RDa")
KKI2009<- loadND("../Data/processed/KKI2009.RDa")
MRN114<- loadND("../Data/processed/MRN114.RDa")
SWU4<- loadND("../Data/processed/SWU4.RDa")

A_list<-list()

A_list[[1]]<- lapply( BNU1, function(x){x$A})
A_list[[2]]<- lapply( BNU3, function(x){x$A})
A_list[[3]]<- lapply( KKI2009, function(x){x$A})
A_list[[4]]<- lapply( MRN114, function(x){x$A})
A_list[[5]]<- lapply( SWU4, function(x){x$A})

label_list<- list()
label_list[[1]]<- unlist(lapply( BNU1, function(x){x$SEX}))
label_list[[2]]<- unlist(lapply( BNU3, function(x){x$SEX}))
label_list[[3]]<- unlist(lapply( KKI2009, function(x){x$SEX}))
label_list[[4]]<- unlist(lapply( MRN114, function(x){x$SEX}))
label_list[[5]]<- unlist(lapply( SWU4, function(x){x$SEX}))


# label_list2<- list()
# label_list2[[1]]<- unlist(lapply( BNU1, function(x){x$SEX}))
# label_list2[[2]]<- unlist(lapply( KKI2009, function(x){x$SEX}))
# label_list2[[3]]<- unlist(lapply( MRN114, function(x){x$SEX}))



#discriminability based on raw data
require("class")

A_list_lowertri<- list()

A_list_lowertri<- lapply(A_list, function(x){
  t(sapply(x,function(x){x[lower.tri(x)]}))
})

KNNtest<-function(k, M,label){
  n<- length(label)
  cl<- sapply(1:length(label), function(i){
    true_cl = label[i]
    true_cl==knn(train = M[-i,],test = M[i,],cl=label[-i],k = k)
  })
  1-sum(cl)/n
}


#get fitted functions
fitted<- function(C_list,F_list,m_j){
  lapply(c(1:m_j),function(j){
    lapply( c(1:ncol(C_list[[j]])), function(x){ 
      F_local = F_list[[j]]
      p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
      p
    } )
  })
}

fitted_adjusted<- function(C_list,F,m_j){
  lapply(c(1:m_j),function(j){
    lapply( c(1:ncol(C_list[[j]])), function(x){ 
      F_local = F
      p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
      p
    } )
  })
}


LDAtest<-function(M,label){
  n<- length(label)
  cl<- sapply(1:length(label), function(i){
    m1 <- colMeans(M[-i,][label[-i]==1,])
    m2 <- colMeans( M[-i,][label[-i]==2,])
    
    if( sum((M[i,]-m1)^2) < sum((M[i,]-m2)^2)){
      choice = 1
    }else{choice =2}
    true_cl = label[i]
    true_cl==choice
  })
  1-sum(cl)/n
}


#individual raw
knn_raw_error<- list()
for(i in 1:length(A_list_lowertri)){
  knn_raw_error[[i]]<- sapply(c(1:20), function(k) KNNtest(k,A_list_lowertri[[i]],label_list[[i]]))
}

# save(knn_raw_error,file="knn_raw_error.RDa")

load("knn_raw_error.RDa")
lda_raw_error<- list()
for(i in 1:length(A_list_lowertri)){
  lda_raw_error[[i]]<-  LDAtest(A_list_lowertri[[i]],label_list[[i]])
}

#center and combine
A_list_lowertri_centered<- lapply(A_list_lowertri, function(x){
  t(t(x) - colMeans(x))
})

A_mat_lowertri_centered<- do.call("rbind",A_list_lowertri_centered)
A_mat_lowertri<- do.call("rbind",A_list_lowertri)

# A_mat_lowertri<-rbind(A_list_lowertri[[1]],A_list_lowertri[[3]],A_list_lowertri[[4]])

#batch removal
load("data_batch_removed.RDa")
load("data_batch_not_removed_rank10.RDa")
# load("data_batch_not_removed_rank20.RDa")

knn_reduced_core_ind_error<- list()
for(i in 1:length(A_list_lowertri)){
  knn_reduced_core_ind_error[[i]]<- sapply(c(1:20), function(k) KNNtest(k,t(batch_removed$MAP$C_list[[i]]* apply(batch_removed$MAP$F,2,sd)),label_list[[i]]))
}


knn_reduced_core_ind_error<- list()
for(i in 1:3){
  knn_reduced_core_ind_error[[i]]<- sapply(c(1:20), function(k) KNNtest(k,t(batch_removed$MAP$C_list[[i]]* apply(batch_removed$MAP$F,2,sd)),label_list2[[i]]))
}

#

shared_factor_C_list<- list()
shared_factor_C_list[[1]]<- batch_not_removed10$MAP$C_list[[1]][,1:ns[1]]
for(i in 2:length(ns)){
  shared_factor_C_list[[i]]<- batch_not_removed10$MAP$C_list[[1]][, (cumsum(ns)[i-1]+1):cumsum(ns)[i]]
}

knn_reduced_shared_factor_core_ind_error<- list()
for(i in 1:length(A_list_lowertri)){
  knn_reduced_shared_factor_core_ind_error[[i]]<- sapply(c(1:20), function(k) KNNtest(k,t(shared_factor_C_list[[i]]* apply(batch_not_removed10$MAP$F,2,sd)),label_list[[i]]))
}





#combine
knn_raw_combined_error<- sapply(c(1:20), function(k) KNNtest(k,A_mat_lowertri,unlist(label_list)))
knn_raw_combined_error_centered<- sapply(c(1:20), function(k) KNNtest(k,A_mat_lowertri_centered,unlist(label_list)))
# save(knn_raw_combined_error,file="knn_raw_combined_error.RDa")
load("knn_raw_combined_error.RDa")

knn_reduced_core_combined_error<-  sapply(c(1:20), function(k) KNNtest(k,t(matrix(unlist(batch_removed$MAP$C_list),10) * apply(batch_removed$MAP$F,2,sd)),unlist(label_list)))

knn_raw_combined_error<- sapply(c(1:20), function(k) KNNtest(k,A_mat_lowertri))


batchAdjustedData<- fitted_adjusted(C_list = batch_removed$C_list, F = batch_removed$F, 5)


p_adjusted_list_lowertri<- lapply(batchAdjustedData, function(x){
  t(sapply(x,function(x){x[lower.tri(x)]}))
})

knn_reduced_p_ind_error<- list()
for(i in 1:length(A_list_lowertri)){
  knn_reduced_p_ind_error[[i]]<- sapply(c(1:20), function(k) KNNtest(k,p_adjusted_list_lowertri[[i]],label_list[[i]]))
}




plot(knn_raw_error[[1]],ylim=c(0,1),type="l")
lines(knn_reduced_core_ind_error[[1]], col="red")
lines(knn_reduced_p_ind_error[[1]], col="blue")
lines(knn_reduced_shared_factor_core_ind_error[[1]],col="green")

plot(knn_raw_error[[2]],ylim=c(0,1),type="l")
lines(knn_reduced_core_ind_error[[2]], col="red")
lines(knn_reduced_p_ind_error[[2]], col="blue")
lines(knn_reduced_shared_factor_core_ind_error[[2]],col="green")

plot(knn_raw_error[[3]],ylim=c(0,1),type="l")
lines(knn_reduced_core_ind_error[[3]], col="red")
lines(knn_reduced_p_ind_error[[3]], col="blue")
lines(knn_reduced_shared_factor_core_ind_error[[3]],col="green")

plot(knn_raw_error[[4]],ylim=c(0,1),type="l")
lines(knn_reduced_core_ind_error[[4]], col="red")
lines(knn_reduced_p_ind_error[[4]], col="blue")
lines(knn_reduced_shared_factor_core_ind_error[[4]],col="green")


plot(knn_raw_error[[5]],ylim=c(0,1),type="l")
lines(knn_reduced_core_ind_error[[5]], col="red")
lines(knn_reduced_p_ind_error[[5]], col="blue")
lines(knn_reduced_shared_factor_core_ind_error[[5]],col="green")

#
hist(batch_removed$MAP$F_list[[1]] - batch_removed$MAP$F)
hist(batch_removed$MAP$F_list[[2]] - batch_removed$MAP$F)
hist(batch_removed$MAP$F_list[[3]]- batch_removed$MAP$F)
hist(batch_removed$MAP$F_list[[4]]- batch_removed$MAP$F)
hist(batch_removed$MAP$F_list[[5]]- batch_removed$MAP$F)

#

p_mat_lowertri<- do.call("rbind",p_adjusted_list_lowertri)
knn_reduced_p_combined_error<-  sapply(c(1:20), function(k) KNNtest(k,p_mat_lowertri,unlist(label_list)))


knn_reduced_core_shared_factor_combined_error<-  sapply(c(1:20), function(k) KNNtest(k,t(batch_not_removed10$MAP$C_list[[1]]),unlist(label_list)))


plot(knn_raw_combined_error,ylim=c(0,1),type="l")
lines(knn_reduced_core_combined_error, col="red")
lines(knn_reduced_p_combined_error,col="blue")
lines(knn_reduced_core_shared_factor_combined_error,col="green")


# lda_random_factor_ind_error<- list()
# for(i in 1:length(A_list_lowertri)){
#   lda_random_factor_ind_error[[i]]<-  LDAtest(t(batch_removed$MAP$C_list[[i]]),label_list[[i]])
# }
# 
# plot(unlist(lda_raw_error))
# lines(unlist(lda_random_factor_ind_error))




# 
# knn_reduced_core_shared_factor<- sapply(c(1:20), function(k) KNNtest(k,t(matrix(unlist(batch_not_removed10$MAP$C_list),10)* apply(batch_not_removed10$MAP$F,2,sd)),unlist(label_list)))
# 
# knn_reduced_core_shared_factor20<- sapply(c(1:20), function(k) KNNtest(k,t(matrix(unlist(batch_not_removed20$MAP$C_list),20)* apply(batch_not_removed10$MAP$F,2,sd)),unlist(label_list)))



# plot(knn_raw_combined_error,ylim=c(0,1),type="l")
# lines(knn_reduced_core_combined_error, col="red")
# lines(knn_reduced_core_shared_factor, col="blue")
# lines(knn_reduced_core_shared_factor20, col="magenta")




#batch removal total error
require("ggplot2")
all_label = unlist(label_list)
ns = sapply(label_list,length)
batchID = unlist(sapply(c(1:length(ns)),function(i){rep(i,ns[i])} ))

plotCore<- function(C_list,F,r){
  df = data.frame("value"=unlist(C_list)*apply(F,MARGIN = 2,sd),"index"=c(1:r),"subject"=rep(c(1:sum(ns)),each=r)
                  ,"batch"=rep(batchID,each=r),"label"=as.factor(rep(all_label,each=r)))
  ggplot(data=df, aes(x=index,y=value))+geom_line(aes(group=subject,col=label))+ facet_grid(~ batch)
}



pdf("5dataCoreRandomFactor.pdf",10,6)
plotCore(batch_removed$MAP$C_list, batch_removed$MAP$F,30)
dev.off()

pdf("5dataCoreSharedFactor.pdf",10,6)
plotCore(batch_not_removed10$MAP$C_list, batch_not_removed10$MAP$F,10)
dev.off()


pdf("5dataCoreSharedFactor20.pdf",10,6)
plotCore(batch_not_removed20$MAP$C_list, batch_not_removed20$MAP$F,20)
dev.off()



###


t(batch_removed$C_list)

plot(knn_raw_error[[1]])
plot(knn_raw_combined_error)

# avg
avgA<- lapply(A_list,function(x){
  A<- matrix(unlist(x),70*70)
  matrix(rowMeans(A),70)
  })


m_j<-5




adjustedFitted<- fitted_adjusted(batch_removed$C_list,batch_removed$F,5)


# avg
avgAadjusted<- lapply(adjustedFitted,function(x){
  A<- matrix(unlist(x),70*70)
  matrix(rowMeans(A),70)
})



pdf("avgA1.pdf",6,6)
image(avgA[[1]],zlim = c(0,1))
dev.off()
pdf("avgA2.pdf",6,6)
image(avgA[[2]],zlim = c(0,1))
dev.off()


pdf("avgA1adjusted.pdf",6,6)
image(avgAadjusted[[1]],zlim = c(0,1))
dev.off()
pdf("avgA2adjusted.pdf",6,6)
image(avgAadjusted[[2]],zlim = c(0,1))
dev.off()




p_random_factor<- fitted(C_list = batch_removed$MAP$C_list, F_list = batch_removed$MAP$F_list,5)
p_shared_factor<- fitted(C_list = batch_not_removed10$MAP$C_list, F_list = batch_not_removed10$MAP$F_list,1)




image(A_list[[2]][[2]],zlim = c(0,1))
image(p_batch_removed[[2]][[2]],zlim = c(0,1))

image(A_list[[3]][[3]],zlim = c(0,1))
image(p_batch_removed[[3]][[3]],zlim = c(0,1))

library(pROC)

auc(c(unlist(A_list[[1]][[20]])),c(unlist(p_random_factor[[1]][[20]])))
auc(c(unlist(A_list[[1]][[20]])),c(unlist(p_shared_factor[[1]][[20]])))



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

# save(testRun,file="resultDataA.RDa")

A_list_flat<- c(A_list[[1]],A_list[[2]],A_list[[3]])
A_list_new<- list()
A_list_new[[1]]<- A_list_flat
testRun2<- runBatchRemoval(A_list_new, r=r, 200)

# save(testRun2,file="resultDataANoBE.RDa")
load("resultDataANoBE.RDa")




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


##dicriminability


label<- unlist(label_list)

C_combined1<- t(matrix(unlist(testRun$C_list),r))
C_combined2<- t(matrix(unlist(testRun2$C_list),r))

require("class")

KNNtest<-function(k, C_combined){
  n<- length(label)
  cl<- sapply(1:length(label), function(i){
    true_cl = label[i]
    true_cl==knn(C_combined[-i,],test = C_combined[i,],cl=label[-i],k = k)
  })
  1-sum(cl)/n
}

KNNtestSubset<-function(k, C_combined, subset){
  label_subset<- label[subset]
  n<- length(label_subset)
  C_combined_subset <- C_combined[subset,]
  cl<- sapply(1:length(label_subset), function(i){
    true_cl = label[i]
    true_cl==knn(C_combined_subset[-i,],test = C_combined_subset[i,],cl=label_subset[-i],k = k)
  })
  1-sum(cl)/n
}


K<- 20




knn1<- sapply(c(1:K),function(k)KNNtest(k, C_combined1))
knn2<- sapply(c(1:K),function(k)KNNtest(k, C_combined2))
knn_raw_combined_error<- sapply(c(1:20), function(k) KNNtest(k,A_mat_lowertri))



##for method with fixed effects, remove group mean
C_combined3<- t(C_combined2)
C_combined3[,batchID==1]<-C_combined3[,batchID==1] - rowMeans(C_combined3[,batchID==1])
C_combined3[,batchID==2]<-C_combined3[,batchID==2] - rowMeans(C_combined3[,batchID==2])
C_combined3[,batchID==3]<-C_combined3[,batchID==3] - rowMeans(C_combined3[,batchID==3])
C_combined3<- t(C_combined3)

knn3<- sapply(c(1:K),function(k)KNNtest(k, C_combined3))


df<-data.frame("k"=c(1:K),"misclassification" = c(knn1,knn2,knn3),"method"=rep(c("random factor model","shared factor model","shared factor model (removing group mean)"),each=K))


pdf("misclassification_data_knn.pdf",8,5)
ggplot(data=df, aes(x=k,y=misclassification))+geom_line(aes(group=method,col=method))
dev.off()


##knn in each group
K<- 40
knn11<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined1, subset = batchID==1))
knn21<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined2, subset = batchID==1))
knn31<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined3, subset = batchID==1))

knn12<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined1, subset = batchID==2))
knn22<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined2, subset = batchID==2))
knn32<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined3, subset = batchID==2))

knn13<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined1, subset = batchID==3))
knn23<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined2, subset = batchID==3))
knn33<- sapply(c(1:K),function(k)KNNtestSubset(k, C_combined3, subset = batchID==3))



df1<-data.frame("k"=c(1:K),"misclassification" = c(knn11,knn21,knn31,knn12,knn22,knn32,knn13,knn23,knn33),"method"=rep(c("random factor model","shared factor model","shared factor model (removing group mean)"),each=K), "batch"=rep(c("BNU1","KKI2009","MRN114"),each= K*3))


pdf("misclassification_data_knn_each_batch.pdf",8,5)
ggplot(data=df1, aes(x=k,y=misclassification))+geom_line(aes(group=method,col=method))+ facet_grid(~ batch)
dev.off()
