setwd("F:\\work/git/semipar-cci/BatchEffectRemoval/")
source("getFitted.r")
source("classification.r")

load("./RDa/5dataset_list.RDa")
load("./RDa/5dataset_sex.RDa")

# load("./RDa/data_batch_removed_1to4.RDa")
load("./RDa/batch_adjusted.RDa")
load("./RDa/batch_not_adjusted.RDa")
# load("./RDa/data_batch_not_removed_rank20.RDa")
A_lowtri<- lapply(A_list, function(y) sapply(y, function(x){x[lower.tri(x)]}))


A_list[[4]]<- NULL
A_list[[2]]<- NULL
label_list[[4]]<- NULL
label_list[[2]]<- NULL
A_lowtri[[4]]<- NULL
A_lowtri[[2]]<- NULL

# knn_raw_combined_error<- knnLOOVaryingK(M = t(do.call("cbind",A_lowtri)), label = unlist(label_list))
# save(knn_raw_combined_error, file="knn_raw_combined_error.RDa")

load("./old/knn_raw_error.RDa")
load("./knn_raw_combined_error.RDa")

knn_raw_error[[4]]<- NULL
knn_raw_error[[2]]<- NULL

m_j = length(A_list)

ns = sapply(label_list,length)
batchID = unlist(sapply(c(1:length(ns)),function(i){rep(i,ns[i])} ))

# pdf("core.pdf")
plotCore(batch_removed$C_list)
# dev.off()

plotCore(batch_not_removed$C_list)

# A_lowtri[[4]]<- NULL
# ts.plot(batch_removed$trace_loss[10:1000])
# ts.plot(batch_removed$trace_sigma2[10:1000])



#############
#get LOO knn error for random and shared factor model 
############

knn_loo_error_random_factor<- sapply(c(1:m_j), function(j) knnLOOVaryingK(M = t(batch_removed$C_list[[j]]), label = label_list[[j]]))

sf_C_list = list()
for(j in 1:m_j){
  sf_C_list[[j]] = batch_not_removed$C_list[[1]][,batchID==j] 
}


knn_loo_error_shared_factor<- sapply(c(1:m_j), function(j) knnLOOVaryingK(M = t(sf_C_list[[j]]), label = label_list[[j]]))

# knnLOOVaryingK(M = t(sf_C_list[[j]] ), label = label_list[[j]])

# knn_loo_error_random_factor1<- sapply(c(1:4), function(j) knnLOOVaryingK(M = t(batch_removed$C_list[[j]]), label = label_list[[j]]))



# df = melt(do.call("cbind",knn_raw_error))
# df = melt(cbind(do.call("cbind",knn_raw_error)[1:5,c(1,3,5)],c(knn_raw_combined_error)))
# 
# colnames(df)<- c("k","dataset","error")
# 
# df$dataset[df$dataset==4]<- "combined"
# df$dataset <- as.factor(df$dataset)
# 
# pdf("selection_dataset_good.pdf",8,6)
# ggplot(data=df, aes(x=k,y=error))+geom_line(aes(group=dataset,col=dataset))+theme_bw()
# dev.off()


knn_raw_error<- sapply(knn_raw_error,function(x)x[1:5])


pdf("knnErrorPerDataset.pdf",10,6)
df<-data.frame("k"=c(1:5),"error" = c(knn_raw_error,knn_loo_error_random_factor,knn_loo_error_shared_factor),"method"=rep(c("adjacency","random factor model","shared factor model"),each=5*3), "batch"=rep(c("1","2","3"),each= 5))
ggplot(data=df, aes(x=k,y=error))+geom_line(aes(group=method,col=method))+ facet_grid(~ batch)
dev.off()



par(mfrow=c(2,2))
for(j in 1:3){
  plot(knn_raw_error[[j]],type="l",ylim=c(0,1), main=j,xlim=c(1,6))
  lines(knn_loo_error_random_factor[,j],col="red")
  lines(knn_loo_error_shared_factor[,j],col="blue")
}
par(mfrow=c(1,1))

# min(knn_raw_error[[4]])
# min(knn_loo_error_random_factor[,4])

###########
#compute the combined LOO knn error
############
combined_knn_loo_error_random_factor<-  knnLOOVaryingK(M = t(do.call("cbind",batch_removed$C_list)), label = unlist(label_list))
combined_knn_loo_error_shared_factor<-  knnLOOVaryingK(M = t(do.call("cbind",batch_not_removed$C_list)), label = unlist(label_list))


pdf("knnErrorCombined.pdf",6,6)
df<-data.frame("k"=c(1:5),"error" = c(knn_raw_combined_error,combined_knn_loo_error_random_factor,combined_knn_loo_error_shared_factor),"method"=rep(c("adjacency","random factor model","shared factor model"),each=5))
ggplot(data=df, aes(x=k,y=error))+geom_line(aes(group=method,col=method))
dev.off()

########
#compute distance metric
####

p_random_factor_model <- getFitted(batch_removed$MAP$C_list,batch_removed$MAP$F_list, batch_removed$MAP$Z_list)
adjusted_p_random_factor_model <- getFitted2(batch_removed$MAP$C_list,batch_removed$MAP$F)


A_lowtri<- lapply(A_list,getLowTri)
p_lowtri<-  lapply(p_random_factor_model,getLowTri)
adjusted_p_lowtri<-  lapply(adjusted_p_random_factor_model,getLowTri)

diffMat_A<- matrix(0,3,3)
diffMat_p<- matrix(0,3,3)
diffMat_adjusted_p<- matrix(0,3,3)

for(i in 1:m_j){
  for(j in 1:m_j){
    diffMat_A[i,j]<- sqrt(mean((rowMeans(A_lowtri[[i]]) - rowMeans(A_lowtri[[j]]))^2))
    diffMat_p[i,j]<- sqrt(mean((rowMeans(p_lowtri[[i]]) - rowMeans(p_lowtri[[j]]))^2))
    diffMat_adjusted_p[i,j]<- sqrt(mean((rowMeans(adjusted_p_lowtri[[i]]) - rowMeans(adjusted_p_lowtri[[j]]))^2))
  }
}


# image(batch_removed$F_list[[1]])
# image(batch_removed$F_list[[2]])
# image(batch_removed$F_list[[3]])
# image(batch_removed$F_list[[4]])
# image(batch_removed$F)
 

image(A_list[[1]][[1]])
image(p_random_factor_model[[1]][[1]])

image(A_list[[3]][[4]])
image(p_random_factor_model[[3]][[4]])



A_flat<- lapply(A_list,function(x) rowMeans(sapply(x, c)))
p_flat<- lapply(p_random_factor_model,function(x) rowMeans(sapply(x, c)))
p_adjusted_flat<-  lapply(adjusted_p_random_factor_model,function(x) rowMeans(sapply(x, c)))


pdf("fitted_dataset.pdf",10,10)
par(mfrow=c(3,3))
for(j in 1:3){
image(matrix(A_flat[[j]],70),main="Adjacency")
image(matrix(p_flat[[j]],70),zlim=c(0,1),main="P")
image(matrix(p_adjusted_flat[[j]],70),zlim=c(0,1),main="P-Adjusted")
}
par(mfrow=c(1,1))
dev.off()


pdf("distance_between_dataset.pdf",10,4)
ggplotHeatmapList(list(diffMat_A,diffMat_p,diffMat_adjusted_p ), c("A","P","P-adjusted"))
dev.off()

barplot(diffMat_A[lower.tri(diffMat_A)])
barplot(diffMat_p[lower.tri(diffMat_A)],add = T,col="blue")
barplot(diffMat_adjusted_p[lower.tri(diffMat_A)],add=T,col="red")



