require("TensorEmbedding")
require("igraph")
require("digest")

logit<- function(x){ 1/(1+exp(-x))}

setwd("~/git/semipar-cci/TestData/openconnecto.me/mrdata/share/dti/ndmg_v0011/")

list_batch<- list.files(path = ".")

list_batch<- list_batch[ list_batch  %in% c("BNU3",    "NKIENH")]


getAforGroup <- function(groupPath) {
  path = paste(groupPath, "/desikan/" , sep = "")
  
  list_graphml <- list.files(path = path, pattern = "*.graphml")
  subject_id <-
    do.call("c", lapply(strsplit(list_graphml, "_"), function(x) {
      x[2]
    }))
  n_matrices <- length(subject_id)
  batch_id <- rep(groupPath, n_matrices)
  
  getA <- function(x) {
    print(x)
    x_path <- paste(path, x, sep = "")
    graph <- read.graph(x_path, format = 'graphml')
    graph1 <- as.undirected(graph, mode = "collapse")
    adjacency <- numeric()
    try(adjacency <-
          (get.adjacency(graph1, attr = 'weight', sparse = FALSE)) > 1)
    return(adjacency * 1)
  }
  
  listA <- lapply(list_graphml, getA)
  
  A <- numeric()
  
  for (i in 1:length(listA)) {
    mat <- listA[[i]]
    if (length(mat) > 0)
      A <- c(A, listA[[i]])
    else{
      subject_id[i] <- NA
      batch_id[i] <- NA
    }
  }
  
  subject_id <- subject_id[!is.na(subject_id)]
  batch_id <- batch_id[!is.na(batch_id)]
  n_matrices <- length(subject_id)
  
  list(
    "A" = A,
    "batch_id" = batch_id,
    "subject_id" = subject_id,
    "n_matrices" = n_matrices
  )
}


listGroupA<- lapply(list_batch,getAforGroup)


A<- do.call("c", lapply(listGroupA, function(x){x$A}))
m<- do.call("sum", lapply(listGroupA, function(x){x$n_matrices}))

subject_id <- do.call("c", lapply(listGroupA, function(x){x$subject_id}))
batch_id <- do.call("c", lapply(listGroupA, function(x){x$batch_id}))

n<- sqrt(length(A)/m)

tensorA <-  array(A, dim = c(n, n, m))

batch_subject_id<- apply(cbind(batch_id, subject_id), MARGIN = 1,function(x){paste(x, collapse="_")})

subject_serial_no <- as.numeric(unclass(factor(batch_subject_id)))


k<- 30

tensorDecomp<- TensorEmbedding::symm_group_tensor_decomp(tensorA, subject_serial_no, n, m, k, 500, 1E-5, 1E-5)

save(tensorDecomp, file = "~/git/semipar-cci/Tests/tensorDecompSubj1.Rda")

# load(file= "tensorDecompSubj.Rda")

L <- tensorDecomp$L
C <- tensorDecomp$C
diagC<- tensorDecomp$diagC

#get batch average:
batch_id2<- tapply(batch_id, batch_subject_id, function(x){x[1]})
batch_subject_id2<- tapply(batch_subject_id, batch_subject_id, function(x){x[1]})


batchMeanC<- numeric()
for(batch in unique(batch_id2)){
  pick <- batch_id2 == batch
  batchMeanC<- rbind(batchMeanC, colMeans(diagC[pick,]))
}

plot(batchMeanC[1,])
lines(batchMeanC[2,])



popMean<- colMeans(batchMeanC)
batchOffset<-t( popMean - t(batchMeanC))

k<-1
correctedDiagC<- numeric()
for(batch in unique(batch_id2)){
  pick <- batch_id2 == batch
  correctedDiagC<- rbind(correctedDiagC, t(t(diagC[pick,])+batchOffset[k,]))
  k<- k+1
}


select<- c(1,2,10,60,100)
setwd("~/git/semipar-cci/Plots/")
pdf("subjects.pdf",8,3)
par(mfrow=c(1,3))
for(i in select)
{
LCL0<- L%*%diag(diagC[i,])%*%t(L)
LCL<- L%*%diag(correctedDiagC[i,])%*%t(L)
pick<- subject_serial_no==i
dataA<- tensorA[,,pick]
if(length(dim(dataA))==3)
  dataA<- dataA[,,1]
image(dataA)
title("Raw Data")
image(logit(LCL0))
title("Fitted")
image(logit(LCL))
title("Batch Intercept Removed")
}
dev.off()



pdf("batch.pdf",8,3)
plot(c(1,30),c(min(diagC),max(diagC)),type="n")
for(i in 1:nrow(diagC)){
  lines(diagC[i,],col=as.factor(batch_id2)[i])
}
dev.off()
