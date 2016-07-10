require("TensorEmbedding")
require("igraph")
require("digest")

logit<- function(x){ 1/(1+exp(-x))}

setwd("~/git/semipar-cci/TestData/openconnecto.me/mrdata/share/dti/ndmg_v0011/")

list_batch<- list.files(path = ".")

list_batch<- list_batch[ list_batch  %in% c("BNU1",    "BNU3")]


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

save(tensorDecomp, file = "tensorDecompSubj.Rda")

L <- tensorDecomp$L
C <- tensorDecomp$C
diagC<- tensorDecomp$diagC


mG<- nrow(diagC)

batch_serial_no<- as.numeric(unclass(factor(batch_id)))

batch_col<- tapply(batch_serial_no, subject_serial_no, function(x){x[1]})

pdf("2batch.pdf",8,6)
plot(c(1:k), seq(0,30, length.out = k), type = "n" )#, ylim = c(0,50))
for(i in 1:mG){
  lines(c(1:k), diagC[i,], col = batch_col[i])
}
dev.off()

pdf("fitted1.pdf",8,4)
par(mfrow=c(1,2))
image(tensorA[,,1],zlim = c(0,1))
image(logit(L%*%C[,,1]%*%t(L)),zlim = c(0,1))
dev.off()


require("pROC")

roc(c(tensorA[,,1:10]), rep(c(logit(L%*%C[,,1]%*%t(L))),10), plot = T)


