require("TensorEmbedding")
require("igraph")

logit<- function(x){ 1/(1+exp(-x))}

setwd("~/git/semipar-cci/TestData/openconnecto.me/mrdata/share/dti/ndmg_v0011/")

list_batch<- list.files(path = ".")

list_batch<- list_batch[ !list_batch  %in% c("HCP500","MRN1313")]

getAforGroup <- function(groupPath) {
  path = paste(groupPath, "/desikan/" , sep = "")
  
  list_graphml <- list.files(path = path, pattern = "*.graphml")
  
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
  
  n <- nrow((listA[[1]]))
  A <- numeric()
  
  for (i in listA) {
    if (length(i) > 0)
      A <- c(A, i)
  }
  
  A <- array(A, dim = c(n, n, length(A) / n / n))
  
  A
}


listGroupA<- lapply(list_batch,getAforGroup)

format(object.size(listGroupA), units="Mb")
save(listGroupA,file="./listA.Rda")


group <- numeric()

idx<- 0

for(As in listGroupA){
  count <- dim(As)[3]
  group <- c(group, rep(idx, count))
  idx <- idx + 1
}

n<- dim(listGroupA[[1]])[1]
m<- length(group)
A <- array( do.call("c", listGroupA), dim=c(n,n, m))

k<- 10

tensorDecomp<- TensorEmbedding::symm_group_tensor_decomp(A, group,n, m, k, 500, 1E-3, 1E-3)

L <- tensorDecomp$L
C <- tensorDecomp$C


out<- seq(range(C)[1], range(C)[2], length.out = k)

mG<- length(unique(group))

plot(c(1:k), out, type = "n" )#, ylim = c(0,50))
for(i in 1:mG){
  lines(c(1:k), diag(C[,,i]), col = floor((i-1)/2))
}

image(A[,,1],zlim = c(0,1))
image(A[,,105],zlim = c(0,1))
image(logit(L%*%C[,,2]%*%t(L)),zlim = c(0,1))


require("pROC")

roc(c(A[,,1:104]), rep(c(logit(L%*%C[,,1]%*%t(L))),104), plot = T)
roc(c(A[,,105:150]), rep(c(logit(L%*%C[,,2]%*%t(L))),46), plot = T)
roc(c(A[,,151:443]), rep(c(logit(L%*%C[,,3]%*%t(L))),293), plot = T)
roc(c(A[,,151:443]), rep(c(logit(L%*%C[,,2]%*%t(L))),293), plot = T)

hist(C[1,1,], breaks = 10)
hist(C[2,2,], breaks = 10)
hist(C[3,3,], breaks = 10)
hist(C[4,4,], breaks = 10)
hist(C[5,5,], breaks = 10)
hist(C[6,6,], breaks = 10)
hist(C[7,7,], breaks = 10)


