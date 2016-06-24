require("TensorEmbedding")
require("igraph")


setwd("~/git/semipar-cci/TestData/openconnecto.me/mrdata/share/dti/ndmg_v0011/KKI2009/desikan/")


list_graphml<- list.files(path = ".", pattern="*.graphml")

getA<- function(x){
  graph <- read.graph(x, format='graphml')
  graph1 <- as.undirected(graph, mode = "collapse")
  adjacency <- (get.adjacency(graph1, attr='weight', sparse=FALSE))>1
  adjacency*1
}

listA<- lapply(list_graphml, getA)


n<- nrow(listA[[1]])
p<- length(listA)


A<- array(0,dim = c(n,n,p))
for(i in 1:p){
  A[,,i]<- listA[[i]]
}

k<-10

testObj<- TensorEmbedding::symmetric_tensor_decomp(A,n,p, k, 500, 1E-3, 1E-3)

L <- testObj$L
C <- testObj$C

# image(A[,,1], zlim=c(-3,3))
# image(A[,,2], zlim=c(-3,3))
# 
# image(A_mean1, zlim=c(-3,3))
# image(L%*% C[,,1]%*%t(L), zlim=c(-3,3))
# 
# 
# hist(A_mean1 - L%*% C[,,1]%*%t(L))
# hist(A_mean2 - L%*% C[,,6]%*%t(L))
# 
# image(C[,,1])
# image(C[,,6])

require(ggplot2)
require(reshape)

cMelted<- numeric()

for(i in 1:p){
  cGroup1 <- melt(C[, , i])
  cMelted<- rbind( cMelted, cbind(cGroup1, "group"=i))
}




out<- seq(range(C)[1], range(C)[2], length.out = k)

plot(c(1:k), out, type = "n" )#, ylim = c(0,50))
for(i in 1:(p)){
  lines(c(1:k), diag(C[,,i]), col = floor((i-1)/2))
}




hist(C[1,1,], breaks = 10)
hist(C[2,2,], breaks = 10)
hist(C[3,3,], breaks = 10)
hist(C[4,4,], breaks = 10)
hist(C[5,5,], breaks = 10)
hist(C[6,6,], breaks = 10)
hist(C[7,7,], breaks = 10)


