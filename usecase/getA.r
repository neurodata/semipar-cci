require("igraph")

getA<- function(x){
  graph <- read.graph(x, format='graphml')
  graph1 <- as.undirected(graph, mode = "collapse")
  adjacency <- (get.adjacency(graph1, attr='weight', sparse=FALSE))>1
  adjacency*1
}