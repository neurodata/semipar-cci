require("class")
knnLOO <-function(k, M,label){
  n<- length(label)
  cl<- sapply(1:length(label), function(i){
    true_cl = label[i]
    true_cl==knn(train = M[-i,],test = M[i,],cl=factor(label[-i]),k = k)
  })
  1-sum(cl)/n
}

knnLOOVaryingK <-function(M,label){
  sapply(c(1:5), function(k)knnLOO(k,M,label))
}


