require("class")
knnLOO <-function(k, M,label){
  n<- length(label)
  cl<- sapply(1:length(label), function(i){
    true_cl = label[i]
    true_cl==knn(train = M[-i,],test = M[i,],cl=factor(label[-i]),k = k)
  })
  1-sum(cl)/n
}

knnLOOVaryingK <-function(M,label,max_k=10){
  sapply(c(1:max_k), function(k)knnLOO(k,M,label))
}


