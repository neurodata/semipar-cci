A_list_lowertri<- list()

A_list_lowertri<- lapply(A_list, function(x){
  t(sapply(x,function(x){x[lower.tri(x)]}))
})

A_list_flat<- lapply(A_list, function(x){
  t(sapply(x,function(x){c(x)}))
})

knn_loo_error_raw_data<- sapply(c(1:5), function(j) knnLOOVaryingK(M = A_list_lowertri[[j]], label = label_list[[j]]))

knn_loo_error_raw_data_flat<- sapply(c(1:4), function(j) knnLOOVaryingK(M = A_list_flat[[j]], label = label_list[[j]]))

