setwd("~/git/semipar-cci/BatchEffectRemoval/")

load("./5dataset_list.RDa")

source("batchremoval.r")

# r=10
# batch_removed<- runBatchRemoval(A_list, r=r, 1,EM = T)

# save(batch_removed, file = "data_batch_removed.RDa")

A_list_flat<- c(A_list[[1]],A_list[[2]],A_list[[3]],A_list[[4]],A_list[[5]])
A_list_new<- list()
A_list_new[[1]]<- A_list_flat

# batch_not_removed10 <- runBatchRemoval(A_list_new, r=10, 1,EM = T)
# save(batch_not_removed10, file = "data_batch_not_removed_rank10.RDa")

batch_not_removed20 <- runBatchRemoval(A_list_new, r=20, 200,EM = T)
save(batch_not_removed20, file = "data_batch_not_removed_rank20.RDa")