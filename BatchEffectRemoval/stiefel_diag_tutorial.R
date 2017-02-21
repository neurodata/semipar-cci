setwd("c://work/git/semipar-cci/BatchEffectRemoval/")

#load the data that contains a list of adjacency matrix
loadData <- function(path) {
  load(path)
  Data1
}
DataList = loadData("../Data/processed/BNU1.RDa")
A_list = DataList$A

source("stiefel_diagonalize.r")

#chooes rank 30 and run stiefel joint embedding
#stiefelDecomp$U0 is the shared factor
#stiefelDecomp$D_list is the list of loadings
r=5
stiefelDecomp = stiefel_diagonalize(A_list, r = r, diagD = T)

#we can then use this function to get a list of estimated p
estP = getPEst(stiefelDecomp$U0,  stiefelDecomp$D_list)

#here is comparison of raw adjacency matrix vs estimated p
image(A_list[[1]])
image(estP[[1]])
