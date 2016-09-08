setwd("~/git/ABadea/")

if(FALSE){
source("~/git/semipar-cci/usecase/getAweighted.r")
listGraphML<- list.files(path = "./graphml/", pattern = "*.graphml")


dataList<- lapply(listGraphML, function(x){
  numID <- unlist(x)
    list(
      "A"= getAweighted(paste("./graphml/",x,sep="")),
      "id" = numID
    )
})

length(dataList)


save(dataList,file="./listData.Rda")
}

load("./listData.Rda")

tensorACount<- do.call("c", lapply( dataList, function(x) length(x$A)))
ids<- do.call("c", lapply( dataList, function(x) x$id))

# ids[tensorACount==110889]
# table(tensorACount)

tensorA<- do.call("cbind", lapply( dataList, function(x) c(x$A)))

n<- sqrt( nrow(tensorA))
p<- ncol(tensorA)
k<- 10

tensorDecomp<- TensorEmbedding::symmetric_tensor_decompEM(tensorA,n,p, k, steps = 10000, 1E-6, 1E-6,loss_type = 2, tol = 1E-12,restrictCoreToDiag = TRUE)

save(tensorDecomp,file="tensorDecomp.Rda")
