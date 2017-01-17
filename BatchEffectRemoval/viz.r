setwd("~/git/semipar-cci/batchRemoval/")

loadNetworkData<- function(path){
  load(path)
 return(dataList) 
}

BNU1<- loadNetworkData("../Data/processed/BNU1.RDa")
BNU3<- loadNetworkData("../Data/processed/BNU3.RDa")


image(BNU1[[1]]$A)
image(BNU3[[1]]$A)
