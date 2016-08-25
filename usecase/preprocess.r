require("TensorEmbedding")
require("base")
setwd("~/git/semipar-cci/")
source("usecase/getA.r")

dataset<- "BNU3"
dataset<- "BNU1"
dataset<- "KKI2009"
dataset<- "MRN114"
dataset<- "NKI1" #only has male
dataset<- "SWU4" #large


covariatePath<- paste("Data/covariates/",dataset,".csv",sep="")
covariates<- read.csv(file = covariatePath)

#BNU1
# covariates$SEX


#KKI2009
# covariates$SEX<- as.numeric(covariates$Sex)
# covariates$SUBID<- covariates$SubjectID


#MRN114
# covariates$SEX<- as.numeric(covariates$Sex)+1
# covariates$SUBID<- covariates$URSI

#NKI1
# covariates$SEX
#SWU4

covariates <- covariates[covariates$SEX%in%c(1,2),]

covariates$SEX

graphmlPath<- paste("Data/desikan/",dataset,sep="")

listGraphML<- list.files(path = graphmlPath, pattern = "*.graphml")

graphMLSubjectID<- (unlist(lapply(listGraphML, function(x){strsplit(x,"_")[[1]][2]})))

graphMLSubjectID



dataList<- lapply(listGraphML, function(x){
  numID <- as.numeric(unlist(strsplit(x,"_"))[2])
  if(numID%in%covariates$SUBID){
    list(
      "A"= getA(paste("Data/desikan/",dataset,"/",x,sep="")),
      "id" = numID,
      "SEX" = min(as.numeric(as.character(covariates$SEX))[covariates$SUBID==numID]),
      "AGE_AT_SCAN_1" = min(as.numeric(as.character(covariates$AGE_AT_SCAN_1))[covariates$SUBID==numID])
    )
  }else{NULL}
})

dataList<- dataList[!sapply(dataList, is.null)]

lapply(dataList, function(x)print(x$id))

pdf(paste("Data/viz/",dataset,".pdf",sep=""),4,8)
par(mfrow=c(4,2))
for(i in 1:length(dataList)){
  it<- dataList[[i]]
  image(it$A)
  title(i)
}
dev.off()

save(dataList,file=paste("Data/processed/",dataset,".RDa",sep=""))

