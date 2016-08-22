require("TensorEmbedding")
require("base")
setwd("~/git/semipar-cci/")

covariates<- read.csv(file = "Data/covariates/BNU3.csv")

covariateSubjectID<- unique(covariates$SUBID)

listGraphML<- list.files(path = "Data/desikan/BNU3", pattern = "*.graphml")

graphMLSubjectID<- as.numeric(unlist(lapply(listGraphML, function(x){strsplit(x,"_")[[1]][2]})))


intersectedID<- intersect( graphMLSubjectID, covariateSubjectID)


covariates2<- covariates[covariates$SUBID %in% intersectedID,]

listGraphML2<- listGraphML[graphMLSubjectID %in% intersectedID]

length(listGraphML2)

x<- listGraphML2[1]

BNU3<- lapply(listGraphML2, function(x){
  numID <- as.numeric(unlist(strsplit(x,"_"))[2])
  list(
    "A"= getA(paste("Data/desikan/BNU3/",x,sep="")),
    "id" = numID,
    "SEX" = min(covariates2$SEX[covariates2$SUBID==numID]),
    "AGE_AT_SCAN_1" = min(covariates2$AGE_AT_SCAN_1[covariates2$SUBID==numID])
  )
  })


pdf("Data/viz/BNU3.pdf",4,8)
par(mfrow=c(4,2))
for(it in BNU3){
  image(it$A)
  title(it$id)
}
dev.off()

save(BNU3,file="BNU3.RDa")


pdf("Data/viz/BNU3.pdf",4,8)
par(mfrow=c(4,2))
sumASex1=matrix(0,n,n)
sumASex2=matrix(0,n,n)
countSex1=0
countSex2=0
for(it in BNU3){
  if(it$SEX==1){
    sumASex1<- sumASex1+it$A
    countSex1<- countSex1+1
  }else{
    sumASex2<- sumASex2+it$A
    countSex2<- countSex2+1
  }
}

image(sumASex1/countSex1)
title("Sex=1")
image(sumASex2/countSex2)
title("Sex=2")

