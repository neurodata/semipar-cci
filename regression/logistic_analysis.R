setwd("~/git/ABadea/")

load("./listData.Rda")

load("./tensorDecompBinary.Rda")

L<- tensorDecomp$L
C<- tensorDecomp$C

n<- dim(L)[1]
p<- dim(C)[1]
k<- dim(C)[2]

diagC<- C

diagC

plot(c(1,k),range(c(C)),type = "n")
for(i in 1:p){
  lines(diagC[i,])
}

covariates<- read.csv("./covariates/predictors.csv",stringsAsFactors = F)

ids <- unlist( lapply(dataList,function(x)strtrim(x$id,6)))

common_id<- intersect(covariates$RUNNO , ids)

diagC_1<- diagC[ids%in% common_id,]
covariates_1<- covariates[covariates$RUNNO %in% common_id,]

diagC_2<- diagC_1[order(ids[ids%in% common_id]),]
covariates_2<- covariates_1[order(covariates_1$RUNNO),]


p<- nrow(diagC_2)
plot(c(1,k),range(c(C)),type = "n")
for(i in 1:p){
  lines(diagC_2[i,],col=covariates_2$GENDER[i])
}

plot(c(1,k),range(c(C)),type = "n")
for(i in 1:p){
  lines(diagC_2[i,],col=covariates_2$GENOTYPE[i])
}
