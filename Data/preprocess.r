

#BNU1
setwd("~/git/semipar-cci/Data/")
source("getA.r")
require("base")
dataset_list<- c("BNU1", "BNU3" , "KKI2009", "MRN114","SWU4")



dataset= dataset_list[[1]]
covariatePath<- paste("covariates/",dataset,".csv",sep="")
covariates<- read.csv(file = covariatePath)
covariates= covariates[covariates$SESSION== "Baseline",]
covariates$SEX = as.numeric(as.character(covariates$SEX))

graphmlPath<- paste("desikan/",dataset,sep="")
listGraphML<- list.files(path = graphmlPath, pattern = "*.graphml")
graphMLSubjectID<- (unlist(lapply(listGraphML, function(x){strsplit(x,"_")[[1]][2]})))

common_id = intersect( covariates$SUBID, as.numeric(graphMLSubjectID))

A_list =list()
Sex_list=numeric()
id_list=numeric()
for(i in 1:length(common_id)){
  x= common_id[[i]]
  A_list[[i]] = getA(paste("desikan/BNU1/BNU1_00",x,"_1_DTI_desikan.graphml",sep=""))
  Sex_list =c (Sex_list, as.numeric(covariates$SEX)[covariates$SUBID==x])
  id_list= c(id_list, x)
}

pdf(paste("viz/",dataset,".pdf",sep=""),4,8)
par(mfrow=c(4,2))
for(A in A_list){
  image(A)
}
dev.off()

Data1 = list("A"=A_list,"SEX"=Sex_list,"id"=id_list)
save(Data1, file="processed/BNU1.RDa")


#BNU3
rm(list=ls())
dataset= dataset_list[[2]]
covariatePath<- paste("covariates/",dataset,".csv",sep="")
covariates<- read.csv(file = covariatePath)
covariates= covariates[covariates$SESSION== "Baseline",]
covariates$SEX = as.numeric(as.character(covariates$SEX))

graphmlPath<- paste("desikan/",dataset,sep="")
listGraphML<- list.files(path = graphmlPath, pattern = "*.graphml")
graphMLSubjectID<- (unlist(lapply(listGraphML, function(x){strsplit(x,"_")[[1]][2]})))

common_id = intersect( covariates$SUBID, as.numeric(graphMLSubjectID))

A_list =list()
Sex_list=numeric()
id_list=numeric()
for(i in 1:length(common_id)){
  x= common_id[[i]]
  A_list[[i]] = getA(paste("desikan/BNU3/BNU3_00",x,"_1_DTI_desikan.graphml",sep=""))
  Sex_list =c (Sex_list, as.numeric(covariates$SEX)[covariates$SUBID==x])
  id_list= c(id_list, x)
}

pdf(paste("viz/",dataset,".pdf",sep=""),4,8)
par(mfrow=c(4,2))
for(A in A_list){
  image(A)
}
dev.off()

Data1 = list("A"=A_list,"SEX"=Sex_list,"id"=id_list)
save(Data1, file="processed/BNU3.RDa")



#KKI2009
rm(list=ls())

setwd("~/git/semipar-cci/Data/")
source("getA.r")
require("base")
dataset_list<- c("BNU1", "BNU3" , "KKI2009", "MRN114","SWU4")


dataset= dataset_list[[3]]
covariatePath<- paste("covariates/",dataset,".csv",sep="")
covariates<- read.csv(file = covariatePath)

covariates = cbind(covariates, "SESSION"=rep(c(1:2)))
covariates= covariates[covariates$SESSION== 1,]
covariates$SEX = (covariates$Sex=="F")*1+(covariates$Sex=="M")*2
covariates$SUBID = covariates$SubjectID

graphmlPath<- paste("desikan/",dataset,sep="")
listGraphML<- list.files(path = graphmlPath, pattern = "*.graphml")
graphMLSubjectID<- (unlist(lapply(listGraphML, function(x){strsplit(x,"_")[[1]][2]})))

common_id = intersect( covariates$SUBID, as.numeric(graphMLSubjectID))

A_list =list()
Sex_list=numeric()
id_list=numeric()
for(i in 1:length(common_id)){
  x= common_id[[i]]
  A_list[[i]] = getA(paste("desikan/KKI2009/KKI2009_",x,"_1_DTI_desikan.graphml",sep=""))
  Sex_list =c (Sex_list, as.numeric(covariates$SEX)[covariates$SUBID==x])
  id_list= c(id_list, x)
}

pdf(paste("viz/",dataset,".pdf",sep=""),4,8)
par(mfrow=c(4,2))
for(A in A_list){
  image(A)
}
dev.off()

Data1 = list("A"=A_list,"SEX"=Sex_list,"id"=id_list)
save(Data1, file="processed/KKI2009.RDa")



#MRN114
rm(list=ls())

setwd("~/git/semipar-cci/Data/")
source("getA.r")
require("base")
dataset_list<- c("BNU1", "BNU3" , "KKI2009", "MRN114","SWU4")


dataset= dataset_list[[4]]
covariatePath<- paste("covariates/",dataset,".csv",sep="")
covariates<- read.csv(file = covariatePath)

covariates$SEX = covariates$Sex+1
covariates$SUBID = as.character(covariates$URSI)

graphmlPath<- paste("desikan/",dataset,sep="")
listGraphML<- list.files(path = graphmlPath, pattern = "*.graphml")
graphMLSubjectID<- (unlist(lapply(listGraphML, function(x){strsplit(x,"_")[[1]][2]})))

common_id = intersect( covariates$SUBID, graphMLSubjectID)

A_list =list()
Sex_list=numeric()
id_list=numeric()
for(i in 1:length(common_id)){
  x= common_id[[i]]
  print(x)
  A_list[[i]] = getA(paste("desikan/MRN114/MRN114_",x,"_1_DTI_desikan.graphml",sep=""))
  Sex_list =c (Sex_list, as.numeric(covariates$SEX)[covariates$SUBID==x])
  id_list= c(id_list, x)
}

pdf(paste("viz/",dataset,".pdf",sep=""),4,8)
par(mfrow=c(4,2))
for(A in A_list){
  image(A)
}
dev.off()

Data1 = list("A"=A_list,"SEX"=Sex_list,"id"=id_list)
save(Data1, file="processed/MRN114.RDa")




#SWU4
rm(list=ls())

setwd("~/git/semipar-cci/Data/")
source("getA.r")
require("base")
dataset_list<- c("BNU1", "BNU3" , "KKI2009", "MRN114","SWU4")


dataset= dataset_list[[5]]
covariatePath<- paste("covariates/",dataset,".csv",sep="")
covariates<- read.csv(file = covariatePath)

covariates= covariates[covariates$SESSION== "Baseline",]


covariates$SEX = as.numeric(as.character(covariates$SEX))
covariates$SUBID

graphmlPath<- paste("desikan/",dataset,sep="")
listGraphML<- list.files(path = graphmlPath, pattern = "*.graphml")
graphMLSubjectID<- as.numeric(unlist(lapply(listGraphML, function(x){strsplit(x,"_")[[1]][2]})))

common_id = intersect( covariates$SUBID, graphMLSubjectID)

A_list =list()
Sex_list=numeric()
id_list=numeric()
for(i in 1:length(common_id)){
  x= common_id[[i]]
  print(x)
  A_list[[i]] = getA(paste("desikan/SWU4/SWU4_00",x,"_1_DTI_desikan.graphml",sep=""))
  Sex_list =c (Sex_list, as.numeric(covariates$SEX)[covariates$SUBID==x])
  id_list= c(id_list, x)
}

pdf(paste("viz/",dataset,".pdf",sep=""),4,8)
par(mfrow=c(4,2))
for(A in A_list){
  image(A)
}
dev.off()

Data1 = list("A"=A_list,"SEX"=Sex_list,"id"=id_list)
save(Data1, file="processed/SWU4.RDa")
