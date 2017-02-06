#logit function
func_logit<- function(x){
  p<- 1/(1+exp(-x))
  p[p==1]<- 1-1E-6
  p[p==0]<- 1E-6
  p
}

#fitted functions
getFitted<- function(C_list,F_list){
  lapply(c(1:length(C_list)),function(j){
    lapply( c(1:ncol(C_list[[j]])), function(x){ 
      F_local = F_list[[j]]
      p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
      p
    } )
  })
}

#batch adjusted fitted functions
getFitted2<- function(C_list,F){
  lapply(c(1:length(C_list)),function(j){
    lapply( c(1:ncol(C_list[[j]])), function(x){ 
      F_local = F
      p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
      p
    } )
  })
}

#get flat lower triagular part
getLowTri<- function(Alist) sapply(Alist,function(x){    x[lower.tri(x)]  })

#plot heat map with ggplot2

require("reshape")
require("ggplot2")

ggplotHeatmap<- function(mat){
  df= melt(mat)
  ggplot(df, aes(X1, X2)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")
}


ggplotHeatmapList<- function(list_mat,mat_names){
  
  df=data.frame()
  for(i in 1:length(list_mat)){
    mat = list_mat[[i]]
    df= rbind(df, cbind(melt(mat), "name"=mat_names[i]))
  }
  
  ggplot(df, aes(X1, X2)) + geom_tile(aes(fill = value), colour = "white") + scale_fill_gradient(low = "white",high = "steelblue")+ facet_grid(~ name)
}



#batch removal total error
# 
# plotCore<- function(C_list,F,r){
#   df = data.frame("value"=unlist(C_list)*apply(F,MARGIN = 2,sd),"index"=c(1:r),"subject"=rep(c(1:sum(ns)),each=r)
#                   ,"batch"=rep(batchID,each=r),"label"=as.factor(rep(unlist(label_list),each=r)))
#   ggplot(data=df, aes(x=index,y=value))+geom_line(aes(group=subject,col=label))+ facet_grid(~ batch)
# }
# 

plotCore<- function(C_list){
  r= nrow(C_list[[1]])
  df = data.frame("value"=unlist(C_list),"index"=c(1:r),"subject"=rep(c(1:sum(ns)),each=r)
                  ,"batch"=rep(batchID,each=r),"label"=as.factor(rep(unlist(label_list),each=r)))
  ggplot(data=df, aes(x=index,y=value))+geom_line(aes(group=subject,col=label))+ facet_grid(~ batch)
}
