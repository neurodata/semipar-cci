setwd("~/git/semipar-cci/BatchEffectRemoval/")

n<- 70
r<- 5
F<- matrix(rnorm(n*r),n)

F0<- F
sigma2<- 0.5^2

F_list<- list()
F_list[[1]]<-  F + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)
F_list[[2]]<-  F + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)
F_list[[3]]<-  F + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)

m_j = 3

image(F_list[[1]])
image(F_list[[2]])
image(F_list[[3]])
image(F)

n1<- 50
n2<- 50
n3<- 50

C_list<- list()
C<- matrix( c(rep( 0.1* c(1:r), n1/2), rep( 1/ c(1:r)^3 ,n1/2)), nrow = r)

C_list[[1]] <- C + matrix( rnorm(r*n1, sd =0.01), r)
C_list[[2]] <- C + matrix( rnorm(r*n2, sd =0.01), r)
C_list[[3]] <- C + matrix( rnorm(r*n3, sd =0.01), r)
C_list0 <- C_list
# image(C_list[[1]])
# image(C_list[[2]])
# image(C_list[[3]])

image(rbind(C_list0[[1]],C_list0[[2]],C_list0[[3]]))


func_logit<-  function(x){1/(1+exp(-x))}
func_ilogit<-  function(x){log(x/(1-log(x)))}

p_list<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list[[j]])), function(x){ 
      
      F_local = F_list[[j]]
      p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
      p
    } )
  })

A_list<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list[[j]])), function(x){ 
    
    p<- p_list[[j]][[x]]
    A<- (p> runif(length(p)))*1
    A[lower.tri(A)]<- t(A)[lower.tri(A)]
    A
  } )
})


rm(C_list)
rm(F_list)

# image(A_list[[1]][[1]])
# image(A_list[[1]][[2]])
# image(A_list[[1]][[3]])
# image(A_list[[1]][[4]])

source("batchremoval.r")

testRun<- runBatchRemoval(A_list, r=5, 100,q=2)

# image(func_logit(testRun$Z_list[[1]]),zlim = c(0,1))

Z_list =  testRun$Z_list

F_list =  testRun$F_list
C_list =  testRun$C_list
F = testRun$F
sigma2 = testRun$sigma2


sqrt(sigma2)/sd(F)


ts.plot(do.call("rbind", C_list))



p_est_list<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list[[j]])), function(x){ 
    F_local = F_list[[j]]
    p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]) + Z_list[[j]])
    p
  } )
})


plot(c(p_list[[1]][[1]]), c(p_est_list[[1]][[1]]), ylab="generating_prob",xlab="fitted")
plot(c(p_list[[2]][[1]]), c(p_est_list[[2]][[1]]), ylab="generating_prob",xlab="fitted")


# plot(c(p_list[[1]][[1]]), c(func_logit(testRun$Z_list[[1]])), ylab="generating_prob",xlab="fitted")
# plot(c(p_list[[3]][[1]]), c(func_logit(testRun$Z_list[[3]])), ylab="generating_prob",xlab="fitted")







p_list_wo_batch<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list[[j]])), function(x){ 
    F_local = F
    p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]) + Z_list[[j]])
    p
  } )
})


p_truth_wo_batch<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list0[[j]])), function(x){ 
    F_local = F0
    p<- func_logit(F_local%*% (t(F_local)*C_list0[[j]][,x])+ Z_list[[j]])
    p
  } )
})

plot(c(p_truth_wo_batch[[1]][[1]]), c(p_list_wo_batch[[1]][[1]]), ylab="truth",xlab="est_wo_be")

plot(c(p_list[[1]][[1]]), c(p_list_wo_batch[[1]][[1]]), ylab="generating_prob",xlab="est_wo_be")
plot(c(p_list_wo_batch[[1]][[2]]),c(p_est_list[[1]][[2]]) - c(p_list_wo_batch[[1]][[2]]), xlab="est_wo_batch",ylab="diff")


# image(A_list[[1]][[1]],zlim = c(0,1))
image(p_est_list[[1]][[1]],zlim = c(0,1))
image(p_list[[1]][[1]],zlim = c(0,1))



image(p_list[[1]][[2]],zlim = c(0,1))
image(p_list_wo_batch[[1]][[2]],zlim = c(0,1))

image(p_est_list[[1]][[1]],zlim = c(0,1))
image(p_est_list[[2]][[1]],zlim = c(0,1))
image(p_est_list[[3]][[1]],zlim = c(0,1))

B1<- (F_list[[1]] -F)
FCB = F%*%diag(C_list[[1]][,1])%*%t(B1)
BCB = B1%*%diag(C_list[[1]][,1])%*%t(B1)

overallB = FCB+t(FCB)+BCB
hist(overallB)
image(overallB)

plot(C_list[[1]][,2])
for(i in 1:20){
  lines(C_list[[1]][,i])
}


require("mclust")
M<- rep(1,50)
M[1:25]<-0



plot(C_list[[1]][,1],type="n", ylim=c(-3,3))
for(j in 1:3){
  for(i in 1:50){
    lines(C_list[[j]][,i]*apply(F,MARGIN = 2,sd),col= (i>25)+1, lty=j)
  }
}



M1<- kmeans(t(C_list[[1]]*apply(F,MARGIN = 2,sd)),2)$cluster
M2<- kmeans(t(C_list[[2]]*apply(F,MARGIN = 2,sd)),2)$cluster
M3<- kmeans(t(C_list[[3]]*apply(F,MARGIN = 2,sd)),2)$cluster

adjustedRandIndex(M1,M)
adjustedRandIndex(M2,M)
adjustedRandIndex(M3,M)

require("ggplot2")

df = data.frame("value"=unlist(C_list)*apply(F,MARGIN = 2,sd),"index"=c(1:5),"subject"=rep(c(1:n1),each=5),"batch"=rep(c(1:3),each=5*n1),"label"=as.factor(rep(c(1,2),each=(5*n1/2))))

pdf("batch_removal.pdf",10,6)
ggplot(data=df, aes(x=index,y=value))+geom_line(aes(group=subject,col=label))+ facet_grid(~ batch)
dev.off()

save(A_list,file="simA.RDa")


A_list_flat<- c(A_list[[1]],A_list[[2]],A_list[[3]])
A_list_new<- list()
A_list_new[[1]]<- A_list_flat
testRun2<- runBatchRemoval(A_list_new, r=5, 200)

(testRun2$trace_sigma2)



C_list =  testRun2$MAP$C_list
F =  testRun2$MAP$F

df = data.frame("value"=unlist(C_list)*apply(F,MARGIN = 2,sd),"index"=c(1:5),"subject"=rep(c(1:n1),each=5),"batch"=rep(c(1:3),each=5*n1),"label"=as.factor(rep(c(1,2),each=(5*n1/2))))

pdf("BE_not_removed.pdf",10,6)
ggplot(data=df, aes(x=index,y=value))+geom_line(aes(group=subject,col=label))+ facet_grid(~ batch)
dev.off()
