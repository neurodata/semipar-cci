setwd("~/git/semipar-cci/BatchEffectRemoval/")

n<- 70
r<- 5
F<- matrix(rnorm(n*r),n)

# i<- 3

# C<- rnorm(r)
# F2<- F%*% diag(C) %*% t(F)
# F[i,] %*%diag(C)%*%t(F[-i,])
# F2[i,-i]
# t(t(F[-i,]) * F[i,]) %*% C
# F[-(1:i),] %*%diag( F[i,]) 
# k<- 1
# X<- F[-k,]
# w<- abs(rnorm(9))


F0<- F





# Binv<- diag(1E-5,r)
# K<- rnorm(9)
# b<- rnorm(r)


sigma2<- 0.5^2

# samplePsi(X, w, K, Binv,b)
# Binv<- diag(1/sigma2, r)


F_list<- list()
F_list[[1]]<-  F + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)
F_list[[2]]<-  F + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)
F_list[[3]]<-  F + matrix(rnorm(n*r,sd=sqrt(sigma2)),n)

m_j = 3

image(F_list[[1]])
image(F_list[[2]])
image(F_list[[3]])
image(F)

n1<- 20
n2<- 30
n3<- 20

C_list<- list()
C<- 1/c(1:r)^3
C_list[[1]] <- C + matrix( rnorm(r*n1, sd =0.1), r)
C_list[[2]] <- C + matrix( rnorm(r*n2, sd =0.1), r)
C_list[[3]] <- C + matrix( rnorm(r*n3, sd =0.1), r)
C_list0 <- C_list

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

testRun<- runBatchRemoval(A_list, r=r, 200)

F_list =  testRun$F_list
C_list =  testRun$C_list
F = testRun$F
sigma2 = testRun$MAP$sigma2
trace_sigma2 = testRun$trace_sigma2

ts.plot(trace_sigma2)

sqrt(sigma2)/sd(F)


image(F)
image(F_list[[1]])
image(F_list[[2]])
image(F_list[[3]])

p_est_list<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list[[j]])), function(x){ 
    F_local = F_list[[j]]
    p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
    p
  } )
})


p_list_wo_batch<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list[[j]])), function(x){ 
    F_local = F
    p<- func_logit(F_local%*% (t(F_local)*C_list[[j]][,x]))
    p
  } )
})


p_truth_wo_batch<- lapply(c(1:m_j),function(j){
  lapply( c(1:ncol(C_list0[[j]])), function(x){ 
    F_local = F0
    p<- func_logit(F_local%*% (t(F_local)*C_list0[[j]][,x]))
    p
  } )
})



plot(c(p_list[[1]][[1]]), c(p_est_list[[1]][[1]]), ylab="generating_prob",xlab="fitted")
plot(c(p_truth_wo_batch[[1]][[1]]), c(p_list_wo_batch[[1]][[1]]), ylab="truth",xlab="est_wo_be")


plot(c(p_list[[1]][[1]]), c(p_list_wo_batch[[1]][[1]]), ylab="generating_prob",xlab="est_wo_be")
# plot(c(p_list_wo_batch[[1]][[2]]),c(p_est_list[[1]][[2]]) - c(p_list_wo_batch[[1]][[2]]), xlab="est_wo_batch",ylab="diff")


image(A_list[[1]][[2]],zlim = c(0,1))
image(p_est_list[[1]][[2]],zlim = c(0,1))
image(p_list_wo_batch[[1]][[2]],zlim = c(0,1))


image(p_list_wo_batch[[1]][[1]],zlim = c(0,1))
image(p_list_wo_batch[[2]][[1]],zlim = c(0,1))
image(p_list_wo_batch[[3]][[1]],zlim = c(0,1))

image(p_est_list[[1]][[1]],zlim = c(0,1))
image(p_est_list[[2]][[1]],zlim = c(0,1))
image(p_est_list[[3]][[1]],zlim = c(0,1))

B1<- (F_list[[1]] -F)
FCB = F%*%diag(C_list[[1]][,1])%*%t(B1)
BCB = B1%*%diag(C_list[[1]][,1])%*%t(B1)

overallB = FCB+t(FCB)+BCB
hist(overallB)
image(overallB)

plot(C_list[[1]][,2],ylim=c(-1,1))
for(i in 1:20){
  lines(C_list[[1]][,i])
}
