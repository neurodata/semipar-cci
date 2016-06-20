require(tensr)

P<- list()
P[[1]]<- c(0.42,0.42)
P[[2]]<- c(0.42,0.5)
rho<- c(0.6,0.4)

n<- 20
m<- 10

C<- rbinom(n,1,rho[1])+1

U<- do.call("rbind",lapply(C,function(x){P[[x]]}))

P<- U%*%t(U)


modes <- c(n, n, m)
tnsr <- array(stats::rnorm(prod(modes)), dim = modes)


for(k in 1:m){
  A<- runif(n*n)<P
  #make it symmetric
  temp <- A[upper.tri(A)]
  A<- t(A)
  A[upper.tri(A)]<-temp
  tnsr[,,k]<-A
}


# avgA<- matrix(0,n,n)
# for(k in 1:m){
#   avgA <- avgA + tnsr[,,k] 
# }
# estP <- avgA/m
# image(P,zlim = c(0,1))
# image(estP,zlim = c(0,1))

d<- 50

R<-matrix( runif(d*n),n,d)
# 
# for(i in 1:10000){
#   L<- t(solve((t(R)%*%R),(t(R)%*%A)))
#   R<- t(solve((t(L)%*%L),(t(L)%*%A)))
# }
# 
# image(L%*%t(R),zlim=c(0,1))
# image(A,zlim=c(0,1))
# 
# cor(c(L%*%t(R)),c(A))

# 

# 
# L0<- cbind(L,matrix(0,n,n-d))
# diff = A- L%*%t(L)
# 
# dL<- -4* diff %*%L
# 
# dLnum<- L*0
# for(i in c(1:length(L))){
#   Lnum <- L
#   Lnum[i] = Lnum[i] +1E-12
#   dLnum[i] =(Q(A,Lnum)- Q(A,L))/1E-12
# }
# 
# plot(dL[1,])
# lines(dLnum[1,])
# 
# cor(dL[1,],dLnum[1,])

L<- matrix( runif(d*n),n,d)

Q<- function(L,A){
   sum((A-L%*%t(L))^2)
}

gradient<- function(L,A){
  diff <- A- L%*%t(L)
  - 4* diff%*%L
}
# 
# 
# line_search<- function(L,A, delta = 0.5){
#   
#   if(delta>1E-8){
#     diff <-  A- L%*%t(L)
#     cur_q <- sum(diff^2)
#     grad <-   -4* diff%*%L
#     
#     new_L <- L - delta * gradient(L,A)
#     new_diff <-  A- new_L%*%t(new_L)
#     new_q <- sum(new_diff^2)
#     
#     if(!is.na(new_q)){
#       if(new_q< cur_q)
#         L <- new_L
#         print(c(new_q))
#         L <- line_search(new_L,A,delta/2)
#     }
#   }   
#   return(L)
# }
# 
# for(i in 1:100)
#  L<-line_search(L,A,0.001)


