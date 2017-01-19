F<- matrix(rnorm(100),10)



i<- 3

C<- rnorm(10)

F2<- F%*% diag(C) %*% t(F)
  
F[i,] %*%diag(C)%*%t(F[-i,])
F2[i,-i]

# length(F[i,])
# dim(t(F[-i,]))

t(t(F[-i,]) * F[i,]) %*% C
