library(devtools)
# .rs.restartR()

setwd("~/git/semipar-cci/")
build('RandomFactorTensor')
install.packages("RandomFactorTensor", repos = NULL, type = "source")

require("RandomFactorTensor")


A<- c(1:100)
g<- c(rep(0,3),rep(1,2),rep(2,20))

a=RandomFactorTensor::random_factor_decomp(A_r = A, n = 2, m = 25 ,g_r = g,rank=2)
