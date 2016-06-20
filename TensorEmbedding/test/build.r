library(devtools)

setwd("~/git/semipar-cci/")
build('TensorEmbedding/')
install.packages("TensorEmbedding_1.0.tar.gz", repos = NULL, type = "source")

require("RcppArmadillo")

require("TensorEmbedding")
