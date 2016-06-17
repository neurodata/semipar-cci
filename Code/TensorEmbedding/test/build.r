library(devtools)

setwd("~/git/semipar-cci/Code")
build('TensorEmbedding/')
install.packages("TensorEmbedding_1.0.tar.gz", repos = NULL, type = "source")

require("RcppArmadillo")

require("TensorEmbedding")
