require("rstiefel")
require("msm")

symmetrize <- function(X)
{
    X = t(X)
    X[lower.tri(X)] <- t(X)[lower.tri(X)]
    X
}

stiefel_diagonalize <- function(Y_list, r, diagD = T)
{

    n = nrow(Y_list[[1]])
    m = length(Y_list)

    ub_list = list()
    lb_list = list()
    for (j in 1:m)
    {
        mat_ub = matrix(Inf, n, n)
        mat_lb = matrix( - Inf, n, n)

        mat_lb[Y_list[[j]] == 1] = 0
        mat_ub[Y_list[[j]] == 0] = 0
        ub_list[[j]] = mat_ub
        lb_list[[j]] = mat_lb
    }

    # k= 10
    tau = 100


    avgZ = matrix(rowMeans(matrix(unlist(Y_list), n * n)), n)

    # Z_stacked = do.call("cbind",Z_list)

    diag(avgZ) = 0.99
    avgZ[avgZ == 1] = 0.99
    avgZ[avgZ == 0] = 0.01
    svd1 <- svd(qnorm(avgZ), nu = r, nv = r)
    U0 <- svd1$u


    zeroB <- matrix(0, r, r)
    zeroA <- matrix(0, n, n)
    k = 1

    U_list = list()
    for (j in 1:m)
    {
        U_list[[j]] = rbmf.matrix.gibbs(zeroA, zeroB, k * U0, U0)
    }

    Z_list = list()
    for (j in 1:m)
    {
        Z_list[[j]] = Y_list[[j]] - 0.5
    }


    trace_D <- numeric()
    for (i in 1:200)
    {

        #update D
        D_list = list()
        for (j in 1:m)
        {
            # vinvm = apply(U_list[[j]], 2, function(x){  sum(x* (Z_list[[j]]%*%x))/2})

            if (diagD)
            {
                vinvm = apply(U0, 2, function(x) { sum(x * (Z_list[[j]] %*% x)) / 2 })
                v = 2 * tau / (2 + tau)
                mean_D = v * vinvm
                D_list[[j]] = diag(rnorm(r, mean_D, sqrt(v)))
            } else
            {

                matD = matrix(0, r, r)
                mat_lowertri_idx = lower.tri(matD)

                vinvm = (t(U0) %*% Z_list[[j]] %*% U0)
                diag(vinvm) = diag(vinvm) / 2
                vij = 1 * tau / (1 + tau)
                mean_Dij = vij * vinvm
                matD[mat_lowertri_idx] = rnorm(sum(mat_lowertri_idx), mean_Dij[mat_lowertri_idx], sqrt(vij))
                matD = symmetrize(matD)

                vii = 2 * tau / (2 + tau)
                mean_Ddiag = vii * diag(vinvm)

                diag(matD) = rnorm(r, mean_Ddiag, sqrt(vii))
                D_list[[j]] = matD

            }
        }

        #update Z

        Z_list = list()
        for (j in 1:m)
        {
            matZ = matrix(0, n, n)

            mean_Z = U0 %*% D_list[[j]] %*% t(U0)


            lowertri_idx = lower.tri(mean_Z)
            mean_Z_lowertri = mean_Z[lowertri_idx]
            lb_lowertri = lb_list[[j]][lowertri_idx]
            ub_lowertri = ub_list[[j]][lowertri_idx]
            z = rtnorm(length(mean_Z_lowertri), mean = mean_Z_lowertri, sd = 1, lb_lowertri, ub_lowertri)

            matZ[lowertri_idx] = z
            matZ = symmetrize(matZ)
            diag(matZ) = rnorm(n, diag(mean_Z), sqrt(2))
            Z_list[[j]] = matZ
        }


        # if(FALSE){
        #   
        #   compLoglik<- function(U0){
        #     cur_loglik = sapply(c(1:m),function(j){
        #     U0ZU0= sapply(1:r, function(i)sum((U0[,i]) * (Z_list[[j]]%*% (U0[,i]))))
        #     sum(U0ZU0* D_list[[j]])/2
        #   })
        #    sum(cur_loglik) 
        #   }
        #   
        #   # compLoglik(U0)
        # 
        #   k=1
        #   #update U
        #   for(j in 1:m){
        #     U_list[[j]] = rbmf.matrix.gibbs(Z_list[[j]]/2, diag(D_list[[j]]), k*U0, U_list[[j]])
        #   }
        #   
        #   #update U0
        #   k=100
        #   U_sum = matrix(rowSums(matrix(unlist(U_list),n*r)),n)
        # 
        #   C= k * U_sum
        #   svdC = svd(C)
        #   propU0 = rbmf.matrix.gibbs(zeroA, zeroB, svdC$u%*%diag(svdC$d), U0%*%svdC$v) %*% t(svdC$v)
        #   
        #   compLoglik(propU0)   -    compLoglik(U0)
        # 
        #   
        #   # U0 = rbmf.matrix.gibbs(zeroA, zeroB, k * U_sum, U0)
        # 
        #   
        # }

        trace_D = rbind(trace_D, unlist(D_list))
        print(i)
    }

    return(list("D_list" = D_list, "trace_D" = trace_D, "U0" = U0))
}

getPEst <- function(U0, D_list)
{
    lapply(D_list, function(d) { pnorm(U0%*%d%*%t(U0)) })
}