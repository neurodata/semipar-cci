require("monomvn")

n = 10
p = 50
X = matrix(rnorm(n * p, 0, 1), n, p)

beta = c(2, rep(0, p - 1))

y = X %*% beta + rnorm(n, 0, 0.1)

fit = bhs(X, y, icept = F)

plot(colMeans(fit$beta))


#compute the constant
X2 = t(X) %*% X
XY = t(X) %*% y
y2 = sum(y^2)

#initialize the parameters
sigma = 1
nu = rep(0.1, p)
lambda = rep(0.1, p)
eps = 0.1
tau = 0.1
trace_sigma =numeric()

for (i in 1:1000)
{
    #update beta
    v = solve(X2 + diag(1 / lambda / tau)) * sigma
    m = v %*% (XY/sigma)
    beta = t(chol(v)) %*% rnorm(p) + m

    #update lambda and nu
    lambda = 1 / rgamma(p, 1, 1 / nu + beta ^ 2 / 2 / tau/sigma)
    nu = 1 / rgamma(p, 1, 1 + 1 / lambda)

    #update tau and eps
    tau = 1 / rgamma(1, (1 + p) / 2, 1 / eps + sum(beta ^ 2 / lambda) / 2 / sigma)
    eps = 1 / rgamma(1, 1, 1 / tau + 1)

    #update sigma
    Xbeta = X%*%beta
    sigma = 1 / rgamma(1, (n + p) / 2 - 1, (sum((y - Xbeta) ^ 2) + sum(beta ^ 2 / tau / lambda)) / 2)
    trace_sigma = c(trace_sigma, sigma)   
}

plot(beta)
ts.plot(trace_sigma[200:1000])