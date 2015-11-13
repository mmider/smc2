source('smc2.R')

N.theta <-  1e4
t.max <- 4 * 1e1
N.x <-  2 * 1e2
dim.x <- 3
sigma1 <- 3
sigma2 <- 0.5
y <- replicate(dim.x,cumsum(rnorm(t.max, sd = sigma1))) + matrix(rnorm(dim.x * t.max, sd = sigma2), ncol = dim.x)

estimated.sigmas <- smc2(N.theta, N.x, y, t.max, ESS)
colMeans(estimated.sigmas$thetas)
