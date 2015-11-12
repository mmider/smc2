source('ibis.R')


# N.theta <- 2 * 1e6
# t.max <- 2 * 1e2
# y <- sample.y(t.max)
# y <- rep(3,t.max)
# 
# estimate.mean <- ibis(gauss.log.likelihood, sample.gauss, kernel.move, N.theta, y, t.max, ESS)
# 
# hist(estimate.mean, freq=F, ylim = c(0,5))
# sd(estimate.mean)
# 3/sqrt((t.max + 1))
# x.axis <- seq(2,4,by=0.01)
# y.axis <- dnorm(x.axis, mean = t.max/(t.max+1) * 3, sd = 3/sqrt((t.max + 1)))
# lines(x.axis, y.axis)

N.theta <- 1e4
t.max <- 1e2
y <- sample.y.m(t.max)
y <- matrix(rep(3,3*t.max), ncol = 3)

estimate.mean <- ibis(gauss.log.likelihood.m, sample.gauss.m, kernel.move.m, N.theta, y, t.max, ESS)
library(rgl)

plot3d(estimate.mean[,1],estimate.mean[,2], estimate.mean[,3])



