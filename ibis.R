update.log.weights <- function(log.likelihood, theta, y, t, log.weights){
  if (t == 1){
    log.weights <- log.likelihood(theta, y[t,])
  }
  else
    log.weights <- log.weights + log.likelihood(theta, y[t,], y[t-1,])
  return(log.weights)
}

ibis <- function(log.likelihood, sample.prior, sample.kernel, N.theta, y, t.max, degen.cond){
  # initialization
  thetas <- as.matrix(sample.prior(N.theta))
  log.weights <- rep(0, N.theta)
  y <- as.matrix(y)
  for (t in 1:t.max){
    # update weight:
    log.weights <- update.log.weights(log.likelihood, thetas, y, t, log.weights)
    
    # check for degeneracy condition:
    if (degen.cond(log.weights)){
      probs <- exp(log.weights)/sum(exp(log.weights))
      sample.indices <- sample(1:N.theta, size = N.theta, prob = probs, replace = T)
      
      thetas <- sample.kernel(thetas[sample.indices,],log.likelihood, y[t,])
      log.weights <- rep(0,N.theta)
    }
  }
  
  # final resampling:
  probs <- exp(log.weights)/sum(exp(log.weights))
  thetas <- thetas[sample(1:N.theta, size = N.theta, prob = probs, replace = T),]
  return(thetas)
}

####
# Multivariate setting:
####
require(mvtnorm)
gauss.log.likelihood.m <- function(theta, y.1, y.2=NULL){
  ll <- dmvnorm(theta, mean = y.1, sigma = diag(rep(3,3)),log = T)
  return(ll)
}

sample.gauss.m <- function(n){
  y <- rmvnorm(n, mean = rep(0,3), sigma = diag(rep(3,3)))
  return(y)
}

sample.y.m <- function(n){
  y <- rmvnorm(n, mean = rep(3,3), sigma = diag(rep(3,3)))
  return(y)
}

kernel.move.m <- function(x,log.lik,y){
  x.prop <- x + rmvnorm(nrow(x), sigma = diag(rep(0.1,3)))
  alpha <- exp(log.lik(x.prop,y) - log.lik(x,y))
  accepted <- runif(nrow(x)) < alpha
  x.new <- x.prop * accepted + x * (1-accepted)
  return(x.new)
}

# #####
# # Univariate setting:
# #####
# gauss.log.likelihood <- function(theta, y.1, y.2=NULL){  
#   ll <- dmvnorm(y.1, mean = theta, sd = 3,log = T)
#   return(ll)
# }
# 
ESS <- function(log.weights){
  N <- length(log.weights)
  weights <- exp(log.weights)
  R <- sum(weights^2)
  ess <- (sum(weights))^2/R
  return((ess < N/2) || (sum(weights^2) < 1e-180))
}
# 
# sample.gauss <- function(n){
#   y <- rnorm(n, mean = 0, sd = 3)
#   return(y)
# }
# 
# sample.y <- function(n){
#   y <- rnorm(n, mean = 3, sd = 3)
#   return(y)
# }
# 
# kernel.move <- function(x,log.lik,y){
#   n <- length(x)
#   x.cur <- x
#   x.prop <- x.cur + rnorm(n, sd = 0.1)
#   alpha <- exp(log.lik(x.prop,y) - log.lik(x.cur,y))
#   accepted <- runif(n) < alpha
#   x.new <- x.prop * accepted + x.cur * (1-accepted)
#   return(x.new)
# }
# 
# 
# 
