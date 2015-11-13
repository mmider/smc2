require(MCMCpack)

log.prior.theta <- function(x, dim.observations = 3){
  log.pr <- log(dinvgamma(x, 1,1)^dim.observations)
  return(log.pr)
}

sample.kernel <- function(thetas, x.t, a.t, weights.t,y.t){
  # no move at all for now, just stay where you are
  return(list(thetas = thetas, x.t = x.t, a.t = a.t, weights.t = weights.t))
}

ESS <- function(log.weights){
  N <- length(log.weights)
  weights <- exp(log.weights)
  R <- sum(weights^2)
  ess <- (sum(weights))^2/R
  return((ess < N/2) || (sum(weights^2) < 1e-180))
}

generate.particle <- function(prev.x, theta){
  new.x <- prev.x + matrix(rnorm(ncol(prev.x) * nrow(prev.x), sd = theta[1]), ncol = ncol(prev.x))
  return(new.x)
}

g.density <- function(y, x ,theta){
  return(apply(dnorm(y-x, sd = theta[2]),1, prod))
}
x.sample.prior <- function(N.x, theta, dim = 3){
  return(matrix(rnorm(N.x * dim, sd = theta[1]), ncol = dim))
}
theta.sample.prior <- function(N.x){
  return(matrix(runif(N.x * 2,min = 1, max = 6), ncol = 2))
}

particle.filter.step <- function(theta, current.y, current.t, previous.x = NULL, N.x = NULL){
  # if current.t == 1, then N.x needs to give the number of x particles to sample
  # else we are continuing the moves, so N.x is whatever - it will not matter -
  # though previous.x - x particles from the previous step must be supplied in this case
  
  if (current.t == 1){
    x <- x.sample.prior(N.x, theta)
    weights <- g.density(current.y, x, theta)
    temp <- list(current.y, x[1:20,], theta)
  }
  else {
    x <- generate.particle(previous.x, theta)
    weights <- g.density(current.y,x, theta)
    temp2 <-list(current.y, x, theta)
  }
  W <- weights/sum(weights)
#   if (any(is.na(W))){
#     print(temp)
#     if (current.t != 1){
#       print("****************************")
#       print(temp2)
#     }
#   }
  a <- sample(1:N.x, prob = W, replace = T)
  return(list(x = x, weights = weights, a = a))
}

smc2 <- function(N.theta, N.x, y, t.max, degen.cond){  
  # initialization
  y <- as.matrix(y)
  thetas <- as.matrix(theta.sample.prior(N.theta))
  log.weights.theta <- rep(0, N.theta)
  x.history <- array(NA, c(N.x, dim.x, N.theta, t.max))
  a.history <- array(NA, c(N.x, N.theta, t.max))
  weights.history <- array(NA, c(N.x,N.theta, t.max))
  for (t in 1:t.max){
    # update weight:
    for (theta.index in 1:N.theta){
      if (t == 1){
        previous.x <- NULL        
      }
      else {
        previous.x <- x.history[a.history[,theta.index,t-1],,theta.index,t-1]
      }
      new.particle <- particle.filter.step(theta = thetas[theta.index,], current.y = y[t, ], current.t = t, previous.x, N.x = N.x)
      x.history[,,theta.index, t] <- new.particle$x
      a.history[,theta.index,t] <- new.particle$a
      weights.history[,theta.index,t] <- new.particle$weights
      
      
      log.evidence <- log(sum(new.particle$weights)/N.x)
      log.weights.theta[theta.index] <- log.weights.theta[theta.index] + log.evidence
    }
    # check for degeneracy condition after all thetas at time t have been updated:
    if (degen.cond(log.weights.theta)){
      probs <- exp(log.weights.theta)/sum(exp(log.weights.theta))
      sample.indices <- sample(1:N.theta, size = N.theta, prob = probs, replace = T)
      
      # this function does nothing at this moment:
      new.parameter.list <- sample.kernel(thetas[sample.indices,], x.history[,,sample.indices,t], a.history[,sample.indices,t], weights.history[,sample.indices,t], y[t,])
      thetas <- new.parameter.list$thetas
      x.history[,,,t] <- new.parameter.list$x.t
      a.history[,,t] <- new.parameter.list$a.t
      weights.history[,,t] <- new.parameter.list$weights.t
      
      log.weights.theta <- rep(0,N.theta)
    }
    
  }
  
  # final resampling:
  return(list(thetas = thetas, log.weights.theta = log.weights.theta,
              x.history = x.history, a.history = a.history,
              weights.history = weights.history))
}
