require(MCMCpack)
require(mvtnorm)
require(Rcpp)
sourceCpp('smc2.cpp')


gaussianModel <- function(){
  true.sigma = 0.5
  true.tau = 1
  
  prior.shape.sigma <- 1
  prior.rate.sigma <- 1
  prior.shape.tau <- 1
  prior.rate.tau <- 1

  dim.param <- 2
  dim.x <- 1
  
  log.prior.theta <- function(x){
    # parameters in matrix (num.obs by 2)
    # possibly need to introduce Jacobian here:
    if (is.vector(x)){
      x <- matrix(x, ncol = 2, byrow = T)
    }
    sigma.lp <- log(dinvgamma(x[,1], prior.shape.sigma, prior.rate.sigma))
    theta.lp <- log(dinvgamma(x[,2], prior.shape.tau, prior.rate.tau))
    return(sigma.lp + theta.lp)
  }
  
  generate.prior.theta <- function(n){
    sigma.samples <- rinvgamma(n, prior.shape.sigma, prior.rate.sigma)
    tau.samples <- rinvgamma(n, prior.shape.tau, prior.rate.tau)
    return(cbind(sigma.samples, tau.samples))
  }
  
  generate.synth.data <- function(t.max){
    noise.x <- rnorm(t.max,mean = 0, sd = true.sigma)
    noise.y <- rnorm(t.max, mean = 0, sd = true.tau)
    x <- cumsum(noise.x)
    y <- as.matrix(x + noise.y)
    return(y)
  }
  
  generate.prior.x <- function(n, theta){
    x <- rnorm(n, mean = 0, sd = theta[1])
    return(x)
  }
  
  g <- function(y,x, theta){
    likelihood <- dnorm(y, mean = x, sd = theta[2])
    return(likelihood)
  }
  
  generate.x.particle <- function(x, theta){
    new.x <- x + as.matrix(rnorm(nrow(x), sd = theta[1]))
    return(new.x)
  }
  
  ESS <- function(log.weights){
    N <- length(log.weights)
    weights <- exp(log.weights)
    R <- sum(weights^2)
    ess <- (sum(weights))^2/R
    return((ess < N/2) || (sum(weights^2) < 1e-180))
  }
  
  sample.kernel.id <- function(thetas, x, a, weights, log.evidence, y,...){
    return(list(thetas = thetas, x = x, a = a, weights = weights, log.evidence = log.evidence))
  }
  
  sample.kernel.gauss <- function(thetas, x, a, weights, log.evidence, y,t){
    Nx <- dim(x)[1]
    mean.estim <- colMeans(thetas)
    cov.estim <- cov(thetas)
    N <- nrow(thetas)
    thetas.prop <- rmvnorm(N, mean = mean.estim, sigma = cov.estim)
    
    final.state <- list(thetas = thetas, x = x, a = a, weights = weights, log.evidence = log.evidence)
        
    for (theta.index in 1:nrow(thetas.prop)){
      theta <- thetas.prop[theta.index,]
      if (any(theta <= 0))
        next
      theta.old <- thetas[theta.index,]
      temp <- particle_filter(theta[1], theta[2], y,Nx)
      x.prop <- temp[,1:t]
      weights.prop <- temp[,(t+1):(2*t)]
      a.prop <- temp[,(2*t+1):(3*t-1)]
      if (any(is.na(log(colMeans(weights.prop))))){
        print("yoyo")
        print(weights.prop)
      }
      log.evidence.prop <- log(colMeans(weights.prop))
      log.Z.estim <- sum(log.evidence.prop)
      log.Z.old.estim <- sum(log.evidence)
      Tbackward <- dmvnorm(theta.old, mean = mean.estim, sigma = cov.estim, log = T)
      Tforward <- dmvnorm(theta, mean = mean.estim, sigma = cov.estim, log = T)
      alpha <- exp(log.prior.theta(theta) + log.Z.estim + Tbackward -
                     log.prior.theta(theta.old) - log.Z.old.estim - Tforward)
      if (runif(1) < alpha){
        final.state$thetas[theta.index,] <- theta
        final.state$x[,theta.index,] <- x.prop
        final.state$a[,theta.index,] <- a.prop
        final.state$weights[,theta.index,]<-weights.prop
        final.state$log.evidence[theta.index,]<-log.evidence.prop
      }
    }
    return(final.state)
  }
  
  return(list(dim.param = dim.param, dim.x = dim.x, log.prior.theta = log.prior.theta, generate.prior.theta = generate.prior.theta,
              generate.synth.data = generate.synth.data, generate.prior.x = generate.prior.x, g = g,
              generate.x.particle = generate.x.particle, degen.cond = ESS, sample.kernel = sample.kernel.gauss))
}

smc2 <- function(N.theta, N.x, t.max, model){  
  # initialization
  if (model == "simpleGaussians"){
    mod.spec <- gaussianModel()
  }
  log.prior.theta <- mod.spec$log.prior.theta
  generate.prior.theta <- mod.spec$generate.prior.theta
  generate.prior.x <- mod.spec$generate.prior.x
  generate.x.particle <- mod.spec$generate.x.particle
  degen.cond <- mod.spec$degen.cond
  sample.kernel <- mod.spec$sample.kernel
  g <- mod.spec$g
  
  y <- mod.spec$generate.synth.data(t.max)
  theta.dim <- mod.spec$dim.param
  dim.x <- mod.spec$dim.x
  thetas <- as.matrix(generate.prior.theta(N.theta))
  log.weights.theta <- rep(0, N.theta)
  
  x.history <- array(NA, c(N.x, dim.x, N.theta, t.max))
  a.history <- array(NA, c(N.x, N.theta, t.max-1))
  weights.history <- array(NA, c(N.x,N.theta, t.max))
  log.evidence.history <- array(NA, c(N.theta, t.max))
  
  # first step:
  for (theta.index in 1:N.theta){
    safety.condition <- T
    while (safety.condition){
      theta <- thetas[theta.index,]
      x <- generate.prior.x(N.x, theta)
      weights <- g(y[1,], x, theta)
      if (sum(weights)>1e-3){
        safety.condition <- F
      }
      log.evidence.history[theta.index,1] <- log(mean(weights))
      x.history[,,theta.index,1] <- x
      weights.history[,theta.index,1] <- weights
    }
  }  
  
  for (t in 2:t.max){
    # update weight:
    for (theta.index in 1:N.theta){
      theta <- thetas[theta.index,]
      weights <- weights.history[,theta.index,t-1]
      W <- weights/sum(weights)
      if (any(is.na(W))){
        print("heyhey")
        print(weights.history[,theta.index,])
      }
      a <- sample(1:N.x, prob = W, replace = T)
      x.resampled <- as.matrix(x.history[a,,theta.index,t-1])
      x <- generate.x.particle(x.resampled, theta)
      weights <- g(y[t,], x, theta)
      
      x.history[,,theta.index, t] <- x
      a.history[,theta.index,t-1] <- a
      weights.history[,theta.index,t] <- weights
      log.evidence <- log(mean(weights))
      log.evidence.history[theta.index,t] <- log.evidence
      log.weights.theta[theta.index] <- log.weights.theta[theta.index] + log.evidence
    }
    # check for degeneracy condition after all thetas at time t have been updated:
    if (degen.cond(log.weights.theta)){
      probs <- exp(log.weights.theta)/sum(exp(log.weights.theta))
      sample.indices <- sample(1:N.theta, size = N.theta, prob = probs, replace = T)
      
      # this function does nothing at this moment:
      new.parameters <- sample.kernel(thetas[sample.indices,], x.history[,,sample.indices,1:t],
                                          a.history[,sample.indices,1:(t-1)], weights.history[,sample.indices,1:t],
                                          log.evidence.history[sample.indices, 1:t], y[1:t,], t)
      thetas <- new.parameters$thetas
      x.history[,,,1:t] <- new.parameters$x
      a.history[,,1:(t-1)] <- new.parameters$a
      weights.history[,,1:t] <- new.parameters$weights
      log.evidence.history[,1:t] <- new.parameters$log.evidence
      
      log.weights.theta <- rep(0,N.theta)
    }
    
  }
  
  return(list(thetas = thetas, log.weights.theta = log.weights.theta,
              x = x.history, a = a.history,
              weights = weights.history, log.evidence = log.evidence.history))
}
