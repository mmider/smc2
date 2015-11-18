# set your working directory
require(ggplot2)
source('smc2v2.R')
priorfunction <- function(x){
  shape <- 1.00000 
  scale <- 1.00000
  return(scale**shape / gamma(shape) * x**(- shape - 1) * exp(-scale / x))
}
N.theta = 500
N.x = 250
t.max = 100
model = "simpleGaussians"

results <- smc2(N.theta, N.x, t.max, model)
data <- data.frame(sigma = results$thetas[,1],
                   tau = results$thetas[,2])
p <- ggplot(data, aes(x = sigma)) +geom_histogram(aes(y=..density..)) + 
  geom_density(fill='blue', alpha=0.5) + geom_vline(xintercept = 2, linetype = 2, size = 1) + 
  stat_function(fun = priorfunction, aes(colour = "prior"), linetype = 1, size = 1)
print(p)

q <- ggplot(data, aes(x = tau)) +geom_histogram(aes(y=..density..)) +
  geom_density(fill='blue', alpha=0.5) + geom_vline(xintercept = 2, linetype = 2, size = 1) +
  stat_function(fun = priorfunction, aes(colour = "prior"), linetype = 1, size = 1)
print(q)
