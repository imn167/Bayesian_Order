library(ArchaeoPhases)
## interval_credible : from the function arkhe::interval_credible
library(truncnorm)
truncnorm <- function(a, n) {
  x = c()
  while (length(x) < n) {
    g = rnorm(1)
    if (abs(g)> a) {
      x = c(x, g)
    }
    
  }
  return(x)
}

X = truncnorm(1, 10000)
kde = density(X)
plot(kde)
sorted = sort(kde$y, decreasing = T, index.return = T)
sorted$ix
N = length(X)
quant <- floor(N*.05)
ind <- sorted$ix[sorted$x>sorted$x[quant]]
sim_HPD = X[ind]
plot(range(sim_HPD), c(0, 1), type = 'n', xlab = "x", ylab = "y")
points(sim_HPD, rep(.5, quant-1), pch = 2, cex = .1)
sort(sim_HPD)
