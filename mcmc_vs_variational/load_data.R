N <- 100  # number of samples
s <- 10   # number of binomial trials

Y <- list()
for (i in 1:opts$ntrials) {
  
  # Sample theta
  theta <- runif(1, min=0.1, max=0.9)
  
  # Sample binomial observations
  Y[[i]] <- rbinom(N, s, theta)
}

