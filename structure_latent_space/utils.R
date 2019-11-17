matrix.please <- function(x) {
  m<-as.matrix(x[,-1])
  rownames(m)<-x[[1]]
  m
}


simulate.iid <- function(n_features = 10, n_samples = 100, Ntotal.max=10) {
  
  D = n_features
  N = n_samples
  
  # Sample observations
  # theta <- rep(NA,D)
  # rho <- rep(NA,D)
  
  # Ntotal <- matrix(NA, nrow=N, ncol=D)
  Nmet <- matrix(NA, nrow=N, ncol=D)
  
  for (j in 1:D) {
    
    # Sample the methylation rate (per feature)
    # theta[i] <- rbeta(n=1, shape1=1, shape2=1)
    theta <- rbeta(n=1, shape1=1, shape2=1)
    
    # Sample the overdispersion
    # rho[i] <- rbeta(n=1, shape1=1, shape2=3)
    # rho[i] <- 0
    
    for (i in 1:N) {
      
      # Sample the total number of reads for feature i in cell j
      # Ntotal[i,j] <- Ntotal.max
      
      # Sample the number of methylated reads for feature i in cell j
      # Nmet[i,j] <- VGAM::rbetabinom(1, Ntotal[i,j], prob = theta[i], rho = rho[i])
      # Nmet[i,j] <- VGAM::rbetabinom(1, Ntotal[i,j], prob = theta[i,j], rho = 0)
      Nmet[i,j] <- rbinom(n=1, size=Ntotal.max, prob=theta)
      
    }
  }
  
  # stopifnot(all(Nmet<=Ntotal))
  
  return(list(Nmet=Nmet, Ntotal=Ntotal))
}

simulate.structured <- function(n_features=100, n_samples = 50, n_factors = 5, Ntotal.max = 10) {
  
  # set dimensionalities
  D = n_features
  N = n_samples
  K = n_factors
  
  # sample Factors
  Z <- matrix(rnorm(N*K, 0, 1), nrow=N, ncol=K)
  
  # sample Weights
  W <- vapply(seq_len(K), function(k) rnorm(D, mean=0, sd=1), numeric(D))

  # sample total number of reads    
  # Ntotal <- matrix(sample(x=Ntotal.max, size=N*D, replace=T), nrow=N, ncol=D)
    
  # sample number of methylated reads according to the binomial likelihood
  Nmet <- matrix(rep(NA,N*D), nrow=N, ncol=D)

  theta <- 1/(1+exp(-(Z%*%t(W))))
  for (i in 1:N) {
    for (j in 1:D) {
      Nmet[i,j] <- rbinom(n=1, size=Ntotal.max, prob=theta[i,j])
    }
  }
  
  # stopifnot(all(Nmet<=Ntotal))
    
  return(list(Nmet=Nmet, Ntotal=Ntotal))
}
