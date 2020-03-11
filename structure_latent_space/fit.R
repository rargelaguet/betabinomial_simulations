##################
## iid data set ##
##################

# Simulate iid binomial data set
simulated.iid <- simulate.iid(n_samples = N, n_features = D, Ntotal.max = Ntotal.max)
Nmet <- simulated.iid$Nmet
Ntotal <- simulated.iid$Ntotal

# Fit Beta-binomial model for every feature
rho.iid <- rep(NA,D)
for (d in 1:D) {
  rho.iid[[d]] <- tryCatch(
    Coef( VGAM::vglm(cbind(Nmet[,d], Ntotal[,d]-Nmet[,d]) ~ 1, betabinomial) )["rho"], error = function(x) return(NA))
}


#########################
## structured data set ##
#########################

# Simulate structured binomial data set
simulated.struct <- simulate.structured(n_samples = N, n_features = D, n_factors = K, Ntotal.max = Ntotal.max)
Nmet <- simulated.struct$Nmet
Ntotal <- simulated.struct$Ntotal

# Fit Beta-binomial model for every feature
rho.struct <- rep(NA,D)
for (d in 1:D) {
  rho.struct[[d]] <- tryCatch(
    Coef( VGAM::vglm(cbind(Nmet[,d], Ntotal[,d]-Nmet[,d]) ~ 1, betabinomial) )["rho"], error = function(x) return(NA))
}

