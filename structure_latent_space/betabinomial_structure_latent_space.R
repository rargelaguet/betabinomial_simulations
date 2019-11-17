########################################################################################################
## Show that the existence of subpopulations in DNA methylation data leads to overdispersed estimates ##
########################################################################################################

library(VGAM)
require(MASS)
require(data.table)
require(purrr)
require(ggpubr)
require(RColorBrewer)

source("/Users/ricard/betabinomial_simulations/structure_latent_space/utils.R")

#####################
## Define settings ##
#####################

D = 50             # Number of features
N = 100              # Number of samples
Ntotal.max <- 25    # Maximum number of counts per feature

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

K = 5 # Number of latent factors

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


################################
## Plot distributions of Nmet ##
################################

hist(simulated.iid$Nmet, col=rgb(0,0,1,1/4))
hist(simulated.struct$Nmet, col=rgb(1,0,0,1/4), add=T)

####################################
## Plot distribution of variances ##
####################################

# # gaussian var
# var.iid <- apply(simulated.iid$Nmet / Ntotal.max,2,var)
# var.struct <- apply(simulated.struct$Nmet / Ntotal.max,2,var)
# dt <- data.table(iid = var.iid, struct = var.struct) %>%
#   melt(variable.name="type", value.name="variance")

# binomial var
rate.iid <- simulated.iid$Nmet / Ntotal.max
rate.struct <- simulated.struct$Nmet / Ntotal.max
var.iid <- rate.iid * (1-rate.iid) * Ntotal.max
var.struct <- rate.struct * (1-rate.struct) * Ntotal.max
dt <- data.table(iid = as.numeric(var.iid), struct = as.numeric(var.struct)) %>%
  melt(variable.name="type", value.name="variance")

p1 <- ggboxplot(dt, x="type", y = "variance", fill = "type") +
  scale_fill_brewer(palette = "Dark2") +
  # coord_cartesian(xlim=c(0,0.15), ylim=c(0,60)) +
  labs(x="", y="Binomial variance")

###################################################
## Plot distribution of overdispersion estimates ##
###################################################

dt <- data.table(
  iid = rho.iid,
  struct = rho.struct
) %>% melt(variable.name="type", value.name="rho")

# ggdensity(dt, x="rho", y = "..density..", fill = "type") +
#   scale_fill_brewer(palette = "Dark2") +
#   # coord_cartesian(xlim=c(0,0.15), ylim=c(0,60)) +
#   labs(x="Overdispersion (rho)", y="Density")

p2 <- ggboxplot(dt, x="type", y = "rho", fill = "type") +
  scale_fill_brewer(palette = "Dark2") +
  # coord_cartesian(xlim=c(0,0.15), ylim=c(0,60)) +
  labs(x="", y="Beta-binomial overdispersion ")

cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
