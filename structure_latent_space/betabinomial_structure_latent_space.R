source("/Users/ricard/betabinomial_simulations/structure_latent_space/utils.R")

#####################
## Define settings ##
#####################

D = 100             # Number of features
N = 100            # Number of samples
Ntotal.max <- 50       # Maximum number of counts per feature
K = 5              # Number of latent factors for the structured data set

################
## Fit models ##
################

source("/Users/ricard/betabinomial_simulations/structure_latent_space/fit.R")

################################
## Plot distributions of Nmet ##
################################

# hist(simulated.iid$Nmet, col=rgb(0,0,1,1/4))
# hist(simulated.struct$Nmet, col=rgb(1,0,0,1/4), add=T)

####################################################
## Plot variance estimates using a gaussian model ##
####################################################

# var.iid <- apply(simulated.iid$Nmet / Ntotal.max,2,var)
# var.struct <- apply(simulated.struct$Nmet / Ntotal.max,2,var)
# 
# dt <- data.table(iid = as.numeric(var.iid[1,]), struct = as.numeric(var.struct[1,])) %>%
#   melt(variable.name="type", value.name="variance")
# 
# p1 <- ggboxplot(dt, x="type", y = "variance", fill = "gray70") +
#   geom_jitter(size=1, color="black", alpha=0.5) +
#   stat_compare_means(method = "wilcox.test") +
#   labs(x="", y="Gaussian variance")

####################################################
## Plot variance estimates using a binomial model ##
####################################################

rate.iid <- simulated.iid$Nmet / Ntotal.max
rate.struct <- simulated.struct$Nmet / Ntotal.max
var.iid <- rate.iid * (1-rate.iid) * Ntotal.max
var.struct <- rate.struct * (1-rate.struct) * Ntotal.max

dt <- data.table(iid = as.numeric(var.iid[1,]), struct = as.numeric(var.struct[1,])) %>%
  melt(variable.name="type", value.name="variance")

p1 <- ggboxplot(dt, x="type", y = "variance", fill = "gray70") +
  geom_jitter(size=1, color="black", alpha=0.5) +
  stat_compare_means(method = "wilcox.test") +
  labs(x="", y="Binomial variance")

###################################
## Plot overdispersion estimates ##
###################################

dt <- data.table(
  iid = rho.iid,
  struct = rho.struct
) %>% melt(variable.name="type", value.name="rho")

p2 <- ggboxplot(dt, x="type", y = "rho", fill = "gray70", outlier.shape=NA) +
  geom_jitter(size=1, color="black", alpha=0.5) +
  stat_compare_means(method = "wilcox.test") +
  labs(x="", y="Beta-binomial overdispersion ")

################
## Joint plot ##
################

p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)

pdf(sprintf("%s/structured_latent_space.pdf",outdir), width = 8, height = 5, useDingbats = F)
print(p)
dev.off()