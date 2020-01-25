##------------------------
# Generate synthetic data and store the results locally.
# NOTE: There is no downstream analysis performed in this script.
##------------------------

library(data.table)
# library(tidyverse)
# library(logitnorm)
# library(gtools)
# library(VGAM)

#########
## I/O ## 
#########

io <- list()
if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  io$lib_dir <- "~/Research/Projects/epigenetic-heterogeneity/code/bbreg/lib/"
  io$out_dir <- "~/datasets/epigenetic-heterogeneity/scalability/synthetic/"
} else if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/bbreg/settings.R")	
  io$out_dir <- "/Users/ricard/data/betabinomial/variational_tolerance/"
}

#############
## Options ##
#############

opts <- list()

# Number of replicates
opts$total_rep <- 1

# Number of cells
opts$N_cells <- c(100)

# Number of features
opts$N_feat <- c(500)

# Number of CpGs (per feature and cell? sampled from binomial?)
opts$N_cpgs <- 15

###################
## Model options ##
###################

# Parameters of the priors
w_mu <- c(-0.5, -1.2)
s_mu <- 0.5
w_gamma <- c(-1.2, -.3, 1.1, -.9)
s_gamma <- 0.3

# Number of basis functions
L <- 4

# Fraction of cells sampled per feature
cells_range <- c(0.4,1)

# Regression of mu (default is linear with negative slope)
use_mu_covariates <- FALSE

# Gamma vs mu regression (default is non-linear)
model_mean_disp <- TRUE

##############
## Simulate ##
##############

for (r in 1:opts$total_rep) {
  
  # Set replication directory
  rep_dir <- paste0(io$out_dir, "rep", r, "/")
  # cat("Directory ", rep_dir, "\n")
  if (!dir.exists(rep_dir)) { dir.create(rep_dir, recursive = TRUE) }
  
  for (c in opts$N_cells) {
    for (f in opts$N_feat) {
      
      cat("Running simulation\n")
      sim_dt <- simulate_bbreg(N_feat = f, N_cells = c, N_cpgs = opts$N_cpgs, use_mu_covariates = use_mu_covariates,
                               w_mu = w_mu, s_mu = s_mu, model_mean_disp = model_mean_disp, L = L, cells_range = cells_range,
                               w_gamma = w_gamma, s_gamma = s_gamma, seed = r)
      
      cat("Storing results\n")
      saveRDS(object = sim_dt, file = paste0(rep_dir, "data_bbreg_feat", f, "_cells", c, "_cpgs", opts$N_cpgs, ".rds"))
      
    }
  }
}


