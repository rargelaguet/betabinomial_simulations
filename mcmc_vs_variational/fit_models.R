##------------------------
# Fit the BB regression model to these two datasets independently and store the results.
# NOTE: There is no downstream analysis performed in this script.
##------------------------

# library(logitnorm)
library(rstan)
# library(gtools)
# library(VGAM)
library(argparse)

###################
## Load settings ##
###################

if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  stop("Fix this")
} else if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/betabinomial_simulations/settings.R")
  io$data_dir <- "/Users/ricard/data/betabinomial/diff_residual_disp/synthetic/"
  io$out_dir <- "/Users/ricard/data/betabinomial/mcmc_vs_vb"
} else{
  stop("Computer not recognised")
}

################################
## Initialise argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('-r',  '--replicates',  type="integer",                help='Number of replicates (integer)')
p$add_argument('-n',  '--ncells',      type="integer",                help='Number of cells (integer)')
p$add_argument('--seed',               type="integer",  default=42,   help='Random seed')
p$add_argument('--test',               action="store_true",           help='Testing mode?')
p$add_argument('--mcmc',               action="store_true",           help='Use MCMC? (default is VB)')

# Parse arguments
args <- p$parse_args(commandArgs(TRUE))
r <- args$replicates    # Number of replicates
N_cells <- args$ncells  # Number of cells

args$test <- TRUE

# Testing mode
if (isTRUE(args$test)) {
  r <- 1
  N_cells <- 100
}

# Set replication directory
io$data_dir <- paste0(io$data_dir, "rep", r, "/")

cat("Replicate:", r, " Cells:", N_cells, "\n")

#########################
## Define STAN options ##
#########################

# Output file
io$stan_model_file <- paste0(io$lib_dir, "stan/bbreg.stan")

# Write Stan model to the hard disk
#rstan::rstan_options(auto_write = TRUE)
Sys.setenv(USE_CXX14 = 1)

###############
## Load data ##
###############

cat("Read simulated data..\n")

N_feat <- 300
N_cpgs <- 15
sim_dt <- readRDS(paste0(io$data_dir, "data_bbreg_feat", N_feat, "_cells", N_cells, "_cpgs", N_cpgs, ".rds"))

#################################
## Initialize model parameters ##
#################################

# Model parameters
m_wmu <- rep(0, NCOL(sim_dt$sim_dt_A$X))
s_wmu <- diag(0.5, NCOL(sim_dt$sim_dt_A$X))
s_mu <- 1
L_fit <- 4
m_wgamma <- rep(0, L_fit)
s_wgamma <- diag(0.5, L_fit)
a_sgamma <- 2
b_sgamma <- 3
rbf_c <- 1

# Number of iterations
if (isTRUE(args$test)) {
  iter <- 1000
} else {
  iter <- 100000
}

# Number of chains
if (isTRUE(args$test)) {
  chains <- 1
} else {
  chains <- 4
}

# Number of cores
#n_cores <- parallel::detectCores() - 1
n_cores <- chains

###############
## Fit model ##
###############

# Fit group A
cat("Doing inference for group A...\n")
fit_bbreg_A <- infer_bbreg(Y = sim_dt$sim_dt_A$Y, stan_model_file = io$stan_model_file, X = sim_dt$sim_dt_A$X,
                           L = L_fit, use_mcmc = args$mcmc, m_wmu = m_wmu, s_wmu = s_wmu, s_mu = s_mu,
                           m_wgamma = m_wgamma, s_wgamma = s_wgamma, a_sgamma = a_sgamma, b_sgamma = b_sgamma,
                           rbf_c = rbf_c, use_eb = TRUE, init_using_eb = TRUE, iter = iter,
                           chains = chains, n_cores = n_cores)

# Fit group B
cat("Doing inference for group B...\n")
fit_bbreg_B <- infer_bbreg(Y = sim_dt$sim_dt_B$Y, stan_model_file = io$stan_model_file, X = sim_dt$sim_dt_B$X,
                           L = L_fit, use_mcmc = args$mcmc, m_wmu = m_wmu, s_wmu = s_wmu, s_mu = s_mu,
                           m_wgamma = m_wgamma, s_wgamma = s_wgamma, a_sgamma = a_sgamma, b_sgamma = b_sgamma,
                           rbf_c = rbf_c, use_eb = TRUE, init_using_eb = TRUE, iter = iter,
                           chains = chains, n_cores = n_cores)

#################
## Save output ##
#################

obj <- list(fit_bbreg_A = fit_bbreg_A, fit_bbreg_B = fit_bbreg_B, sim_dt = sim_dt)
outfile <- paste0(io$out_dir, "bbreg_feat", N_feat, "_cells", N_cells, "_cpgs", N_cpgs, "_mcmc", args$mcmc, ".rds")
saveRDS(object = obj, file = outfile)
