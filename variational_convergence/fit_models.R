suppressPackageStartupMessages(library(logitnorm))
suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(gtools))
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(argparse))

###################
## Load settings ##
###################

if (grepl("S34-R31YLVDR",Sys.info()['nodename'])) {
  stop("Fix this")
} else if (grepl("ricard",Sys.info()['nodename'])) {
  source("/Users/ricard/bbreg/settings.R")
  io$data_dir <- "/Users/ricard/data/betabinomial/variational_tolerance/data"
  io$outdir <- "/Users/ricard/data/betabinomial/variational_tolerance"
} else{
  stop("Computer not recognised")
}

dir.create(io$outdir, showWarnings = F)

################################
## Initialise argument parser ##
################################

p <- ArgumentParser(description='')
p$add_argument('--replicate',   type="integer",                help='replicates ID (integer)')
p$add_argument('--ncells',      type="integer",                help='Number of cells (integer)')
p$add_argument('--nfeatures',   type="integer",                help='Number of features (integer)')
p$add_argument('--test',        action="store_true",           help='Testing mode?')
p$add_argument('--tolerance',   type="double",                 help='Tolerance for the VB algorithm')
p$add_argument('--seed',        type="integer",  default=42,   help='Random seed')
args <- p$parse_args(commandArgs(TRUE))

## TEST ##
# args <- list()
# args$replicate <- 1
# args$ncells <- 100
# args$nfeatures <- 500
# args$test <- TRUE
# args$tolerance <- 0.01
# args$seed <- 42
## TEST ##

# Set replication directory
io$data_dir <- paste0(io$data_dir, "/rep", args$replicate, "/")

cat("Replicate:", args$replicate)
cat(", Cells:", args$ncells)
cat(", Features:", args$nfeatures)
cat(", Tolerance:",args$tolerance)

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

N_cpgs <- 15
sim_dt <- readRDS(paste0(io$data_dir, "data_bbreg_feat", args$nfeatures, "_cells", args$ncells, "_cpgs", N_cpgs, ".rds"))

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

###############
## Fit model ##
###############

cat("Doing inference...\n")

fit_bbreg <- infer_bbreg(Y = sim_dt$Y, stan_model_file = io$stan_model_file, X = sim_dt$sim_dt_A$X,
                           L = L_fit, use_mcmc = FALSE, m_wmu = m_wmu, s_wmu = s_wmu, s_mu = s_mu,
                           m_wgamma = m_wgamma, s_wgamma = s_wgamma, a_sgamma = a_sgamma, b_sgamma = b_sgamma,
                           rbf_c = rbf_c, use_eb = TRUE, init_using_eb = TRUE, iter = iter, tol_rel_obj = args$tolerance)

obj <- list(fit_bbreg = fit_bbreg, sim_dt = sim_dt)

#################
## Save output ##
#################

outfile <- paste0(io$outdir, "/bbreg_feat", args$nfeatures, "_cells", args$ncells, "_cpgs", N_cpgs, ".rds")
saveRDS(object = obj, file = outfile)
