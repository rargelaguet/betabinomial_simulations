library(data.table)
library(purrr)
library(rstan)
rstan_options(auto_write = TRUE)

######################
## Define settings  ##
######################

## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/betabinomial_simulations/mcmc_vs_variational/complex_model"
  io$model <- "/Users/ricard/Ecker_2017/variability/bayesian_bb/lib/bb_hierarchical.stan"
  io$outdir <- "/Users/ricard/data/betabinomial/validation"
} else {
  io$basedir <- "/homes/ricard/betabinomial_simulations/mcmc_vs_variational/complex_model"
  io$model <- "/homes/ricard/Ecker_2017/variability/bayesian_bb/lib/bb_hierarchical.stan"
  io$outdir <- "/Users/ricard/data/betabinomial/validation"
}

## Options ##
opts <- list()

# Number of trials
opts$ntrials <- 25

# Maximum number of iterations
opts$iterations <- 4000

# Convergence criteria for VB
opts$tol_rel_obj <- 0.001

# Number of cores
opts$cores <- 1
options(mc.cores = opts$cores)

# Simulation settings
opts$N_feat = 50
opts$N_cells = 50
opts$N_cpgs = 15
opts$cells_range = c(0.5, 1.0)       # min and max from runif
opts$cpgs_range = c(1.0, 1.0)        # min and max from runif
opts$mean_met_range = c(0.05, 0.95)  # min and max from runif
opts$disp_met_range = c(2, 8)        # shape1 and shape2 from beta

# Initialisations for the hyperparameters
opts$a0_mu = 1
opts$b0_mu = 1
opts$a0_delta = 1
opts$b0_delta = 2
opts$a0_rho = 20
opts$b0_rho = 100

###################################
## Load data-generating function ##
###################################

source(paste0(io$basedir,"/load_data.R"))

#################
## Build model ##
#################

st_model <- stan_model(file = io$model, model_name="beta_binomial")

#########################
## Fit Bayesian models ##
#########################

fit.mcmc <- list()
fit.vb.meanfield <- list()
fit.vb.fullrank <- list()

true_params <- list()

for (i in 1:opts$ntrials) {
  
  # Simulate data
  tmp <- simulate_met(
    N_feat = opts$N_feat,
    N_cells = opts$N_cells,
    N_cpgs = opts$N_cpgs,
    cells_range = opts$cells_range,
    cpgs_range = opts$cpgs_range,
    mean_met_range = opts$mean_met_range,
    disp_met_range = opts$disp_met_range,
    seed = i
  )

  met <- tmp[[1]]
  true_params[[i]] <- tmp[[2]]
  
  # Create data object for Stan
  n_obs_cells <- met[,.N, by = c("Feature")]$N
  
  dat <- list(
    N = nrow(met), 
    J = length(n_obs_cells), 
    y = met[ ,met_reads],
    n = met[ ,total_reads], 
    s = n_obs_cells, 
    a0_mu = opts$a0_mu, 
    b0_mu = opts$b0_mu,
    a0_delta = opts$a0_delta, 
    b0_delta = opts$b0_delta, 
    a0_rho = opts$a0_rho, 
    b0_rho = opts$b0_rho
  )
  
  # Perform inference using HMC sampling
  fit.mcmc[[i]] <- sampling(object = st_model, data = dat, chains = 1,
                            pars = c("mu","gamma", "delta", "rho"),
                            iter = opts$iterations, cores = opts$cores)
  
  # Perform inference using mean-field ADVI
  fit.vb.meanfield[[i]] <- vb(object = st_model, data = dat,
                              pars = c("mu","gamma", "delta", "rho"),
                              init = 'random',
                              algorithm = "meanfield",
                              iter = opts$iterations,
                              tol_rel_obj = opts$tol_rel_obj
  )

  # Perform inference using full-rank ADVI
  fit.vb.fullrank[[i]] <-  vb(object = st_model, data = dat,
                              pars = c("mu","gamma", "delta", "rho"),
                              init = 'random',
                              algorithm = "fullrank",
                              iter = opts$iterations,
                              tol_rel_obj = opts$tol_rel_obj)
  
}

##################
## Save results ##
##################

saveRDS(list("MCMC"=fit.mcmc, "VB.meanfield"=fit.vb.meanfield, "VB.fullrank"=fit.vb.fullrank), paste0(io$outdir,"/fitted_models.rds"))
saveRDS(true_params, paste0(io$outdir,"/true_params.rds"))
