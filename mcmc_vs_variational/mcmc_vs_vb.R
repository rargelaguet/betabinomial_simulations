library(data.table)
library(purrr)
library(rstan)
rstan_options(auto_write = TRUE)

######################
## Define settings  ##
######################

## I/O ##
io <- list()
io$basedir <- "/Users/ricard/stan/mcmc_vs_variational/betabinomial"
io$model <- paste0(io$basedir,"/model.R")
io$outdir <- paste0(io$basedir,"/out"); dir.create(io$outdir)

## Options ##
opts <- list()

# Number of trials
opts$ntrials <- 25

# Number of cores
opts$cores <- 2
options(mc.cores = opts$cores)

###############
## Load data ##
###############

# Load model
source(io$model)

# Load data set
source(paste0(io$basedir,"/load_data.R"))

# Utils
source("/Users/ricard/stan/mcmc_vs_variational/utils.R")

#########################
## Fit Bayesian models ##
#########################

# Create stan model
st_model <- stan_model(model_code = beta_binomial, model_name="beta_binomial")

fit.mcmc <- list()
fit.vb.meanfield <- list()
fit.vb.fullrank <- list()
for (i in 1:opts$ntrials) {
  
  # Create data object for Stan
  data <- list(y=Y[[i]], N=N, s=s)
  
  # Perform inference using HMC sampling
  fit.mcmc[[i]] <- sampling(st_model,  data = data, chains = 1, iter=3000)
  
  # Perform inference using mean-fieldADVI
  fit.vb.meanfield[[i]] <- vb(st_model,  data = data, algorithm="meanfield", tol_rel_obj=0.001)

  # Perform inference using full-rank ADVI
  fit.vb.fullrank[[i]] <- vb(st_model,  data = data, algorithm="fullrank", tol_rel_obj=0.001)
  
}

##################
## Save results ##
##################

io$outfile <- paste0(io$outdir,"/fitted_models.rds")
saveRDS(list("MCMC"=fit.mcmc, "VB.meanfield"=fit.vb.meanfield, "VB.fullrank"=fit.vb.fullrank), io$outfile)
