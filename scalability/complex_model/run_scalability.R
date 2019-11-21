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
  io$basedir <- "/Users/ricard/betabinomial_simulations/scalability/complex_model"
  io$outdir <- "/Users/ricard/data/betabinomial/scalability"
  io$model <- "/Users/ricard/Ecker_2017/variability/bayesian_bb/lib/bb_hierarchical.stan"
} else {
  io$basedir <- "/homes/ricard/betabinomial_simulations/scalability/complex_model"
  io$model <- "/homes/ricard/Ecker_2017/variability/bayesian_bb/lib/bb_hierarchical.stan"
  io$outdir <- "/hps/nobackup2/stegle/users/ricard/betabinomial/scalability"
}

## Options ##
opts <- list()

# Number of trials per setting
opts$ntrials <- 10

# Define ranges 
opts$N <- seq(10,150, by=10)
opts$D <- seq(10,150, by=10)

# Number of cores
opts$cores <- 1
options(mc.cores = opts$cores)

# Simulation settings
opts$N_feat = 50     # Default number of features
opts$N_cells = 50    # Default number of cells
opts$N_cpgs = 15     # Number of CpGs per genomic feature
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

# Maximum number of iterations
# opts$iterations <- 3000

# Convergence criteria for VB
# opts$tol_rel_obj <- 0.001

# Load data-generating function
source(paste0(io$basedir,"/../../load_data.R"))

#######################
## Create stan model ##
#######################

st_model <- stan_model(file = io$model, model_name="beta_binomial")

#########################
## Test influence of N ##
#########################

advi <- expand.grid(opts$N, 1:opts$ntrials) %>% as.data.table %>%
  setnames(c("N","trial")) %>%.[,time:=-Inf] %>% .[,inference:="ADVI"]

mcmc <- expand.grid(opts$N, 1:opts$ntrials) %>% as.data.table %>%
  setnames(c("N","trial")) %>%.[,time:=-Inf] %>% .[,inference:="MCMC"]

sink(tempfile())
for (n in opts$N) {
  for (i in 1:opts$ntrials) {

    # Simulate data
    met <- simulate_met(
      N_feat = opts$N_feat,
      N_cells = n,
      N_cpgs = opts$N_cpgs,
      cells_range = opts$cells_range,
      cpgs_range = opts$cpgs_range,
      mean_met_range = opts$mean_met_range,
      disp_met_range = opts$disp_met_range,
      seed = i
    )[[1]]

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

    # Fit MCMC
    ptm <- proc.time()
    fit.mcmc <- sampling(
      object = st_model, data = dat, chains = 1,
      pars = c("mu","gamma", "delta", "rho"),
      # iter = opts$iterations,
      cores = opts$cores
    )
    mcmc[N==n & trial==i, time:=(proc.time()-ptm)["elapsed"]]
  
    # Fit ADVI
    ptm <- proc.time()
    fit.vb <- vb(
      object = st_model, data = dat,
      pars = c("mu","gamma", "delta", "rho"),
      init = 'random',
      algorithm = "meanfield"
      # iter = opts$iterations,
      # tol_rel_obj = opts$tol_rel_obj
    )
    advi[N==n & trial==i, time:=(proc.time()-ptm)["elapsed"]]
    
  }
}
sink()

df <- rbind(mcmc,advi)
fwrite(df, paste0(io$outdir, "/N.txt.gz"), col.names=T, quote=F, sep="\t")

#########################
## Test influence of D ##
#########################

advi <- expand.grid(opts$D, 1:opts$ntrials) %>% as.data.table %>% 
  setnames(c("D","trial")) %>%.[,time:=-Inf] %>% .[,inference:="ADVI"]

mcmc <- expand.grid(opts$D, 1:opts$ntrials) %>% as.data.table %>% 
  setnames(c("D","trial")) %>%.[,time:=-Inf] %>% .[,inference:="MCMC"]

sink(tempfile())
for (d in opts$D) {
  for (i in 1:opts$ntrials) {
    
    # Simulate data
    met <- simulate_met(
      N_feat = d,
      N_cells = opts$N_cells,
      N_cpgs = opts$N_cpgs,
      cells_range = opts$cells_range,
      cpgs_range = opts$cpgs_range,
      mean_met_range = opts$mean_met_range,
      disp_met_range = opts$disp_met_range,
      seed = i
    )[[1]]

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


    # Fit MCMC
    ptm <- proc.time()
    fit.mcmc <- sampling(
      object = st_model, data = dat, chains = 1,
      pars = c("mu","gamma", "delta", "rho"),
      # iter = opts$iterations, 
      cores = opts$cores
    )
    mcmc[D==d & trial==i, time:=(proc.time()-ptm)["elapsed"]]
  
    # Fit ADVI
    ptm <- proc.time()
    fit.vb <- vb(
      object = st_model, data = dat,
      pars = c("mu","gamma", "delta", "rho"),
      init = 'random',
      algorithm = "meanfield"
      # iter = opts$iterations,
      # tol_rel_obj = opts$tol_rel_obj
    )
    advi[D==d & trial==i, time:=(proc.time()-ptm)["elapsed"]]
    
  }
}
sink()

df <- rbind(mcmc,advi)
fwrite(df, paste0(io$outdir, "/D.txt.gz"), col.names=T, quote=F, sep="\t")
