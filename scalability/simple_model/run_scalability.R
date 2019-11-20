library(data.table)
library(purrr)
library(rstan)
rstan_options(auto_write = TRUE)

simulate_binomial <- function(N, D, S) {
  theta <- runif(N*D, min=0, max=1)
  matrix(rbinom(N*D, S, theta), nrow=N, ncol=D)
}

######################
## Define settings  ##
######################

## I/O ##
io <- list()
io$basedir <- "/Users/ricard/betabinomial_simulations/scalability"
io$model <- paste0(io$basedir,"/model.R")
io$outdir <- paste0(io$basedir,"/out"); dir.create(io$outdir, showWarnings = F)

## Options ##
opts <- list()

# Define default dimensions
opts$default.N = 100
opts$default.D = 100
opts$default.S = 25

# Number of trials per setting
opts$ntrials <- 5

# Define ranges 
opts$N <- seq(10,100, by=10)
opts$D <- seq(10,100, by=10)

# Number of cores
opts$cores <- 1
options(mc.cores = opts$cores)

#######################
## Create stan model ##
#######################

# Load model
source(io$model)

st_model <- stan_model(model_code = beta_binomial, model_name="beta_binomial")

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
    
    Y <- simulate_binomial(N=n, D=opts$default.D, opts$default.S)
    
    # MCMC
    ptm <- proc.time()
    for (d in 1:opts$default.D) {
      data <- list(y=Y[,d], N=n, s=opts$default.S)
      sampling(st_model,  data = data, chains = 1)
    }
    mcmc[N==n & trial==i, time:=(proc.time()-ptm)["elapsed"]]

    # ADVI
    ptm <- proc.time()
    for (d in 1:opts$default.D) {
      data <- list(y=Y[,d], N=n, s=opts$default.S)
      vb(st_model,  data = data)
    }
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
for (D in opts$D) {
  for (i in 1:opts$ntrials) {
    
    Y <- simulate_binomial(N=opts$default.N, D=D, opts$default.S)
    
    # MCMC
    ptm <- proc.time()
    for (d in 1:D) {
      data <- list(y=Y[,d], N=opts$default.N, s=opts$default.S)
      sampling(st_model,  data = data, chains = 1)
    }
    mcmc[D==d & trial==i, time:=(proc.time()-ptm)["elapsed"]]
    
    # ADVI
    ptm <- proc.time()
    for (d in 1:D) {
      data <- list(y=Y[,d], N=opts$default.N, s=opts$default.S)
      vb(st_model,  data = data)
    }
    advi[D==d & trial==i, time:=(proc.time()-ptm)["elapsed"]]
    
  }
}
sink()

df <- rbind(mcmc,advi)
fwrite(df, paste0(io$outdir, "/D.txt.gz"), col.names=T, quote=F, sep="\t")

