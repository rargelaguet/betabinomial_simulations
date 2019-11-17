library(ggplot2)
library(ggpubr)
library(data.table)
library(purrr)

######################
## Define settings  ##
######################

## I/O ##
io <- list()
io$file <- "/Users/ricard/stan/mcmc_vs_variational/betabinomial/out/fitted_models.rds"
io$outdir <- "/Users/ricard/stan/mcmc_vs_variational/betabinomial/pdf"

## Options ##
opts <- list()

# Number of trials
opts$ntrials <- 25

############################
## User-defined functions ##
############################

data_summary <- function(data){
  mean <- colMeans(data)
  sd <- apply(data, 2, sd)
  return(data.table("pred_mean" = mean, "pred_sd" = sd))
}

extract_summary_statistics <- function(fit, params) {
  lapply(params, function(j) {
    data_summary(as.data.frame(rstan::extract(fit)[[j]])) %>% 
      as.data.table %>% setnames(c("mean","sd")) %>%
      .[,parameter:=paste0(j,1:.N)]
  }) %>% rbindlist
}

###############################
## Load pre-computed results ##
###############################

tmp <- readRDS(io$file)
fit.mcmc <- tmp$MCMC
fit.advi.meanfield <- tmp$VB.meanfield
fit.advi.fullrank <- tmp$VB.fullrank

################################
## Compare summary statistics ##
################################

params <- c("theta")

# MCMC: Extract summary statistics for the posterior distributions
dt.mcmc <- lapply(1:opts$ntrials, function(i) {
	extract_summary_statistics(fit.mcmc[[i]], params) %>%
   .[,inference:="MCMC"] %>%
   .[,trial:=i]
  }) %>% rbindlist


# ADVI mean-field: Extract summary statistics for the posterior distributions
dt.advi_meanfield <- lapply(1:opts$ntrials, function(i) {
  extract_summary_statistics(fit.advi.meanfield[[i]], params) %>%
    .[,inference:="(mean-field)"] %>%
    .[,trial:=i]
}) %>% rbindlist

# ADVI full-rank: Extract summary statistics for the posterior distributions
dt.advi_fullrank <- lapply(1:opts$ntrials, function(i) {
  extract_summary_statistics(fit.advi.fullrank[[i]], params) %>%
    .[,inference:="(full-rank)"] %>%
    .[,trial:=i]
}) %>% rbindlist

dt <- rbindlist(list(dt.advi_meanfield, dt.advi_fullrank, dt.mcmc)) %>%
  melt(id.vars=c("parameter","inference","trial"), variable.name="estimate")

################################################
## Scatter plots comparing summary statistics ##
################################################

to.plot <- dt %>% dcast(parameter+trial+estimate ~ inference) %>%
  melt(id.vars=c("parameter","trial","estimate","MCMC"), variable.name="VB") %>%
  .[,vb_parameter:=paste(parameter,VB, sep=" ")]

# posterior mean
p1 <- ggscatter(to.plot[estimate=="mean"], x="value", y="MCMC", facet="vb_parameter") +
  labs(x="Posterior mean (Variational)", y="Posterior mean (MCMC)") +
  coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
  geom_abline(slope=1, intercept=0, linetype="solid")
# p1


# posterior sd
p2 <- ggscatter(to.plot[estimate=="sd"], x="value", y="MCMC", facet="vb_parameter") +
  labs(x="Posterior sd (Variational)", y="Posterior sd (MCMC)") +
  # coord_cartesian(xlim=c(0,0.025), ylim=c(0,0.025)) +
  geom_abline(slope=1, intercept=0, linetype="solid")
# p2

cowplot::plot_grid(plotlist=list(p1,p2), nrow=2)


##################################################
## Histograms comparing posterior distributions ##
##################################################

dt.mcmc <- list()
dt.advi_meanfield <- list()
dt.advi_fullrank <- list()
for (i in head(1:opts$ntrials)) {
  
  dt.mcmc[[i]] <- rstan::extract(fit.mcmc[[i]])[["theta"]] %>%
    as.data.frame %>% tibble::rownames_to_column("iteration") %>% as.data.table %>% 
    setnames(".","theta") %>%
    melt(id.vars="iteration", variable.name="id") %>%
    .[,inference:="MCMC"] %>% .[,trial:=i]
  
  dt.advi_meanfield[[i]] <- rstan::extract(fit.advi.meanfield[[i]])[["theta"]] %>%
    as.data.frame %>% tibble::rownames_to_column("iteration") %>% as.data.table %>% 
    setnames(".","theta") %>%
    melt(id.vars="iteration", variable.name="id") %>%
    .[,inference:="ADVI (mean-field)"] %>% .[,trial:=i]
  
  dt.advi_fullrank[[i]] <- rstan::extract(fit.advi.fullrank[[i]])[["theta"]] %>%
    as.data.frame %>% tibble::rownames_to_column("iteration") %>% as.data.table %>% 
    setnames(".","theta") %>%
    melt(id.vars="iteration", variable.name="id") %>%
    .[,inference:="ADVI (full-rank)"] %>% .[,trial:=i]
}

dt.mcmc <- rbindlist(dt.mcmc)
dt.advi_meanfield <- rbindlist(dt.advi_meanfield)
dt.advi_fullrank <- rbindlist(dt.advi_fullrank)

to.plot <- rbindlist(list(dt.advi_meanfield, dt.advi_fullrank, dt.mcmc)) %>%
  .[iteration>500] 

ggdensity(to.plot[id=="theta"], x="value", y="..density..", fill="inference", color="black", palette = "jco", 
          facet="trial", scales="free_x") +
  labs(x="", y="Density") +
  theme(
    legend.title = element_blank()
  )
