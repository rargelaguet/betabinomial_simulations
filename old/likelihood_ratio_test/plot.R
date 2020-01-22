library(ggplot2)
library(ggpubr)
library(data.table)
library(purrr)

######################
## Define settings  ##
######################

## I/O ##
io <- list()
io$fitted.models <- "/Users/ricard/data/betabinomial/validation/fitted_models.rds"
io$ground.truth <- "/Users/ricard/data/betabinomial/validation/true_params.rds"
io$outdir <- "/Users/ricard/data/betabinomial/validation/pdf"

## Options ##
opts <- list()

# Number of trials
opts$ntrials <- 5

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

tmp <- readRDS(io$fitted.models)
fit.mcmc <- tmp$MCMC
fit.advi.meanfield <- tmp$VB.meanfield
fit.advi.fullrank <- tmp$VB.fullrank

################################
## Compare summary statistics ##
################################

# params <- c("mu","gamma", "delta", "rho")
params <- c("mu","gamma")

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
  .[,parameter_feature:=parameter] %>%
  # .[,vb_parameter:=paste(parameter,VB, collapse=" ")] %>%
  .[,feature:=unlist(regmatches(parameter, gregexpr("[[:digit:]]+", parameter)))] %>%
  .[,parameter:=gsub('[[:digit:]]+', '', parameter)]

features.to.plot <- as.character( 1:12 )
  
for (i in params) {
  for (j in unique(to.plot$VB)) {
  
    # posterior mean
    p1 <- ggscatter(to.plot[parameter==i & feature%in%features.to.plot & estimate=="mean" & VB==j], x="value", y="MCMC", 
                    facet="parameter_feature", add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=FALSE, fullrange=TRUE) +
      labs(x="Posterior mean (Variational)", y="Posterior mean (MCMC)") +
      # facet_wrap(~feature) +
      coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
      geom_abline(slope=1, intercept=0, linetype="solid")
    # p1
    
    
    # posterior sd
    p2 <- ggscatter(to.plot[parameter==i & feature%in%features.to.plot & estimate=="sd" & VB==j], x="value", y="MCMC", 
                    facet="parameter_feature", add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=FALSE, fullrange=TRUE) +
      labs(x="Posterior stdev (Variational)", y="Posterior stdev (MCMC)") +
      # facet_wrap(~) +
      coord_cartesian(xlim=c(0,0.2), ylim=c(0,0.2)) +
      geom_abline(slope=1, intercept=0, linetype="solid")
    # p2
    
    p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=2)
    
    
    pdf(sprintf("%s/mcmc_vs_vb_%s_%s.pdf",io$outdir,i,j), useDingbats = F, width=12, height=15)
    print(p)
    dev.off()
  }
  
}


##################################################
## Histograms comparing posterior distributions ##
##################################################

# dt.mcmc <- list()
# dt.advi_meanfield <- list()
# dt.advi_fullrank <- list()
# for (i in head(1:opts$ntrials)) {
#   
#   dt.mcmc[[i]] <- rstan::extract(fit.mcmc[[i]])[["theta"]] %>%
#     as.data.frame %>% tibble::rownames_to_column("iteration") %>% as.data.table %>% 
#     setnames(".","theta") %>%
#     melt(id.vars="iteration", variable.name="id") %>%
#     .[,inference:="MCMC"] %>% .[,trial:=i]
#   
#   dt.advi_meanfield[[i]] <- rstan::extract(fit.advi.meanfield[[i]])[["theta"]] %>%
#     as.data.frame %>% tibble::rownames_to_column("iteration") %>% as.data.table %>% 
#     setnames(".","theta") %>%
#     melt(id.vars="iteration", variable.name="id") %>%
#     .[,inference:="ADVI (mean-field)"] %>% .[,trial:=i]
#   
#   dt.advi_fullrank[[i]] <- rstan::extract(fit.advi.fullrank[[i]])[["theta"]] %>%
#     as.data.frame %>% tibble::rownames_to_column("iteration") %>% as.data.table %>% 
#     setnames(".","theta") %>%
#     melt(id.vars="iteration", variable.name="id") %>%
#     .[,inference:="ADVI (full-rank)"] %>% .[,trial:=i]
# }
# 
# dt.mcmc <- rbindlist(dt.mcmc)
# dt.advi_meanfield <- rbindlist(dt.advi_meanfield)
# dt.advi_fullrank <- rbindlist(dt.advi_fullrank)
# 
# to.plot <- rbindlist(list(dt.advi_meanfield, dt.advi_fullrank, dt.mcmc)) %>%
#   .[iteration>500] 
# 
# ggdensity(to.plot[id=="theta"], x="value", y="..density..", fill="inference", color="black", palette = "jco", 
#           facet="trial", scales="free_x") +
#   labs(x="", y="Density") +
#   theme(
#     legend.title = element_blank()
#   )


##########################################
## Compute difference with ground truth ##
##########################################

ground.truth <- readRDS(io$ground.truth)
for (i in 1:opts$ntrials) {
  ground.truth[[i]] <- ground.truth[[i]] %>% as.data.table %>% 
    setnames(c("mu","gamma")) %>%
    .[,feature:=1:.N] %>% 
    melt(id.vars="feature", variable.name="parameter", value.name="value") %>%
    .[,trial:=rep(i,.N)] %>%
    .[,parameter:=paste0(parameter,feature)]
}
ground.truth <- rbindlist(ground.truth)

to.plot <- merge(
  dt[estimate=="mean", c("parameter","inference","trial","value")],
  ground.truth[, c("parameter","trial","value")],
  by=c("parameter","trial")
) %>% setnames(c("parameter", "trial", "inference", "estimate", "true")) %>%
  .[,feature:=unlist(regmatches(parameter, gregexpr("[[:digit:]]+", parameter)))] %>%
  .[,parameter:=gsub('[[:digit:]]+', '', parameter)]


for (i in params) {
  to.plot2 <- to.plot[parameter==i]
  
  max.value = max(c(to.plot2$true,to.plot2$estimate))
  
  p <- ggscatter(to.plot2, x="estimate", y="true", 
                  facet="inference", add="reg.line", add.params = list(color="blue", fill="lightgray"), conf.int=FALSE, fullrange=TRUE) +
    coord_cartesian(xlim=c(0,max.value), ylim=c(0,max.value)) +
    geom_abline(slope=1, intercept=0, linetype="solid") +
    labs(x="Infered", y="True")
  
  pdf(sprintf("%s/bayesian_vs_truth_%s.pdf",io$outdir,i), useDingbats = F, width=11, height=7.5)
  print(p)
  dev.off()
  
}
