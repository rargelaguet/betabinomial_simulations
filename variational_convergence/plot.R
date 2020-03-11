suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))

#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/betabinomial/variational_tolerance"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$basedir <- "/hps/nobackup2/research/stegle/users/ricard/betabinomial/variational_tolerance/models"
} else {
  stop("Computer not recognised")
}
io$outdir <- paste0(io$basedir,"/pdf")
dir.create(io$outdir, showWarnings = F)

#############
## Options ##
#############

opts <- list()

# Tolerance for convergence for the VB algorithm
# opts$tolerance_thresholds <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
opts$tolerance_thresholds <- c(1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1)
# opts$tolerance_thresholds <- c(1e-1)

# Number of replicates
opts$replicates <- c(1)

# Number of cells
opts$ncells <- seq(100,500,by=100)

# Number of features
opts$nfeatures <- seq(100,500,by=100)

##############################
## Load pre-computed models ##
##############################

dt <- opts$ncells %>% map(function(i) {
        opts$nfeatures %>% map(function(j) {
          opts$replicates %>% map(function(r) {
            opts$tolerance_thresholds %>% map(function(k) {
			        fread(sprintf("%s/models/%s_%s_%s_%s.txt.gz",io$basedir, i, j, r, k))
            }) %>% rbindlist
          }) %>% rbindlist
        }) %>% rbindlist
      }) %>% rbindlist %>%
  melt(id.vars=c("ncells","nfeatures","replicate","tolerance","id","time"), variable.factor=FALSE) %>%
  .[,estimate:=strsplit(variable,"_") %>% map_chr(2)] %>%
  .[,variable:=strsplit(variable,"_") %>% map_chr(1)]

#######################
## Load ground truth ##
#######################

dt.truth <- opts$ncells %>% map(function(i) {
  opts$nfeatures %>% map(function(j) {
    opts$replicates %>% map(function(r) {
        readRDS(sprintf("%s/data/rep%s/data_bbreg_feat%s_cells%s_cpgs15.rds",io$basedir, r, j, i))[["theta_true"]] %>%
          tibble::rownames_to_column("id") %>%
          as.data.table %>% .[,c("ncells","nfeatures","replicate"):=list(i,j,r)]
    }) %>% rbindlist
  }) %>% rbindlist
}) %>% rbindlist %>%
  melt(id.vars=c("ncells","nfeatures","replicate","id"), 
       variable.factor=FALSE,
       value.name="truth")

# Merge estimates and ground truth
dt2 <- dt[estimate=="mean"] %>% .[,estimate:=NULL] %>%
  merge(dt.truth, by=c("ncells","nfeatures","replicate","id","variable")) %>%
  .[,tolerance:=factor(format(tolerance,scientific=T), levels = format(opts$tolerance_thresholds, scientific = T))] %>%
  .[,ncells:=as.factor(ncells)] %>%
  .[,error:=(value-truth)**2]

#####################################################
## Plot tolerance vs error for different variables ##
#####################################################

to.plot <- dt2 %>% copy %>%
  .[,nfeatures:=as.factor(sprintf("Number of features: %s", nfeatures))]

limits <- c("gamma"=0.03, "mu"=0.006)

for (i in unique(to.plot$variable)) {
  
  p <- ggboxplot(to.plot[variable==i], x="ncells", y="error", fill="tolerance", outlier.shape=NA) +
    coord_cartesian(ylim=c(0,limits[[i]])) +
    facet_wrap(~nfeatures, scales="fixed", nrow=1) +
    scale_fill_brewer( palette="YlOrRd" ) +
    labs(x="Number of cells", y="Mean squared error") +
    theme(
      axis.text = element_text(size=rel(0.8)),
      legend.direction = "horizontal"
    )
  
  pdf(sprintf("%s/%s_tolerance_vs_error.pdf",io$outdir,i), width=13, height=5) 
  print(p)
  dev.off()
}

########################################################
## Plot tolerance vs standard deviation of posteriors ##
########################################################

to.plot <- dt[estimate=="sd"] %>% .[,estimate:=NULL] %>%
  .[,tolerance:=factor(format(tolerance,scientific=T), levels = format(opts$tolerance_thresholds, scientific = T))]

limits <- c("mu"=0.05, "gamma"=0.12, "epsilon"=0.8)

for (i in unique(to.plot$variable)) {
  
  p <- ggboxplot(to.plot[variable==i], x="ncells", y="value", fill="tolerance", outlier.shape=NA) +
    coord_cartesian(ylim=c(0,limits[[i]])) +
    facet_wrap(~nfeatures, scales="fixed", nrow=1) +
    scale_fill_brewer( palette="YlOrRd" ) +
    labs(x="Number of cells", y="Standard deviation of posterior distribution") +
    theme(
      axis.text = element_text(size=rel(0.8)),
      legend.direction = "horizontal"
    )
  
  pdf(sprintf("%s/%s_tolerance_vs_sd.pdf",io$outdir,i), width=13, height=5)
  print(p)
  dev.off()
}

############################################
## Plot tolerance vs time for convergence ##
############################################

to.plot <- dt2 %>%
  .[,.(time=unique(time)),c("ncells","nfeatures","replicate","tolerance")] %>%
  .[,ncells_nfeatures:=as.factor(sprintf("Cells: %s, Features: %s", ncells, nfeatures))]

to.plot[,log2_time:=log2(time)]

max.time <- 3000
to.plot[time>max.time,time:=max.time]

p <- ggbarplot(to.plot, x="tolerance", y="time", fill="tolerance", width=1) +
  facet_wrap(~ncells_nfeatures, scales="fixed") +
  labs(x="Tolerance", y="Time (sec)") +
  scale_fill_brewer( palette="YlOrRd" ) +
  theme(
    axis.text.y = element_text(size=rel(0.8)),
    axis.text.x = element_text(size=rel(0.8), angle=40, hjust=1, vjust=1),
    legend.position = "right"
  )

pdf(sprintf("%s/tolerance_vs_time.pdf",io$outdir), width=12, height=12) 
print(p)
dev.off()


#########################
## Plot training curve ##
#########################

# TO-DO...