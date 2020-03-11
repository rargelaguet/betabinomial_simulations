
#########
## I/O ##
#########

io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  # io$basedir <- "/Users/ricard/data/betabinomial"
  io$script <- "/Users/ricard/betabinomial_simulations/variational_convergence/fit_models.R"
} else if (grepl("ebi",Sys.info()['nodename'])) {
  io$tmpdir <- "/hps/nobackup2/research/stegle/users/ricard/betabinomial/variational_tolerance/tmp"
  io$script <- "/homes/ricard/betabinomial_simulations/variational_convergence/fit_models.R"
} else {
  stop("Computer not recognised")
}
# io$tmpdir <- paste0(io$basedir,"/tmp")
# dir.create(io$outdir, showWarnings = F)
dir.create(io$tmpdir, showWarnings = F)

#############
## Options ##
#############

opts <- list()

# Tolerance for convergence for the VB algorithm
# opts$tolerance_thresholds <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)
opts$tolerance_thresholds <- c(5e-5)
# opts$tolerance_thresholds <- c(0.1)

# Testing mode (subsetting the number of features)
opts$test <- FALSE

# Number of replicates
opts$replicates <- 1

# Number of cells
opts$ncells <- seq(100,500,by=100)

# Number of features
opts$nfeatures <- seq(100,500,by=100)

##############
## Run jobs ##
##############

for (i in opts$ncells) {
  for (j in opts$nfeatures) {
    for (r in opts$replicates) {
      for (k in opts$tolerance_thresholds) {

        # lsf <- ""
        lsf <- sprintf("bsub -M 4096 -n 1 -q research-rh74 -o %s/%s_%s_%s_%s.txt", io$tmpdir, i,j,r,k)

        # Run
        cmd <- sprintf("%s Rscript %s --ncells %s --nfeatures %s --replicate %s --tolerance %s --seed %s", 
          lsf, io$script, i, j, r, k, r)
        if (opts$test) cmd <- paste0(cmd," --test")
        system(cmd)
      }
    }
  }
}
