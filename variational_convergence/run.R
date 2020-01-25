
## I/O ##
io <- list()
if (grepl("ricard",Sys.info()['nodename'])) {
  io$basedir <- "/Users/ricard/data/betabinomial_simulations"
  io$script <- "/Users/ricard/betabinomial_simulations/variational_convergence/fit_models.R"
} else if (grepl("yoda",Sys.info()['nodename'])) {
  stop()
  # io$basedir <- "/hps/nobackup2/stegle/users/ricard/betabinomial_simulations"
  # io$script <- "/homes/ricard/betabinomial_simulations/variational_convergence/fit_models.R"
} else {
  stop("Computer not recognised")
}
# io$tmpdir <- paste0(io$basedir,"/tmp")
# dir.create(io$outdir, showWarnings = F)
# dir.create(io$tmpdir, showWarnings = F)

## Options ##
opts <- list()

# Tolerance for convergence for the VB algorithm
opts$tolerance_thresholds <- c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1)

# Testing mode (subsetting the number of features)
opts$test <- FALSE

# Number of replicates
opts$replicates <- c(1)

# Number of cells
opts$ncells <- c(100)

# Number of features
opts$nfeatures <- c(500)

for (i in opts$ncells) {
  for (j in opts$nfeatures) {
    for (r in opts$replicates) {
      for (k in opts$tolerance_thresholds) {


        lsf <- ""
        # lsf <- sprintf("bsub -M 5120 -n 1 -q standard -o %s/%s.txt", io$tmpdir, i)

        # Run
        seed <- r
        cmd <- sprintf("%s Rscript %s --ncells %s --nfeatures %s --replicate %s --tolerance %s --seed %s", 
          lsf, io$script, i, j, r, k, seed)
        if (opts$test) cmd <- paste0(cmd," --test")
        system(cmd)
      }
    }
  }
}
