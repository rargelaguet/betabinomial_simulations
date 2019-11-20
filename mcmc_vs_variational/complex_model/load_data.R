library(data.table)
library(purrr)

simulate_met <- function(N_feat = 100, N_cells = 50, N_cpgs = 15,
                         cells_range = c(0.4, 0.8), cpgs_range = c(0.4, 0.8),
                         mean_met_range = c(0.04, 0.96), disp_met_range = c(2, 8), seed = 1) {
  # Set random seed
  set.seed(seed)
  
  # Total number of cells
  cells_vec <- rbinom(N_feat, N_cells, runif(1, min = cells_range[1], max = cells_range[2]))
  
  # Total number of CpGs for each cell
  cpgs_list <- lapply(cells_vec, function(x) {
    n <- rbinom(x, N_cpgs, runif(1, min = cpgs_range[1], max = cpgs_range[2]))
    # If we have features with less than 3 CpGs, ...
    idx <- which(n < 3)
    # ... sample values from [3, 5].
    if (length(idx) > 0) {
      n[idx] <- sample(3:5, length(idx), replace = TRUE)
    }
    return(n)
  })
  
  # Generate number of methylated CpGs from Beta Binomial
  bb_true_params <- t(sapply(cells_vec, function(x) {
    return(c(runif(1, min = mean_met_range[1], max = mean_met_range[2]),
             rbeta(1, shape1 = disp_met_range[1], shape2 = disp_met_range[2])))
  }))
  
  # Sample data
  met_cpgs_list <- lapply(X = 1:N_feat, function(n)
    VGAM::rbetabinom(length(cpgs_list[[n]]), cpgs_list[[n]], prob = bb_true_params[n, 1],
                     rho = bb_true_params[n, 2]))
  # Prepare output data structure
  met <- data.table("Feature" = unlist(sapply(1:N_feat, function(n) rep(paste0("Feature_", n), length(cpgs_list[[n]])))),
                    "Cell" = unlist(sapply(1:N_feat, function(n) paste0("Cell_", 1:length(cpgs_list[[n]])))) )
  met <- met %>% .[, c("total_reads", "met_reads") := list(unlist(cpgs_list), unlist(met_cpgs_list)) ]

  return(list(met = met, bb_true_params = bb_true_params))
}