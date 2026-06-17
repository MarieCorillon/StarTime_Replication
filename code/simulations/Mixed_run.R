# ========================================================
# Main script for executing Mixed-Frequency model simulations
# ========================================================

#' Main Simulation Runner for MIDAS Models
#'
#' @description
#' Orchestrates the complete simulation pipeline for MIDAS models,
#' including parallel execution, result processing, and output generation.
#'
#' @details
#' This script:
#' 1. Configures parallel processing
#' 2. Defines Data Generating Processes (DGPs)
#' 3. Runs Monte Carlo simulations
#' 4. Generates performance metrics and visualizations
#' 5. Produces publication-ready tables

source(file = "code/functions/sim_results_functions.R")

# Configure Parallel Processing
if (USE_PARALLEL) {
  cluster <- parallel::makeCluster(NUM_CORES)
  doParallel::registerDoParallel(cluster)
} else {
  foreach::registerDoSEQ()
}

# Define DGP Files
all_dgp_files_Mixed <- c(
  "code/simulations/DGPs/Mixed/Mixed_DGP_1_n100.R",
  "code/simulations//DGPs/Mixed/Mixed_DGP_1_n200.R",
  "code/simulations//DGPs/Mixed/Mixed_DGP_2_n100.R",
  "code/simulations//DGPs/Mixed/Mixed_DGP_2_n200.R",
  "code/simulations//DGPs/Mixed/Mixed_DGP_3_n100.R",
  "code/simulations//DGPs/Mixed/Mixed_DGP_3_n200.R"
)


names(all_dgp_files_Mixed) <- c(
  "Mixed_DGP_1_n100", "Mixed_DGP_1_n200", "Mixed_DGP_2_n100",
  "Mixed_DGP_2_n200", "Mixed_DGP_3_n100", "Mixed_DGP_3_n200"
)


dgp_files_Mixed <- all_dgp_files_Mixed[selected_dgps_Mixed]


registerDoRNG(20250612)

#' Execute Mixed Model Simulation for a Single DGP
#'
#' @param dgp_file Path to DGP configuration script
#' @param M Number of Monte Carlo repetitions (default = 500)
#' @return Data frame containing averaged performance metrics
run_dgp_simulation_Mixed <- function(dgp_file, M = 500) {
  results <- foreach(i = 1:M,
                     .packages = c("Matrix", "midasr", "midasml", "StarTime"),
                     .combine = c) %dopar% {
                       # Seed inside each worker
                       set.seed(i + 20000)

                       source(file = "code/functions/sim_estimation.R")
                       source(file = "code/functions/sim_data_generation.R")
                       source(file = "code/functions/metrics.R")
                       source(dgp_file)
                       list(run_Mixed_simulation(i))
                     }
  
  # Aggregate metrics across simulations
  all_metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  averaged_metrics <- calculate_average_metrics(all_metrics)
  
  all_betas <- lapply(results, `[[`, "beta_estimates")
  
  return(list(metrics = averaged_metrics, betas = all_betas))
}

# Run All DGPs
final_results_Mixed_pre <- lapply(dgp_files_Mixed, 
                                  function(f) run_dgp_simulation_Mixed(f, M))
names(final_results_Mixed_pre) <- tools::file_path_sans_ext(basename(dgp_files_Mixed))

# Clean up parallel resources
if (USE_PARALLEL) {
  parallel::stopCluster(cluster)
}


methods_to_keep_Mixed <- c(
  "ss", "ss.1", "ss.2", "ss.3",  # StarTime variants
  "ols_basic",                   # OLS
  "midas_ml_ic",                 # MIDAS-ML (IC)
  "midas"                        # MIDAS (non-penalized)
)

method_labels_Mixed <- c(
  ss = "StarTime (Simple, Oracle)",
  ss.1 = "StarTime (Simple, IC)",
  ss.2 = "StarTime (Post, Oracle)",
  ss.3 = "StarTime (Post, IC)",
  ols_basic = "OLS",
  midas_ml_ic = "MIDAS-ML (IC)",
  midas = "MIDAS"
)

metric_labels <- c(
  betadiff = "MSE",
  ari = "ARI",
  f = "F1"
)


final_results_Mixed <- lapply(final_results_Mixed_pre, function(res) {
  relabel_result_table(
    res$metrics,
    methods_to_keep = methods_to_keep_Mixed,
    method_labels   = method_labels_Mixed,
    metric_labels   = metric_labels
  )
})
