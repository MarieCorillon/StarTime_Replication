# ========================================================
# Main script for executing AR model simulations
# ========================================================

#' Main Simulation Runner for AR Models
#' 
#' @description
#' Orchestrates the complete simulation pipeline for autoregressive models,
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
all_dgp_files_AR <- c(
  "code/simulations/DGPs/AR/AR_DGP_1_n100.R",
  "code/simulations/DGPs/AR/AR_DGP_1_n200.R",
  "code/simulations/DGPs/AR/AR_DGP_2_n100.R",
  "code/simulations/DGPs/AR/AR_DGP_2_n200.R",
  "code/simulations/DGPs/AR/AR_DGP_3_n100.R",
  "code/simulations/DGPs/AR/AR_DGP_3_n200.R"
)

names(all_dgp_files_AR) <- c(
  "AR_DGP_1_n100", "AR_DGP_1_n200", "AR_DGP_2_n100",
  "AR_DGP_2_n200", "AR_DGP_3_n100", "AR_DGP_3_n200"
)

dgp_files_AR <- all_dgp_files_AR[selected_dgps_AR]


registerDoRNG(20250610)

#' Execute AR Model Simulation for a Single DGP
#'
#' @param dgp_file Path to DGP configuration script
#' @param M Number of Monte Carlo repetitions (default = 500)
#' 
#' @return Data frame containing averaged performance metrics
#' 
#' @details
#' 1. Sources required functions in parallel workers
#' 2. Generates AR data according to DGP specifications
#' 3. Runs estimation methods in parallel
#' 4. Computes and stores performance metrics
run_dgp_simulation_AR <- function(dgp_file, M) {
  results <- foreach(i = 1:M,
                     .packages = c("glmnet", "Matrix", "StarTime"),
                     .combine = c) %dopar% {
                       # Seed inside each worker
                       set.seed(i + 10000)
                       
                       source(file = "code/functions/sim_estimation.R")
                       source(file = "code/functions/sim_data_generation.R")
                       source(file = "code/functions/metrics.R")
                       source(dgp_file)
                       list(run_AR_simulation(i))
                     }
  
  all_metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  averaged_metrics <- calculate_average_metrics(all_metrics)
  all_betas <- lapply(results, `[[`, "beta_estimates")
  
  return(list(metrics = averaged_metrics, betas = all_betas)) 
}




# Run All DGPs
final_results_AR_pre <- lapply(dgp_files_AR, 
                               function(f) run_dgp_simulation_AR(f, M))

names(final_results_AR_pre) <- tools::file_path_sans_ext(basename(dgp_files_AR))


# Clean up parallel resources
if (USE_PARALLEL) {
  parallel::stopCluster(cluster)
}

methods_to_keep_AR <- c(
  "ss", "ss.1", "ss.2", "ss.3",     # StarTime variants
  "lasso.1", "lasso.3",             # Lasso IC variants
  "ols_basic", "ridge_ic"           # Baseline methods
)

method_labels_AR <- c(
  ss = "StarTime (Simple, Oracle)",
  ss.1 = "StarTime (Simple, IC)",
  ss.2 = "StarTime (Post, Oracle)",
  ss.3 = "StarTime (Post, IC)",
  lasso.1 = "Lasso (Simple, IC)",
  lasso.3 = "Lasso (Post, IC)",
  ols_basic = "OLS",
  ridge_ic = "Ridge (IC)"
)

metric_labels <- c(
  betadiff = "MSE",
  ari = "ARI",
  f = "F1"
)


final_results_AR <- lapply(final_results_AR_pre, function(res) {
  relabel_result_table(
    res$metrics,
    methods_to_keep = methods_to_keep_AR,
    method_labels   = method_labels_AR,
    metric_labels   = metric_labels
  )
})
