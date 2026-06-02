# Dedicated script for generating all simulation plots

message("=========================================================")
message("Starting Simulation Plot Generation...")

source("code/functions/sim_results_functions.R")

# Detect whether to look at paper or user results
dir_path <- ifelse(USE_PAPER_RESULTS, "results/data/paper/simulations", "results/data/user/simulations")

# Look for AR data
ar_results_path <- file.path(dir_path, "AR", "final_results_AR.RData")
if (file.exists(ar_results_path)) {
  ar_env <- new.env()
  load(ar_results_path, envir = ar_env)
  final_results_AR <- ar_env$final_results_AR
  dgp_title_map_AR <- c(
    "AR_DGP_1_n100" = "AR DGP 1, T = 100", "AR_DGP_1_n200" = "AR DGP 1, T = 200",
    "AR_DGP_2_n100" = "AR DGP 2, T = 100", "AR_DGP_2_n200" = "AR DGP 2, T = 200",
    "AR_DGP_3_n100" = "AR DGP 3, T = 100", "AR_DGP_3_n200" = "AR DGP 3, T = 200"
  )
  plot_all_sim_results(
    final_results = final_results_AR, 
    dgp_title_map = dgp_title_map_AR, 
    plot_fun = generate_AR_barplot,
    save_dir = file.path(FIGURES_DIR, "simulations/AR")
  )
  message("   ✓ AR plots saved to figs dir")
}

# Look for Mixed data  
mixed_results_path <- file.path(dir_path, "Mixed", "final_results_Mixed.RData")
if (file.exists(mixed_results_path)) {
  mixed_env <- new.env()
  load(mixed_results_path, envir = mixed_env)
  final_results_Mixed <- mixed_env$final_results_Mixed
  dgp_title_map_Mixed <- c(
    "Mixed_DGP_1_n100" = "Mixed DGP 1, T = 100", "Mixed_DGP_1_n200" = "Mixed DGP 1, T = 200",
    "Mixed_DGP_2_n100" = "Mixed DGP 2, T = 100", "Mixed_DGP_2_n200" = "Mixed DGP 2, T = 200",
    "Mixed_DGP_3_n100" = "Mixed DGP 3, T = 100", "Mixed_DGP_3_n200" = "Mixed DGP 3, T = 200"
  )
  plot_all_sim_results(
    final_results = final_results_Mixed, 
    dgp_title_map = dgp_title_map_Mixed, 
    plot_fun = generate_Mixed_barplot,
    save_dir = file.path(FIGURES_DIR, "simulations/Mixed")
  )
  message("   ✓ Mixed plots saved to figs dir")
}

message("Simulation Plot Generation Complete!")
message("=========================================================")
