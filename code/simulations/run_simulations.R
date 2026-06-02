message("=========================================================")
message("Starting Simulation Pipeline...")

# 1. Define the Full Universe of DGPs
all_DGPs_AR <- c(
  "AR_DGP_1_n100", "AR_DGP_1_n200",
  "AR_DGP_2_n100", "AR_DGP_2_n200",
  "AR_DGP_3_n100", "AR_DGP_3_n200"
)

all_DGPs_Mixed <- c(
  "Mixed_DGP_1_n100", "Mixed_DGP_1_n200",
  "Mixed_DGP_2_n100", "Mixed_DGP_2_n200",
  "Mixed_DGP_3_n100", "Mixed_DGP_3_n200"
)

# 2. Map Plot Titles
dgp_title_map_AR <- c(
  "AR_DGP_1_n100" = "AR DGP 1, T = 100", "AR_DGP_1_n200" = "AR DGP 1, T = 200",
  "AR_DGP_2_n100" = "AR DGP 2, T = 100", "AR_DGP_2_n200" = "AR DGP 2, T = 200",
  "AR_DGP_3_n100" = "AR DGP 3, T = 100", "AR_DGP_3_n200" = "AR DGP 3, T = 200"
)

dgp_title_map_Mixed <- c(
  "Mixed_DGP_1_n100" = "Mixed DGP 1, T = 100", "Mixed_DGP_1_n200" = "Mixed DGP 1, T = 200",
  "Mixed_DGP_2_n100" = "Mixed DGP 2, T = 100", "Mixed_DGP_2_n200" = "Mixed DGP 2, T = 200",
  "Mixed_DGP_3_n100" = "Mixed DGP 3, T = 100", "Mixed_DGP_3_n200" = "Mixed DGP 3, T = 200"
)

# 3. Handle Paper Replication Override
if (REPLICATE_PAPER_RESULTS) {
  message("-> REPLICATE_PAPER_RESULTS = TRUE: overriding DGP selections to run ALL.")
  selected_dgps_AR    <- all_DGPs_AR
  selected_dgps_Mixed <- all_DGPs_Mixed
  RESULTS_DIR <- "results/data/user"
} else {
  # 4. Read User Input from user_run_file.R
  # Map the global uppercase variables to the lowercase ones expected by AR_run/Mixed_run
  selected_dgps_AR    <- SELECTED_DGPS_AR
  selected_dgps_Mixed <- SELECTED_DGPS_MIXED
}

# Map global M_SIM to the internal 'M' variable expected by AR_run/Mixed_run
M <- M_SIM 

# ==============================================================================
# AR SIMULATIONS
# ==============================================================================
if (RUN_SIM_AR) {
  if (length(selected_dgps_AR) > 0) {
    message(sprintf("-> Running AR simulations for: %s", paste(selected_dgps_AR, collapse = ", ")))
    source(file = "code/simulations/AR_run.R")
    
    dir.create(file.path(RESULTS_DIR, "simulations/AR"), recursive = TRUE, showWarnings = FALSE)
    save(final_results_AR_pre,file = file.path(RESULTS_DIR, "simulations/AR/final_results_AR_pre.RData"))
    save(final_results_AR, file = file.path(RESULTS_DIR, "simulations/AR/final_results_AR.RData"))
    
    message("   -> Generating AR plots via dedicated script...")
    source("code/simulations/run_simulation_plots.R")
  }
}

# ==============================================================================
# MIXED-FREQUENCY SIMULATIONS
# ==============================================================================
if (RUN_SIM_MIXED) {
  if (length(selected_dgps_Mixed) > 0) {
    message(sprintf("-> Running Mixed simulations for: %s", paste(selected_dgps_Mixed, collapse = ", ")))
    source(file = "code/simulations/Mixed_run.R")
    
    dir.create(file.path(RESULTS_DIR, "simulations/Mixed"), recursive = TRUE, showWarnings = FALSE)
    save(final_results_Mixed_pre, file = file.path(RESULTS_DIR, "simulations/Mixed/final_results_Mixed_pre.RData"))
    save(final_results_Mixed, file = file.path(RESULTS_DIR, "simulations/Mixed/final_results_Mixed.RData"))
    
    message("   -> Generating Mixed plots via dedicated script...")
    source("code/simulations/run_simulation_plots.R")
  }
}

message("Simulation Pipeline Complete!")
message("=========================================================")
