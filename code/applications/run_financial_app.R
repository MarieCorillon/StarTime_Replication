# 1. Load Data and Source Functions
Variances_10min <- suppressWarnings(suppressMessages(
  readr::read_csv("data/financial/Variances_10min.csv", show_col_types = FALSE, name_repair = "minimal")
))
source("code/functions/financial_functions.R")

# 2. Determine Parameter Grids Based on Replication Mode
if (REPLICATE_PAPER_RESULTS) {
  # These are the exact combinations used in your paper
  grid_lags    <- c(20, 40)
  grid_windows <- c(1000, 250, 125)
  grid_h       <- c(1,5,20)
  RESULTS_DIR <- "results/data/user"
} else {
  # These are pulled from Section 3B of user_run_file.R
  grid_lags    <- c(FIN_MAX_LAG)
  grid_windows <- c(FIN_WINDOW_SIZE)
  grid_h       <- c(FIN_H)
}

# Fixed algorithmic parameters
rho <- 1; K_in <- 1000; stop_thresh <- 1e-05; ic_type <- "BIC" 
thresh <- 1; nlambda <- 10; up <- 1

# 3. Configure Execution Backend
if (USE_PARALLEL) {
  cl <- parallel::makeCluster(NUM_CORES)
  doParallel::registerDoParallel(cl)
} else {
  foreach::registerDoSEQ() 
}

doRNG::registerDoRNG(1234) 

# 4. Execute the Loops
out_folder <- file.path(RESULTS_DIR, "applications/financial")
if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

for (current_window in grid_windows) {
  for (current_lag in grid_lags) {
    for (current_h in grid_h) {
      
      message(sprintf("-> Running Financial App: Lag=%d, Window=%d, h=%d", 
                      current_lag, current_window, current_h))
      
      # Dynamically set tree structures based on the current lag in the loop
      if (current_lag == 20) {
        P_list  <- list(c(20, 4, 1))
        Ki_list <- list(c(5, 4, 1))
      } else if (current_lag == 40) {
        P_list  <- list(c(40, 8, 2, 1))
        Ki_list <- list(c(5, 4, 2, 1))
      }
      
      # Re-export variables to cluster for each new loop iteration
      if (USE_PARALLEL) {
        parallel::clusterExport(cl, varlist = c(
          "P_list", "Ki_list", "current_window", "current_h", "up",
          "rho", "K_in", "stop_thresh", "ic_type", "thresh", "nlambda",
          "prepare_window_data", "startime_rolling_predict", "ols_rolling_predict", 
          "random_walk_rolling_predict", "qlike" 
        ), envir = environment())
      }
      
      # Run the core function
      results <- parallel_run_all_stocks(
        variances_df = Variances_10min, 
        max_lag      = P_list[[1]][1], 
        window_size  = current_window, 
        h            = current_h, 
        up           = up,
        P_list       = P_list, 
        Ki_list      = Ki_list, 
        rho          = rho,
        K_in         = K_in, 
        stop_thresh  = stop_thresh, 
        ic_type      = ic_type,
        thresh       = thresh, 
        nlambda      = nlambda
      )
      
      # Dynamically Save Results for this specific combination
      out_filename <- file.path(out_folder, sprintf("Financial_lag%d_h%d_window%d.RData", 
                                                    current_lag, current_h, current_window))
      save(results, file = out_filename)
    }
  }
}

# 5. Clean up
if (USE_PARALLEL) {
  parallel::stopCluster(cl)
}
message("Financial application batch processing completed!")