source("code/functions/macro_functions.R")

# 1. Define the Configuration Grid based on mode
if (REPLICATE_PAPER_RESULTS) {
  RESULTS_DIR <- "results/data/user"
  configs <- expand.grid(
    dataset   = c("nowcast", "forecast"),
    reduction = c("full", "reduced"),
    window    = c(66, 105),
    h         = 0,
    up        = 1,
    b2_final  = FALSE,
    stringsAsFactors = FALSE
  ) } else {
  configs <- expand.grid(
    dataset   = MACRO_DATASET,
    reduction = MACRO_REDUCTION,
    window    = MACRO_WINDOW,
    stringsAsFactors = FALSE
  )
  # Add fixed columns that don't vary across combinations
  configs$h <- 0
  configs$up <- 1
  configs$b2_final <- FALSE
  }


# } else {
#   configs <- data.frame(
#     dataset   = MACRO_DATASET,
#     reduction = MACRO_REDUCTION,
#     window    = MACRO_WINDOW,
#     h         = 0,
#     up        = 1,
#     b2_final  = FALSE,
#     stringsAsFactors = FALSE
#   )
# }


# 2. Execution
out_folder <- file.path(RESULTS_DIR, "applications/macro")
if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

message(sprintf("Starting Macro batch processing (%d configurations)...", nrow(configs)))

# Initialize list to store summaries
all_summaries <- list()

for (i in seq_len(nrow(configs))) {
  
  cfg <- configs[i, ]
  
  file_prefix <- ifelse(cfg$dataset == "nowcast", "nowcast", "forecast")
  label <- paste0(file_prefix, "_", cfg$reduction, "_w", cfg$window)
  out_filename <- file.path(out_folder, paste0(label, ".RData"))
  
  # --- THE RESUME-ABLE CHECK ---
  if (file.exists(out_filename)) {
    message(sprintf("-> SKIPPING Macro App: %s (Already computed)", label))
    env <- new.env()
    load(out_filename, envir = env)
    all_summaries[[i]] <- env$summary_row
    next
  }
  
  message(sprintf("-> RUNNING Macro App: %s", label))
  
  tryCatch({
    # Load Data & Build Trees
    dat <- get_Xy(cfg$dataset, cfg$reduction)
    groups <- build_groups(cfg$dataset, cfg$reduction)
    
    # Run Nowcast/Forecast
    res <- macro_run(X = dat$X, y = dat$y, window_size = cfg$window, h = cfg$h, up = cfg$up,
                     P_list = groups$P_list, Ki_list = groups$Ki_list, rho = 1,
                     K_in = 1000, stop_thresh = 1e-5, ic_type = "BIC", thresh = 0.3,
                     nlambda = 10, b2_final = cfg$b2_final, dataset = cfg$dataset)
    
    # Extract Metrics
    summary_row <- data.frame(
      dataset = cfg$dataset, reduction = cfg$reduction, window = cfg$window, h = cfg$h,
      mse_star_post = res$penalized$evaluation$mse_star_post,
      mse_star_simple = res$penalized$evaluation$mse_star_simple,
      mse_midasml = res$penalized$evaluation$mse_midasml,
      mse_rw = res$random_walk$evaluation$mse,
      mse_ar1 = res$ar1$evaluation$mse,
      stringsAsFactors = FALSE
    )
    
    # Save individual result safely
    save(res, summary_row, file = out_filename)
    
    # Add to our master list
    all_summaries[[i]] <- summary_row
    
  }, error = function(e) {
    message(sprintf("Task %d failed! ERROR: %s", i, e$message))
  })
}

# Save Master Summary Table
results_summary <- do.call(rbind, all_summaries)
if (!is.null(results_summary)) {
  write.csv(results_summary, file.path(out_folder, "macro_summary_table.csv"), row.names = FALSE)
}
message("Macro Application batch processing completed!")