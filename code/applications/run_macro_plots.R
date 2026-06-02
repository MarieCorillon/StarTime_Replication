source("code/functions/macro_plot_functions.R")

message("=========================================================")
message("Starting Macro Plot Generation...")

in_folder  <- file.path(RESULTS_DIR, "applications/macro")
out_folder <- file.path(FIGURES_DIR, "applications/macro")
dates_path <- "data/macro/dates_y.RData"

if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

if (!file.exists(dates_path)) stop("Dates file missing: ", dates_path)
cat("Loading dates from:", dates_path, "\n")
date_env <- new.env()
load(dates_path, envir = date_env)

if ("dates_y" %in% ls(date_env)) {
  full_dates_vector <- date_env$dates_y
} else {
  obj_names <- ls(date_env)
  full_dates_vector <- date_env[[obj_names[1]]]
  warning(paste("Variable 'dates_y' not found. Using:", obj_names[1]))
}

# 1. Determine Configurations
if (REPLICATE_PAPER_RESULTS) {
  plot_configs <- expand.grid(
    dataset   = c("nowcast", "forecast"),
    reduction = c("full", "reduced"),
    window    = c(66, 105),
    stringsAsFactors = FALSE
  )
  plot_configs$h <- 0 
} else {
  plot_configs <- expand.grid(
    dataset   = MACRO_DATASET,
    reduction = MACRO_REDUCTION,
    window    = MACRO_WINDOW,
    h         = 0,
    stringsAsFactors = FALSE
  )
}

# 2. Execute Plotting
for (i in 1:nrow(plot_configs)) {
  ds  <- plot_configs$dataset[i]
  red <- plot_configs$reduction[i]
  w   <- plot_configs$window[i]
  h   <- plot_configs$h[i] 
  
  file_prefix <- ifelse(ds == "nowcast", "nowcast", "forecast")
  fname <- paste0(file_prefix, "_", red, "_w", w)
  f <- file.path(in_folder, paste0(fname, ".RData"))
  
  
  if (!file.exists(f)) {
    warning(sprintf("File not found, skipping: %s", f))
    next
  }
  
  cat("Processing:", basename(f), "...\n")
  
  config <- list(label=fname, dataset=ds, reduction=red, window=w, h=h)
  
  e <- new.env()
  tryCatch({
    load(f, envir = e)
  }, error = function(err) {
    warning(paste("Failed to load", basename(f)))
  })
  
  target_obj <- NULL
  if ("res" %in% ls(e)) {
    target_obj <- e$res
  } else {
    for (nm in ls(e)) {
      obj <- e[[nm]]
      if (is.list(obj) && "penalized" %in% names(obj)) {
        target_obj <- obj
        break
      }
    }
  }
  
  if (is.null(target_obj)) {
    warning(paste("Skipping", basename(f), "- No valid result object found"))
    next
  }
  
  tryCatch({
    plots <- plot_macro_objects(target_obj, full_dates_vector, config, 
                                exclude_benchmarks = FALSE, 
                                y_zoom_limit = 1.0)
    
    combined_plot <- (plots$full | plots$zoom) +
      plot_layout(guides = "collect") & 
      theme(legend.position = "right")
    
    fn_combined <- file.path(out_folder, paste0(config$label, "_MSE_Combined.pdf"))
    suppressWarnings( # Warning due to the zoomed version, where the lines go out of frame
      ggsave(fn_combined, combined_plot, width = 16, height = 5, bg = "white")
    )
    
    cat("  ✓ Saved combined plots for", config$label, "\n")
    
  }, error = function(err) {
    warning(paste("Error plotting", config$label, ":", err$message))
  })
}

message("Macro Plot Generation Complete!")
message("=========================================================")