source("code/functions/financial_plot_functions.R")
source("code/functions/metrics.R")

message("=========================================================")
message("Starting Financial Plot Generation...")

# Set to English to guarantee axis dates match the paper exactly ("Jan 2012", etc.)
Sys.setlocale("LC_TIME", "C")

in_folder  <- file.path(RESULTS_DIR, "applications/financial")
out_folder_lines <- file.path(FIGURES_DIR, "applications/financial/lines")
out_folder_coefs <- file.path(FIGURES_DIR, "applications/financial/coefficients")
data_path <- "data/financial/Variances_10min.csv"

if (!dir.exists(out_folder_lines)) dir.create(out_folder_lines, recursive = TRUE)
if (!dir.exists(out_folder_coefs)) dir.create(out_folder_coefs, recursive = TRUE)

if (!file.exists(data_path)) stop("Data file missing: ", data_path)
Variances_10min <- readr::read_csv(data_path, show_col_types = FALSE, name_repair = "minimal")

# Determine Configurations
if (REPLICATE_PAPER_RESULTS) {
  plot_configs <- expand.grid(window = c(125, 250, 1000), lag = c(20, 40), h = c(1, 20), stringsAsFactors = FALSE)
} else {
  plot_configs <- expand.grid(window = FIN_WINDOW_SIZE, lag = FIN_MAX_LAG, h = FIN_H, stringsAsFactors = FALSE)
}

coef_stocks <- c("JPM", "NKE")

for (i in 1:nrow(plot_configs)) {
  w   <- plot_configs$window[i]
  lag <- plot_configs$lag[i]
  h   <- plot_configs$h[i]
  
  file_name <- sprintf("Financial_lag%d_h%d_window%d.RData", lag, h, w)
  file_path <- file.path(in_folder, file_name)
  
  if (!file.exists(file_path)) {
    warning(sprintf("Skipping Financial Plots for W=%d, Lag=%d, H=%d (File missing: %s)", w, lag, h, file_name))
    next
  }
  
  message(sprintf("-> Loading results for W=%d, Lag=%d, H=%d...", w, lag, h))
  env <- new.env()
  load(file_path, envir = env)
  results <- env$results
  
  # 1. GENERATE STACKED LINE PLOTS (MSE + ARI)
  title_win_str <- if (w == 125) "Half-Year" else if (w == 250) "Year" else paste(w, "Day")
  
  tryCatch({
    mse_plots <- plot_metric_refined(
      results, Variances_10min, w, h, 
      metric_type = "mse", center_stat = "median", smooth_daily = TRUE, 
      y_limit_quantile = NULL, smoothing_window = 40, relative_to_rw = FALSE
    )
    
    ari_plots <- plot_metric_refined(
      results, Variances_10min, w, h, 
      metric_type = "ari", center_stat = "median", band_type = "deciles", 
      smooth_daily = TRUE, y_limit_quantile = NULL, smoothing_window = 40
    )
    
    final_plot <- stack_plots(
      plot_top = mse_plots,
      plot_bottom = ari_plots,
      overall_title = sprintf("Forecast Performance: %s Window", title_win_str),
      overall_subtitle = "Smoothed Daily MSFE vs Smoothed Daily ARI",
      legend_title = "Models"
    )
    
    pdf_line_name <- file.path(out_folder_lines,
                               sprintf("stacked_lines_lag%d_window%d_h%d.pdf", lag, w, h))
    suppressWarnings( # Lines extend beyond axis limits in zoomed plots
      ggsave(filename = pdf_line_name, plot = final_plot, width = 10.26, height = 5.98, bg = "white")
    )
    message("   ✓ Saved Stacked Lines: ", basename(pdf_line_name))
  }, error = function(e) {
    warning("   ! Failed to generate stacked lines: ", e$message)
  })
  
  # GENERATE COEFFICIENT MATRIX PLOTS (JPM & NKE Only)
  for (stock in coef_stocks) {
    if (is.null(results[[stock]])) {
      warning(sprintf("   ! Skipping Coefficient Plot for %s (Stock not found in results)", stock))
      next
    }
    
    tryCatch({
      start_index <- w + h
      end_index <- w + h + length(results[[stock]]$startime_coefficients) - 1
      dates_vec <- as.Date(Variances_10min[[1]][start_index:end_index])
      
      tidy_coef <- prepare_stime_coef_matrix(results[[stock]]$startime_coefficients, dates_vec = dates_vec)
      
      p_coef <- plot_stime_coef_matrix(
        tidy_coef,
        title = sprintf("Grouping Performance: %s Window", title_win_str),
        subtitle = stock,
        y_label_format = "year"
      )
      
      pdf_coef_name <- file.path(out_folder_coefs,
                                 sprintf("%s_lag%d_window%d_h%d.pdf", tolower(stock), lag, w, h))
      suppressWarnings( # Lines extend beyond axis limits in zoomed plots
        ggsave(pdf_coef_name, plot = p_coef, width = 10, height = 8, dpi = 200, bg = "white")
      )
      message("   ✓ Saved Coefficient Matrix: ", basename(pdf_coef_name))
    }, error = function(e) {
      warning("   ! Failed to generate coefficient plot for ", stock, ": ", e$message)
    })
  }
}

message("Financial Plot Generation Complete!")
message("=========================================================")
