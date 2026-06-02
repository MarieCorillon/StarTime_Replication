
# 1. SHARED THEME
my_theme <- function(base_size = 16) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor = element_line(color = "grey92", linewidth = 0.2),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y  = element_text(size = 16),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 16, face = "bold"),
      plot.title   = element_text(size = 22, face = "bold", hjust = 0.5),
      legend.position = "bottom"
    )
}


# 2. SHARED COLOR PALETTE & LABELS
macro_colors <- c(
  "AR1"             = "#e41a1c",
  "HAR"             = "#377eb8",
  "Random_Walk"     = "#4daf4a",
  "Simple_StarTime" = "#984ea3",
  "StarTime"        = "#ff7f00",
  "MIDASML"         = "#a65628"
)

macro_labels <- c(
  "AR1"             = "AR(1)",
  "Random_Walk"     = "Random Walk",      
  "Simple_StarTime" = "Simple StarTime",   
  "StarTime"        = "Post StarTime",
  "MIDASML"         = "MIDAS-ML"
)


# 3. PLOT MACRO OBJECTS
plot_macro_objects <- function(res_object, full_dates, config, 
                               exclude_benchmarks = FALSE,
                               y_zoom_limit = 1.0) {
  
  y_true <- as.numeric(res_object$penalized$evaluation$y_true)
  
  preds_list <- list(
    StarTime        = as.numeric(res_object$penalized$evaluation$pred_star_post),
    Simple_StarTime = as.numeric(res_object$penalized$results$simple_predictions),
    MIDASML         = as.numeric(res_object$penalized$results$midasml_predictions)
  )
  
  # Only add benchmarks if we are looking at raw forecasts
  if (config$dataset == "raw"|| config$dataset == "forecast") {
    preds_list$Random_Walk <- as.numeric(res_object$random_walk$evaluation$pred)
    preds_list$AR1         <- as.numeric(res_object$ar1$evaluation$pred)
  }
  
  w <- config$window
  h <- config$h
  
  test_indices <- seq(from = w + h, to = length(full_dates), by = 1)
  if (length(y_true) != length(test_indices)) {
    test_dates <- tail(full_dates, length(y_true))
  } else {
    test_dates <- full_dates[test_indices]
  }
  test_dates <- tryCatch(as.Date(test_dates), error = function(e) test_dates)
  
  if (config$dataset == "raw" || config$dataset == "forecast") {
    y_true <- y_true[-length(y_true)]
    test_dates <- test_dates[-length(test_dates)]
    for (nm in names(preds_list)) {
      preds_list[[nm]] <- preds_list[[nm]][-length(preds_list[[nm]])]
    }
  }
  
  df_list <- lapply(names(preds_list), function(nm) {
    data.frame(
      date = test_dates,
      model = nm,
      value = ((y_true - preds_list[[nm]])^2) * 1000, 
      stringsAsFactors = FALSE
    )
  })
  metric_df <- do.call(rbind, df_list)
  
  metric_df$model <- factor(metric_df$model, levels = names(macro_colors))
  
  if (exclude_benchmarks) {
    metric_df <- metric_df %>% 
      filter(!model %in% c("Random_Walk", "AR1"))
  }
  
  # COMMON LAYERS
  common_layers <- list(
    my_theme(),
    scale_color_manual(values = macro_colors, labels = macro_labels, name = "Models"), 
    scale_x_date(date_breaks = "2 years", date_labels = "%Y", expand = c(0, 0)),
    theme(
      legend.position = "right",
      axis.text.x = element_text(angle = 0, hjust = 0.5) 
    ) 
  )
  
  # PLOT 1: Full History
  p_full <- ggplot(metric_df, aes(x = date, y = value, color = model)) +
    geom_line(linewidth = 1) +
    labs(x = NULL, y = "Squared Forecast Errors") + 
    common_layers
  
  # PLOT 2: Censored Zoom
  p_zoom <- ggplot(metric_df, aes(x = date, y = value, color = model)) +
    geom_line(linewidth = 1) +
    labs(x = NULL, y = NULL) + 
    common_layers +
    coord_cartesian(ylim = c(0, y_zoom_limit))
  
  return(list(full = p_full, zoom = p_zoom))
}