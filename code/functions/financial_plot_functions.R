get_har_labels <- function(coef_length) {
  if (coef_length == 20) {
    labs <- rep(NA, 20)
    labs[1]      <- "lag1"
    labs[2:5]    <- "week"
    labs[6:20]   <- "month"
    return(labs)
  }
  if (coef_length == 40) {
    labs <- rep(NA, 40)
    labs[1]      <- "lag1"
    labs[2:5]    <- "week"
    labs[6:20]   <- "month"
    labs[21:40]  <- "zero"
    return(labs)
  }
  stop("HAR labels only defined for 20 or 40")
}

compute_ari_grid <- function(results) {
  stocks <- names(results)
  out <- list()
  for (stock in stocks) {
    coef_list <- results[[stock]]$startime_coefficients
    for (step in seq_along(coef_list)) {
      coefs <- as.numeric(coef_list[[step]])
      true_labels <- get_har_labels(length(coefs))
      pred_labels <- coefs
      ari_val <- adjusted_rand_index(true_labels, pred_labels)
      out[[length(out) + 1]] <- data.frame(stock = stock, step  = step, ari   = ari_val)
    }
  }
  do.call(rbind, out)
}

add_ari_dates <- function(ari_grid, variances_df, max_lag, window_size, h) {
  all_dates <- as.Date(variances_df[[1]])
  n_steps <- nrow(subset(ari_grid, stock == ari_grid$stock[1])) 
  if (n_steps > length(all_dates)) stop("n_steps > length(all_dates); check that Dates match the series.")
  pred_dates <- tail(all_dates, n_steps)
  ari_grid$Date <- pred_dates[ari_grid$step]
  ari_grid
}

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

compute_ari_vector <- function(coef_list, zero_threshold = 1e-6, value_tol = 1e-8) {
  ari_vals <- numeric(length(coef_list))
  for (i in seq_along(coef_list)) {
    coefs <- as.numeric(coef_list[[i]])
    true_labels <- get_har_labels(length(coefs))
    pred_labels <- coefs
    ari_vals[i] <- adjusted_rand_index(true_labels, pred_labels)
  }
  ari_vals
}

prepare_metric_df <- function(results, variances_df, window_size, h, metric_type = c("mse", "qlike", "ari", "freq"), models_to_plot = NULL) {
  metric_type <- match.arg(metric_type)
  if (metric_type %in% c("mse", "qlike")) {
    model_keys   <- c("predictions_ar1", "predictions_har", "predictions_rw", "predictions_simple_startime", "predictions_startime")
    pretty_names <- c("AR(1)", "HAR", "Random Walk", "Simple StarTime", "Post StarTime")
  } else {
    model_keys   <- c("startime_coefficients") 
    pretty_names <- c("Post StarTime")
  }
  if (!is.null(models_to_plot)) {
    keep_idx <- pretty_names %in% models_to_plot
    if (!any(keep_idx)) stop("No matching models found in models_to_plot.")
    model_keys   <- model_keys[keep_idx]
    pretty_names <- pretty_names[keep_idx]
  }
  df_list <- list()
  for (stock in names(results)) {
    true_vals <- results[[stock]]$true_values
    if (is.null(true_vals)) next 
    n_preds   <- length(true_vals)
    base_date_vec  <- variances_df[[1]][(window_size + h):(window_size + h + n_preds - 1)]
    for (j in seq_along(model_keys)) {
      current_date_vec <- base_date_vec
      if (metric_type %in% c("mse", "qlike")) {
        preds <- results[[stock]][[model_keys[j]]]
        if (is.null(preds)) next
        value <- if (metric_type == "mse") (true_vals - preds)^2 else (true_vals / preds) - log(true_vals / preds) - 1
      } else {
        coef_list <- results[[stock]][[model_keys[j]]]
        if (is.null(coef_list)) next
        ari_vals  <- compute_ari_vector(coef_list)
        value <- if (metric_type == "ari") ari_vals else as.numeric(ari_vals >= 0.90)
        len_diff <- length(value) - length(current_date_vec)
        if (len_diff > 0) value <- tail(value, length(current_date_vec)) else if (len_diff < 0) current_date_vec <- tail(current_date_vec, length(value))
      }
      df_list[[length(df_list) + 1]] <- data.frame(date = as.Date(current_date_vec), stock = stock, model = pretty_names[j], value = value)
    }
  }
  if (length(df_list) == 0) stop("No data frame created. Check results structure.")
  dplyr::bind_rows(df_list)
}

smoothing_function <- function(data, smoothing_window){
  n <- length(data)
  smoothed_data <- rep(NA_real_, n)
  half_window <- floor(smoothing_window / 2)
  start_loop <- half_window + 1
  end_loop   <- n - half_window
  if (n < smoothing_window) return(smoothed_data)
  for (i in start_loop:end_loop){
    start_index <- i - half_window
    end_index   <- i + half_window
    window_data <- data[start_index:end_index]
    smoothed_data[i] <- mean(window_data, na.rm = TRUE)
  }
  smoothed_data
}

plot_metric_refined <- function(results, variances_df, window_size, h, metric_type = c("mse", "qlike", "ari", "freq"), models_to_plot = NULL, plot_title = NULL, center_stat = c("mean", "median"), band_type = c("none", "minmax", "iqr", "deciles"), y_limit_quantile = 0.99, smooth_daily = FALSE, smoothing_window = 20, relative_to_rw = FALSE) { 
  metric_type <- match.arg(metric_type)
  center_stat <- match.arg(center_stat)
  band_type   <- match.arg(band_type)
  fetch_models <- models_to_plot
  if (relative_to_rw && !is.null(models_to_plot) && !"Random Walk" %in% models_to_plot) fetch_models <- c(models_to_plot, "Random Walk")
  df <- prepare_metric_df(results, variances_df, window_size, h, metric_type, fetch_models)
  if (relative_to_rw) {
    df <- df |> tidyr::pivot_wider(names_from = model, values_from = value) |> dplyr::mutate(dplyr::across(where(is.numeric), ~ .x / `Random Walk`)) |> tidyr::pivot_longer(cols = -c(date, stock), names_to = "model", values_to = "value") |> dplyr::filter(is.finite(value)) |> dplyr::filter(model != "Random Walk")
  }
  daily_agg <- df |> dplyr::group_by(date, model) |> dplyr::summarise(center_val = if (center_stat == "mean") mean(value, na.rm=TRUE) else median(value, na.rm=TRUE), lower = switch(band_type, "none" = NA_real_, "minmax" = min(value, na.rm=TRUE), "iqr" = quantile(value, 0.25, na.rm=TRUE), "deciles" = quantile(value, 0.10, na.rm=TRUE)), upper = switch(band_type, "none" = NA_real_, "minmax" = max(value, na.rm=TRUE), "iqr" = quantile(value, 0.75, na.rm=TRUE), "deciles" = quantile(value, 0.90, na.rm=TRUE)), .groups = "drop")
  if (smooth_daily) {
    daily_agg <- daily_agg |> dplyr::group_by(model) |> dplyr::mutate(smoothed_val = smoothing_function(center_val, smoothing_window = smoothing_window)) |> dplyr::ungroup()
  } else {
    daily_agg$smoothed_val <- daily_agg$center_val 
  }
  calc_ylim <- function(data_table) {
    if (is.null(y_limit_quantile) || nrow(data_table) == 0) return(NULL) 
    vals <- if (band_type == "none") data_table$smoothed_val else data_table$upper
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) return(NULL)
    c(0, quantile(vals, y_limit_quantile, na.rm = TRUE) * 1.05)
  }
  y_label <- paste(switch(metric_type, mse="MSFE", qlike="QLIKE", ari="ARI", freq="Frequency"))
  if (is.null(plot_title)) plot_title <- paste(y_label, "- W", window_size, "H", h)
  cols <- c("AR(1)" = "#e41a1c", "HAR" = "#377eb8", "Random Walk" = "#4daf4a", "Simple StarTime" = "#984ea3", "Post StarTime" = "#ff7f00")
  y_limits <- calc_ylim(daily_agg)
  p <- ggplot(daily_agg, aes(x = date, group = model, color = model)) + labs(title = plot_title, y = y_label, x = "Date")
  if (relative_to_rw) p <- p + geom_hline(yintercept = 1, linetype = "dashed", color = "black", alpha = 0.5)
  if (band_type != "none") p <- p + geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.15, color = NA) + scale_fill_manual(values = cols, guide = "none")
  if (smooth_daily) p <- p + geom_line(aes(y = smoothed_val), linewidth = 1.0) else p <- p + geom_line(aes(y = center_val), linewidth = 0.8)
  p <- p + scale_color_manual(values = cols) + scale_x_date(date_labels = "%b %Y", date_breaks = "6 months") + my_theme()
  if (!is.null(y_limits)) p <- p + coord_cartesian(ylim = y_limits)
  invisible(list(plot = p, data = daily_agg))
}

stack_plots <- function(plot_top, plot_bottom, overall_title = NULL, overall_subtitle = NULL, legend_title = "Models") {
  if (is.list(plot_top) && !is_ggplot(plot_top)) plot_top <- plot_top$plot
  if (is.list(plot_bottom) && !is_ggplot(plot_bottom)) plot_bottom <- plot_bottom$plot
  plot_top <- plot_top + labs(color = legend_title, fill = legend_title, title = NULL, x = NULL) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), legend.position = "right")
  plot_bottom <- plot_bottom + labs(title = NULL) + theme(legend.position = "none")
  combined <- plot_top / plot_bottom + plot_layout(heights = c(1, 1))
  if (!is.null(overall_title) || !is.null(overall_subtitle)) {
    combined <- combined + plot_annotation(title = overall_title, subtitle = overall_subtitle, theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10))))
  }
  combined
}

prepare_stime_coef_matrix <- function(coef_list, dates_vec = NULL, tol = 1e-10) {
  mat <- do.call(rbind, lapply(coef_list, function(v) as.numeric(v[, 1])))
  df <- as.data.frame(mat)
  colnames(df) <- paste0("Lag", seq_len(ncol(df)))
  df$Time <- seq_len(nrow(df))
  df$Date <- if (!is.null(dates_vec)) as.Date(dates_vec) else as.Date(NA)
  long <- pivot_longer(df, cols = starts_with("Lag"), names_to = "Lag", values_to = "Value", names_prefix = "Lag")
  long$Lag <- as.integer(long$Lag)
  long <- mutate(long, is_zero = abs(Value) < tol)
  long
}

plot_stime_coef_matrix <- function(tidy_coef, title = NULL, subtitle = NULL, y_label_format = c("year", "ym")) {
  y_label_format <- match.arg(y_label_format)
  tidy_coef <- tidy_coef %>% dplyr::group_by(Time) %>% dplyr::mutate(RowGroup = ifelse(is_zero, 0L, as.integer(factor(Value)))) %>% dplyr::ungroup()
  max_group <- max(tidy_coef$RowGroup, na.rm = TRUE)
  levels_all <- as.character(0:max_group)
  nz_n <- max(1L, max_group)
  # nz_cols <- if (nz_n <= 12) {
  #   c(RColorBrewer::brewer.pal(min(8, nz_n), "Set2"), if (nz_n > 8) RColorBrewer::brewer.pal(nz_n - 8, "Dark2") else NULL)
  # } else {
  #   grDevices::hcl.colors(nz_n, "Dark3")
  # }
  nz_cols <- if (nz_n <= 12) {
    if (nz_n <= 8) {
      RColorBrewer::brewer.pal(8, "Set2")[seq_len(nz_n)]
    } else {
      c(
        RColorBrewer::brewer.pal(8, "Set2"),
        RColorBrewer::brewer.pal(8, "Dark2")[seq_len(nz_n - 8)]
      )
    }
  } else {
    grDevices::hcl.colors(nz_n, "Dark3")
  }
  pal <- c("grey70", nz_cols)
  names(pal) <- levels_all
  lab_df <- tidy_coef %>% dplyr::distinct(Time, Date) %>% dplyr::arrange(Time) %>% dplyr::mutate(YM = ifelse(!is.na(Date), format(Date, ifelse(y_label_format == "ym", "%Y-%m", "%Y")), NA_character_))
  if (any(!is.na(lab_df$YM))) {
    first_of_period <- lab_df %>% dplyr::filter(!is.na(YM)) %>% dplyr::group_by(YM) %>% dplyr::slice_head(n = 1) %>% dplyr::ungroup()
    y_breaks <- first_of_period$Time
    y_labels <- first_of_period$YM
    y_axis_title <- NULL
  } else {
    y_breaks <- pretty(lab_df$Time, n = 10)
    y_labels <- y_breaks
    y_axis_title <- "Time Step"
  }
  ggplot(tidy_coef, aes(x = Lag, y = Time, fill = factor(RowGroup, levels = levels_all))) + geom_tile(alpha = 0.95) + scale_fill_manual(values = pal, guide = "none") + scale_x_continuous(breaks = seq(1, max(tidy_coef$Lag, na.rm = TRUE), by = 1), labels = seq(1, max(tidy_coef$Lag, na.rm = TRUE), by = 1)) + scale_y_continuous(breaks = y_breaks, labels = y_labels) + labs(title = title, subtitle = subtitle, x = "Lag", y = y_axis_title) + coord_cartesian(expand = FALSE, clip = "off") + theme_minimal() + theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5), plot.subtitle = element_text(size = 14, hjust = 0.5, margin = margin(b = 15)), axis.text.x = element_text(angle = 90, vjust = 0.5, size = 14, face = "bold"), axis.text.y = element_text(size = 14, face = "bold"), axis.title.x = element_text(size = 16, face = "bold"), axis.title.y = element_text(size = 16, face = "bold"))
}