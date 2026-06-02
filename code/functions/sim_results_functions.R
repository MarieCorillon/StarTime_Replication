
# Utility functions for summarizing and visualizing simulation results

#' Calculate Average Metrics Across Simulations
calculate_average_metrics <- function(data) {
  metric_names <- rownames(data)
  
  base_metrics_full <- gsub("\\.?[0-9]+$", "", metric_names)
  base_metrics <- unique(base_metrics_full)
  
  methods <- colnames(data)[colnames(data) != "sim_id"]
  
  results <- matrix(NA, nrow = length(base_metrics), ncol = length(methods))
  rownames(results) <- base_metrics
  colnames(results) <- methods
  
  for (metric in base_metrics) {
    metric_rows <- which(base_metrics_full == metric)
    
    for (method in methods) {
      results[metric, method] <- mean(data[metric_rows, method], na.rm = TRUE)
    }
  }
  
  results_df <- as.data.frame(results)
  
  return(results_df)
}

relabel_result_table <- function(result_table, methods_to_keep, method_labels, metric_labels) {
  methods_present <- intersect(methods_to_keep, colnames(result_table))
  result_table <- result_table[, methods_present, drop = FALSE]
  
  colnames(result_table) <- method_labels[methods_present]
  
  metrics_present <- intersect(rownames(result_table), names(metric_labels))
  result_table <- result_table[metrics_present, , drop = FALSE]
  
  rownames(result_table) <- metric_labels[metrics_present]
  
  return(result_table)
}

# Computes shared y-axis upper limits for an n100+n200 pair.
# Returns a named list with 'ylims' and 'allow_top_bars' flag.
compute_shared_ylims <- function(res_list, is_mixed_dgp3 = FALSE) {
  metric_names <- rownames(res_list[[1]])
  ylims <- setNames(lapply(metric_names, function(metric) {
    
    # 1. ARI and F1: 0 to 1.05 (gives visual headroom so bars don't touch top)
    if (metric %in% c("ARI", "F1")) {
      return(c(0, 1.05))
    }
    
    # 2. MSE limits
    all_vals <- unlist(lapply(res_list, function(df) as.numeric(df[metric, ])))
    all_vals <- all_vals[is.finite(all_vals)]
    
    # Only use the 3rd-highest logic for Mixed DGP3, NOT AR DGP3
    if (metric == "MSE" && is_mixed_dgp3) {
      # For Mixed DGP3: use 3rd highest value to allow top 2 bars to reach ceiling
      max_val <- if(length(all_vals) >= 3) sort(all_vals, decreasing = TRUE)[3] else max(all_vals)
    } else {
      # Normal max (for AR DGP3 and all other DGPs)
      max_val <- if(length(all_vals) > 0) max(all_vals) else 1
    }
    
    c(0, max_val)  
  }), metric_names)
  
  # Determine if this DGP pair is Mixed DGP3 n=100 (special case for bar height)
  pair_names <- names(res_list)
  allow_top_bars <- any(pair_names == "Mixed_DGP_3_n100")
  
  list(ylims = ylims, allow_top_bars = allow_top_bars)
}

#' Generate Barplot of Averaged Simulation Metrics
generate_AR_barplot <- function(averaged_metrics, plot_title, dgp_name, ylims = NULL, allow_top_bars = FALSE) {  
  
  averaged_metrics$metric <- rownames(averaged_metrics)
  
  plot_data <- averaged_metrics %>%
    pivot_longer(cols = -metric, names_to = "method", values_to = "value") %>%
    mutate(
      metric = factor(metric, levels = c("MSE", "ARI", "F1")),
      method = factor(method, levels = colnames(averaged_metrics)[colnames(averaged_metrics) != "metric"]),
      value = ifelse(is.nan(value) | is.na(value), 0, value)
    )
  
  color_palette <- c("#648FFF", "#DC267F", "#FE6100", "#FFB000",
                     "#785EF0", "#7F7F7F", "#4ECDC4", "#FF6F61")
  color_palette <- setNames(
    color_palette[seq_along(colnames(averaged_metrics)[colnames(averaged_metrics) != "metric"])],
    colnames(averaged_metrics)[colnames(averaged_metrics) != "metric"]
  )
  
  if (!is.null(ylims)) {
    lim_df <- data.frame(
      metric = factor(names(ylims), levels = levels(plot_data$metric)),
      limit = sapply(ylims, `[`, 2, USE.NAMES = FALSE)
    )
    
    mse_mult <- if(allow_top_bars) 1 else 1.05
    
    plot_data <- plot_data %>%
      left_join(lim_df, by = "metric") %>%
      mutate(
        panel_top = ifelse(metric %in% c("ARI", "F1"), limit, limit * mse_mult),
        value = ifelse(value > panel_top, panel_top, value)
      ) %>%
      select(-limit, -panel_top)
    
    lim_data <- lim_df %>%
      mutate(
        method = levels(plot_data$method)[1],
        value  = ifelse(metric %in% c("ARI", "F1"), limit, limit * mse_mult)
      )
    blank_layer <- geom_blank(data = lim_data)
  } else {
    blank_layer <- NULL
  }
  
  p <- ggplot(plot_data, aes(x = method, y = value, fill = method)) +
    blank_layer +                                     
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(rows = vars(metric), scales = "free_y") +
    scale_fill_manual(values = color_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) + 
    labs(title = plot_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title       = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      axis.line.x      = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 16),
      axis.text.y      = element_text(size = 16),
      axis.title.y     = element_text(size = 18, face = "bold"),
      panel.grid.major = element_line(color = "grey70"),
      strip.text       = element_text(face = "bold", size = 20),
      plot.background  = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.spacing    = unit(1.5, "lines") 
    )
  
  invisible(p)
}

#' Generate bar plot for Mixed simulation results
generate_Mixed_barplot <- function(averaged_metrics, plot_title, dgp_name, ylims = NULL, allow_top_bars = FALSE) {  
  
  averaged_metrics$metric <- rownames(averaged_metrics)
  
  plot_data <- averaged_metrics %>%
    pivot_longer(cols = -metric, names_to = "method", values_to = "value") %>%
    mutate(
      metric = factor(metric, levels = c("MSE", "ARI", "F1")),
      method = factor(method, levels = colnames(averaged_metrics)[colnames(averaged_metrics) != "metric"]),
      value = ifelse(is.nan(value) | is.na(value), 0, value)
    )
  
  color_palette <- c("#648FFF", "#DC267F", "#FE6100", "#FFB000",
                     "#4ECDC4", "#785EF0", "#7F7F7F", "#FF6F61")
  color_palette <- setNames(
    color_palette[seq_along(colnames(averaged_metrics)[colnames(averaged_metrics) != "metric"])],
    colnames(averaged_metrics)[colnames(averaged_metrics) != "metric"]
  )
  
  if (!is.null(ylims)) {
    lim_df <- data.frame(
      metric = factor(names(ylims), levels = levels(plot_data$metric)),
      limit = sapply(ylims, `[`, 2, USE.NAMES = FALSE)
    )
    
    mse_mult <- if(allow_top_bars) 1 else 1.05
    
    plot_data <- plot_data %>%
      left_join(lim_df, by = "metric") %>%
      mutate(
        panel_top = ifelse(metric %in% c("ARI", "F1"), limit, limit * mse_mult),
        value = ifelse(value > panel_top, panel_top, value)
      ) %>%
      select(-limit, -panel_top)
    
    lim_data <- lim_df %>%
      mutate(
        method = levels(plot_data$method)[1],
        value  = ifelse(metric %in% c("ARI", "F1"), limit, limit * mse_mult)
      )
    blank_layer <- geom_blank(data = lim_data)
  } else {
    blank_layer <- NULL
  }
  
  p <- ggplot(plot_data, aes(x = method, y = value, fill = method)) +
    blank_layer +                                      
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_grid(rows = vars(metric), scales = "free_y") +
    scale_fill_manual(values = color_palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) + 
    labs(title = plot_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title       = element_text(size = 22, face = "bold", hjust = 0.5),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      axis.line.x      = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 16),
      axis.text.y      = element_text(size = 16),
      axis.title.y     = element_text(size = 18, face = "bold"),
      panel.grid.major = element_line(color = "grey70"),
      strip.text       = element_text(face = "bold", size = 20),
      plot.background  = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.spacing    = unit(1.5, "lines") 
    )
  
  invisible(p)
}

plot_all_sim_results <- function(final_results, dgp_title_map, plot_fun = generate_AR_barplot, save_dir = NULL) {
  plot_list <- list()
  
  base_names <- unique(sub("_n(100|200)$", "", names(final_results)))
  
  for (dgp_base in base_names) {
    pair_names <- grep(paste0("^", dgp_base, "_n(100|200)$"),
                       names(final_results), value = TRUE)
    
    ylim_info <- compute_shared_ylims(
      final_results[pair_names],
      is_mixed_dgp3 = grepl("Mixed_DGP_3", dgp_base)
    )
    ylims <- ylim_info$ylims
    
    for (res_name in pair_names) {
      # Only allow top bars to touch ceiling for Mixed DGP 3 n=100.
      # n=200 and all other DGPs get normal 5% headroom.
      current_allow_top_bars <- (res_name == "Mixed_DGP_3_n100")
      
      p <- plot_fun(final_results[[res_name]], dgp_title_map[[res_name]],
                    res_name, ylims = ylims, allow_top_bars = current_allow_top_bars)
      plot_list[[res_name]] <- p
      
      if (!is.null(save_dir)) {
        if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
        file_path <- file.path(save_dir, paste0(res_name, "_metrics.pdf"))
        
        ggplot2::ggsave(
          filename = file_path, 
          plot = p, 
          width = 10, 
          height = 8, 
          dpi = 300,
          bg = "white"
        )
      }
    }
  }
  
  invisible(plot_list)
}

