# Extract and summarize StarTime coefficients from macro applications

message("=========================================================")
message("Starting Macro Coefficient Summaries...")

# 1. Core function: summarize beta activity and grouping
summarize_beta_activity_grouping <- function(betas_list, feature_names, tol = .Machine$double.eps^0.5) {
  betas_named <- lapply(betas_list, function(beta_mat) {
    beta_vec <- as.numeric(beta_mat)
    names(beta_vec) <- feature_names
    beta_vec
  })
  betas_mat <- do.call(cbind, betas_named)
  
  base_var <- sub("_(lag|wlag|mlag|qlag)[0-9]+$", "", rownames(betas_mat))
  vars <- unique(base_var)
  n_steps <- ncol(betas_mat)
  
  all_lags_grouped_mat <- matrix(FALSE, nrow = length(vars), ncol = n_steps,
                                 dimnames = list(vars, NULL))
  some_lags_grouped_mat <- matrix(FALSE, nrow = length(vars), ncol = n_steps,
                                  dimnames = list(vars, NULL))
  
  for(j in 1:n_steps) {
    col_j <- betas_mat[, j]
    grouped_by_regressor <- tapply(col_j, base_var, function(x) {
      nz <- abs(x) >= tol
      vals <- x[nz]
      n_nz <- length(vals)
      if(n_nz < 2) return(c(all = FALSE, some = FALSE)) 
      
      all_eq <- all(abs(vals - vals[1]) < tol)
      tab <- table(round(vals / tol))
      some_eq <- any(tab >= 2)
      
      c(all = all_eq, some = some_eq & !all_eq)
    })
    
    if(is.list(grouped_by_regressor)) {
      res_mat <- do.call(rbind, grouped_by_regressor)
      common <- intersect(rownames(all_lags_grouped_mat), rownames(res_mat))
      all_lags_grouped_mat[common, j]  <- res_mat[common, "all"]
      some_lags_grouped_mat[common, j] <- res_mat[common, "some"]
    } else if (is.matrix(grouped_by_regressor)) {
      if(ncol(grouped_by_regressor) == 2) {
        common <- intersect(rownames(all_lags_grouped_mat), rownames(grouped_by_regressor))
        all_lags_grouped_mat[common, j]  <- grouped_by_regressor[common, "all"]
        some_lags_grouped_mat[common, j] <- grouped_by_regressor[common, "some"]
      }
    }
  }
  
  grouped_counts_all <- rowSums(all_lags_grouped_mat)
  proportion_grouped_all <- grouped_counts_all / n_steps
  grouped_counts_some <- rowSums(some_lags_grouped_mat)
  proportion_grouped_some <- grouped_counts_some / n_steps
  
  coef_zero_all_steps <- apply(betas_mat, 1, function(x) all(abs(x) < tol))
  nonzero_counts <- rowSums(abs(betas_mat) >= tol)
  
  coef_summary_lag <- data.frame(
    variable = rownames(betas_mat),
    base_var = base_var,
    always_zero = coef_zero_all_steps,
    nonzero_steps = nonzero_counts,
    total_steps = n_steps,
    proportion_nonzero = nonzero_counts / n_steps,
    row.names = NULL
  )
  
  coef_summary_regressor <- data.frame(
    base_var = vars,
    all_lags_grouped_steps = grouped_counts_all,
    proportion_all_grouped = proportion_grouped_all,
    some_lags_grouped_steps = grouped_counts_some,
    proportion_some_grouped = proportion_grouped_some,
    total_steps = n_steps,
    row.names = NULL
  )
  
  list(
    lag_level = coef_summary_lag,
    regressor_level = coef_summary_regressor,
    all_lags_grouped_matrix = all_lags_grouped_mat,
    some_lags_grouped_matrix = some_lags_grouped_mat
  )
}

# 2. Batch wrapper
parse_config <- function(fname) {
  base <- sub("\\.RData$", "", fname)
  parts <- strsplit(base, "_")[[1]]
  h_part <- grep("^h[0-9]+", parts, value = TRUE)
  w_part <- grep("^w[0-9]+", parts, value = TRUE)
  h <- if(length(h_part)>0) as.numeric(sub("h", "", h_part[1])) else NA
  w <- if(length(w_part)>0) as.numeric(sub("w", "", w_part[1])) else NA
  
  # Distinguish "nowcast" vs "forecast"
  type <- if("nowcast" %in% parts) "nowcast" else "forecast"
  dataset <- if("full" %in% parts) "full" else if("reduced" %in% parts) "reduced" else "other"
  
  list(filename = fname, type = type, dataset = dataset, h = h, w = w, label=base)
}

summarize_beta_batch <- function(save_dir = file.path(RESULTS_DIR, "applications/macro/coefficient_summaries")) {
  
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  all_regressor_summaries <- list()
  all_lag_summaries <- list()
  
  # Dynamically discover files based on current mode
  files <- list.files(
    path = file.path(RESULTS_DIR, "applications/macro"),
    pattern = "\\.RData$",
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    warning("No macro application RData files found in: ", file.path(RESULTS_DIR, "applications/macro"))
    return(invisible(NULL))
  }
  
  for (f in files) {
    if (!file.exists(f)) {
      warning(sprintf("File not found, skipping: %s", f))
      next
    }
    
    fname <- basename(f)
    cat("Processing:", fname, "...\n")
    
    env <- new.env()
    tryCatch({ load(f, envir = env) }, error = function(e) {
      cat("  -> Skipping (load failed)\n")
      return(NULL)
    })
    
    # Check for result object structure
    res <- NULL
    if ("res" %in% ls(env)) {
      res <- env$res
    } else {
      for (nm in ls(env)) {
        obj <- env[[nm]]
        if (is.list(obj) && "penalized" %in% names(obj)) {
          res <- obj
          break
        }
      }
    }
    
    if (is.null(res) || is.null(res$penalized$results$betas)) {
      cat("  -> Skipping (Invalid structure)\n")
      next
    }
    
    cfg <- parse_config(fname)
    betas_list <- res$penalized$results$betas
    
    # Trim last observation for forecast files
    if (cfg$type == "forecast") {
      if (length(betas_list) > 0) {
        betas_list <- betas_list[-length(betas_list)]
      }
    }
    
    first_beta <- betas_list[[1]]
    if (is.null(dim(first_beta))) {
      feat_names <- names(first_beta)
    } else {
      feat_names <- rownames(first_beta)
      if (is.null(feat_names)) feat_names <- names(first_beta[,1])
    }
    
    if(is.null(feat_names)) {
      warning("  -> Feature names missing for ", fname)
      next
    }
    
    summary_out <- tryCatch({
      summarize_beta_activity_grouping(betas_list, feature_names = feat_names)
    }, error = function(e) {
      warning("  -> Failed: ", e$message)
      return(NULL)
    })
    
    if (is.null(summary_out)) next
    
    out_file <- file.path(save_dir, paste0(cfg$label, "_beta_summary.RData"))
    save(summary_out, cfg, file = out_file)
    
    reg_df <- summary_out$regressor_level %>%
      mutate(filename = cfg$filename, type = cfg$type, dataset = cfg$dataset, h = cfg$h, w = cfg$w)
    
    lag_df <- summary_out$lag_level %>%
      mutate(filename = cfg$filename, type = cfg$type, dataset = cfg$dataset, h = cfg$h, w = cfg$w)
    
    all_regressor_summaries[[fname]] <- reg_df
    all_lag_summaries[[fname]]       <- lag_df
  }
  
  if (length(all_regressor_summaries) > 0) {
    master_regressor <- do.call(rbind, all_regressor_summaries)
    master_lag       <- do.call(rbind, all_lag_summaries)
    
    write.csv(master_regressor, file.path(save_dir, "batch_regressor_summary.csv"), row.names = FALSE)
    write.csv(master_lag, file.path(save_dir, "batch_lag_summary.csv"), row.names = FALSE)
    
    cat("\nSaved individual .RData files and master CSVs to:", save_dir, "\n")
  }
}

# 3. Generate overviews
generate_overview <- function(reg_df, lag_df, dataset_filter, window_filter = NULL, type_filter = NULL) {
  
  sub_reg <- reg_df %>% dplyr::filter(dataset == dataset_filter)
  sub_lag <- lag_df %>% dplyr::filter(dataset == dataset_filter)
  
  lbl <- paste0(dataset_filter, "_set")
  
  if (!is.null(type_filter)) {
    sub_reg <- sub_reg %>% dplyr::filter(type == type_filter)
    sub_lag <- sub_lag %>% dplyr::filter(type == type_filter)
    lbl <- paste0(lbl, "_", type_filter)
  } else {
    lbl <- paste0(lbl, "_combined_types")
  }
  
  if (!is.null(window_filter)) {
    sub_reg <- sub_reg %>% dplyr::filter(w == window_filter)
    sub_lag <- sub_lag %>% dplyr::filter(w == window_filter)
    lbl <- paste0(lbl, "_w", window_filter)
  } else {
    lbl <- paste0(lbl, "_ALL_windows")
  }
  
  if(nrow(sub_reg) == 0) {
    warning("No data found for filter: ", lbl)
    return(NULL)
  }
  
  grouping_overview <- sub_reg %>%
    dplyr::ungroup() %>%
    dplyr::group_by(base_var) %>% 
    dplyr::summarise(
      grand_total_steps = sum(total_steps),
      total_all_grouped = sum(all_lags_grouped_steps),
      pct_all_grouped   = (sum(all_lags_grouped_steps) / sum(total_steps)) * 100,
      total_some_grouped = sum(some_lags_grouped_steps),
      pct_some_grouped   = (sum(some_lags_grouped_steps) / sum(total_steps)) * 100,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(pct_all_grouped))
  
  selection_overview <- sub_lag %>%
    dplyr::ungroup() %>%
    dplyr::group_by(variable, base_var) %>% 
    dplyr::summarise(
      grand_total_steps    = sum(total_steps),
      total_nonzero_steps  = sum(nonzero_steps),
      pct_selected         = (sum(nonzero_steps) / sum(total_steps)) * 100,
      is_always_zero       = sum(nonzero_steps) == 0,
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(pct_selected))
  
  list(
    label = lbl,
    selection_per_lag = selection_overview,
    grouping_per_var  = grouping_overview
  )
}

save_overview_rdata <- function(ov_list, output_dir) {
  if(is.null(ov_list)) return()
  fname <- file.path(output_dir, paste0("Overview_", ov_list$label, ".RData"))
  save(ov_list, file = fname)
  cat("Saved RData:", fname, "\n")
}

# 4. Execution
save_dir <- file.path(RESULTS_DIR, "applications/macro/coefficient_summaries")
message("-> Running batch extraction in: ", save_dir)
summarize_beta_batch(save_dir)

# Load master CSVs and generate overviews
reg_df <- NULL
lag_df <- NULL

master_reg_path <- file.path(save_dir, "batch_regressor_summary.csv")
master_lag_path <- file.path(save_dir, "batch_lag_summary.csv")

if (file.exists(master_reg_path)) {
  reg_df <- read.csv(master_reg_path)
}
if (file.exists(master_lag_path)) {
  lag_df <- read.csv(master_lag_path)
}

if (!is.null(reg_df) && !is.null(lag_df)) {
  out_dir <- file.path(save_dir, "overviews")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  message("-> Generating overviews...")
  
  # --- NOWCASTS ---
  save_overview_rdata(generate_overview(reg_df, lag_df, "reduced", type_filter = "nowcast"), out_dir)
  save_overview_rdata(generate_overview(reg_df, lag_df, "reduced", window_filter = 66, type_filter = "nowcast"), out_dir)
  save_overview_rdata(generate_overview(reg_df, lag_df, "reduced", window_filter = 105, type_filter = "nowcast"), out_dir)
  
  # --- FORECASTS ---
  save_overview_rdata(generate_overview(reg_df, lag_df, "reduced", type_filter = "forecast"), out_dir)
  save_overview_rdata(generate_overview(reg_df, lag_df, "reduced", window_filter = 66, type_filter = "forecast"), out_dir)
  save_overview_rdata(generate_overview(reg_df, lag_df, "reduced", window_filter = 105, type_filter = "forecast"), out_dir)
  
  cat("\n✓ Distinct beta summary and overview pipeline complete!\n")
} else {
  warning("Master CSVs not found. Skipping overviews.")
}

message("Macro Coefficient Summaries Complete!")
message("=========================================================")
