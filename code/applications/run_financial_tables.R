if (!requireNamespace("forecast", quietly = TRUE)) install.packages("forecast")

message("=========================================================")
message("Starting Financial Table Generation...")

in_folder  <- file.path(RESULTS_DIR, "applications/financial")
tests_dir  <- file.path(RESULTS_DIR, "applications/financial/tests")
tables_dir <- file.path(FIGURES_DIR, "applications/financial/tables")

if (!dir.exists(tests_dir))  dir.create(tests_dir,  recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)


# CORE STATISTICAL TEST FUNCTIONS
run_dm_for_results <- function(results, h, model_names, main_models, loss_metric = "MSE", alternative = "two.sided") {
  all_stocks <- names(results)
  dm_all     <- list()
  idx_stock  <- 1
  is_qlike   <- (loss_metric == "QLIKE")
  dm_power   <- if (loss_metric == "MAE") 1 else if (loss_metric == "MSE") 2 else 1
  
  for (stock in all_stocks) {
    res    <- results[[stock]]
    true_y <- res$true_values
    models <- list()
    for (m in model_names) if (!is.null(res[[m]])) models[[m]] <- res[[m]]
    if (length(models) < 2) next
    
    dm_list <- list()
    idx     <- 1
    for (main_name in names(models)) {
      for (base_name in names(models)) {
        if (main_name == base_name) next
        if (!is.null(main_models) && !(main_name %in% main_models)) next
        
        f1 <- models[[main_name]]; f2 <- models[[base_name]]
        n  <- min(length(true_y), length(f1), length(f2))
        y_a  <- true_y[1:n]; f1_a <- f1[1:n]; f2_a <- f2[1:n]
        
        if (is_qlike) {
          input1     <- (exp(y_a) / exp(f1_a)) - (y_a - f1_a) - 1
          input2     <- (exp(y_a) / exp(f2_a)) - (y_a - f2_a) - 1
          test_power <- 1
        } else {
          input1     <- y_a - f1_a
          input2     <- y_a - f2_a
          test_power <- dm_power
        }
        
        dm_obj <- tryCatch(
          dm.test(input1, input2, h = h, power = test_power, alternative = alternative),
          error = function(e) return(NULL)
        )
        if (is.null(dm_obj)) next
        
        dm_list[[idx]] <- data.frame(
          stock      = stock,      main_model = main_name, base_model = base_name,
          loss       = loss_metric, statistic  = as.numeric(dm_obj$statistic),
          p_value    = as.numeric(dm_obj$p.value), stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
    if (length(dm_list) > 0) { dm_all[[idx_stock]] <- do.call(rbind, dm_list); idx_stock <- idx_stock + 1 }
  }
  if (length(dm_all) == 0) return(NULL)
  do.call(rbind, dm_all)
}

run_mcs_for_all_stocks <- function(results, model_names, alpha, B, statistic = "TR", loss = "MSE") {
  all_stocks <- names(results)
  in_mcs     <- matrix(FALSE, nrow = length(model_names), ncol = length(all_stocks),
                       dimnames = list(model_names, all_stocks))
  stat_col   <- if (statistic == "TR") "MCS_R" else "MCS_M"
  
  for (j in seq_along(all_stocks)) {
    stock <- all_stocks[j]
    res   <- results[[stock]]
    y     <- res$true_values
    Fmat  <- sapply(model_names, function(m) {
      if (is.null(res[[m]])) return(rep(NA, length(y)))
      return(res[[m]])
    })
    
    Tn       <- min(length(y), nrow(Fmat))
    y_cut    <- y[1:Tn]; Fmat_cut <- Fmat[1:Tn, , drop = FALSE]
    valid    <- complete.cases(Fmat_cut) & !is.na(y_cut)
    y_final  <- y_cut[valid]; Fmat_final <- Fmat_cut[valid, , drop = FALSE]
    if (length(y_final) < 10) next
    
    if (loss == "MSE") {
      loss_mat <- (y_final - Fmat_final)^2
    } else {
      loss_mat <- (exp(y_final) / exp(Fmat_final)) - (y_final - Fmat_final) - 1
    }
    colnames(loss_mat) <- model_names
    
    mcs_res <- try(suppressWarnings(
      MCSprocedure(Loss = loss_mat, alpha = alpha, B = B, statistic = statistic, verbose = FALSE)
    ), silent = TRUE)
    if (inherits(mcs_res, "try-error")) next
    
    tab   <- mcs_res@show
    p_val <- tab[, stat_col]
    inc   <- p_val >= alpha
    in_mcs[names(inc), stock] <- inc
  }
  list(in_mcs = in_mcs, perc = rowMeans(in_mcs, na.rm = TRUE) * 100)
}

run_all_for_results <- function(results, config, model_names, main_models, summary_benchmark,
                                mcs_B, mcs_alpha, mcs_statistic, dm_alternative) {
  mcs_mse <- run_mcs_for_all_stocks(results, model_names, alpha = mcs_alpha, B = mcs_B,
                                    statistic = mcs_statistic, loss = "MSE")
  dm_mse  <- run_dm_for_results(results, h = config$h, model_names = model_names,
                                main_models = main_models, loss_metric = "MSE",
                                alternative = dm_alternative)
  
  dm_mse_better <- NA; dm_mse_worse <- NA
  if (!is.null(dm_mse)) {
    st_ben_mse    <- subset(dm_mse, base_model == summary_benchmark & main_model %in% main_models)
    dm_mse_better <- mean(st_ben_mse$p_value < mcs_alpha & st_ben_mse$statistic < 0, na.rm = TRUE)
    dm_mse_worse  <- mean(st_ben_mse$p_value < mcs_alpha & st_ben_mse$statistic > 0, na.rm = TRUE)
  }
  
  mcs_qlike <- run_mcs_for_all_stocks(results, model_names, alpha = mcs_alpha, B = mcs_B,
                                      statistic = mcs_statistic, loss = "QLIKE")
  dm_qlike  <- run_dm_for_results(results, h = config$h, model_names = model_names,
                                  main_models = main_models, loss_metric = "QLIKE",
                                  alternative = dm_alternative)
  
  dm_qlike_better <- NA; dm_qlike_worse <- NA
  if (!is.null(dm_qlike)) {
    st_ben_qlike    <- subset(dm_qlike, base_model == summary_benchmark & main_model %in% main_models)
    dm_qlike_better <- mean(st_ben_qlike$p_value < mcs_alpha & st_ben_qlike$statistic < 0, na.rm = TRUE)
    dm_qlike_worse  <- mean(st_ben_qlike$p_value < mcs_alpha & st_ben_qlike$statistic > 0, na.rm = TRUE)
  }
  
  dm_combined <- rbind(dm_mse, dm_qlike)
  summary_row <- data.frame(
    label             = config$label,    lags   = config$lags, h      = config$h,
    window            = config$window,
    mcs_mse_tr_main   = mcs_mse$perc[main_models[1]],
    mcs_mse_tr_bench  = mcs_mse$perc[summary_benchmark],
    dm_mse_better_bench  = dm_mse_better,  dm_mse_worse_bench  = dm_mse_worse,
    mcs_qlike_tr_main    = mcs_qlike$perc[main_models[1]],
    mcs_qlike_tr_bench   = mcs_qlike$perc[summary_benchmark],
    dm_qlike_better_bench = dm_qlike_better, dm_qlike_worse_bench = dm_qlike_worse,
    stringsAsFactors = FALSE
  )
  list(mcs_mse = mcs_mse, mcs_qlike = mcs_qlike, dm_results = dm_combined, summary_row = summary_row)
}

parse_config_from_filename <- function(fname) {
  lags   <- as.numeric(sub(".*lag([0-9]+).*",    "\\1", fname))
  h      <- as.numeric(sub(".*_h([0-9]+)_.*",    "\\1", fname))
  window <- as.numeric(sub(".*window([0-9]+).*", "\\1", fname))
  list(lags = lags, h = h, window = window, label = sub("\\.RData$", "", basename(fname)))
}

run_batch <- function(files, save_dir, model_names, main_models, summary_benchmark,
                      mcs_B, mcs_alpha, mcs_statistic, dm_alternative) {
  
  all_summaries <- list()
  files <- sort(files)
  
  for (k in seq_along(files)) {
    f <- files[k]
    set.seed(271 + k)
    cat("  -> Calculating tests for:", basename(f), "\n")
    cfg <- parse_config_from_filename(basename(f))
    env <- new.env(); load(f, envir = env)
    if (!"results" %in% ls(env)) next
    
    out <- run_all_for_results(
      env$results, cfg, model_names, main_models, summary_benchmark,
      mcs_B, mcs_alpha, mcs_statistic, dm_alternative
    )
    save(out, file = file.path(save_dir, paste0(cfg$label, "_tests.RData")))
    all_summaries[[k]] <- out$summary_row
  }
  if (length(all_summaries) == 0) return(NULL)
  summary_df <- do.call(rbind, all_summaries)
  write.csv(summary_df, file = file.path(save_dir, "batch_summary.csv"), row.names = FALSE)
  return(summary_df)
}


# LATEX PRESENTATION FUNCTIONS
extract_one_metric <- function(out, loss_key, alpha, model_order) {
  cfg     <- out$summary_row %>% dplyr::select(lags, window, h) %>%
    dplyr::rename(Lags = lags, Window = window, Horizon = h)
  mcs_obj <- if (loss_key == "MSE") out$mcs_mse else out$mcs_qlike
  mcs_vec <- mcs_obj$perc
  available_models <- intersect(model_order, names(mcs_vec))
  mcs_vec <- mcs_vec[available_models]
  
  mcs_row       <- as.data.frame(t(mcs_vec / 100))
  names(mcs_row) <- paste0("MCS__", names(mcs_row))
  
  dm     <- out$dm_results
  dm_row <- data.frame()
  
  if (!is.null(dm) && nrow(dm) > 0) {
    dm_sub <- dm %>% dplyr::filter(loss == loss_key, main_model == "predictions_startime")
    if (nrow(dm_sub) > 0) {
      dm_wide_df <- dm_sub %>%
        dplyr::group_by(base_model) %>%
        dplyr::summarise(
          win  = mean(p_value < alpha & statistic < 0, na.rm = TRUE),
          loss = mean(p_value < alpha & statistic > 0, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(combined = sprintf("%.3f / %.3f", win, loss)) %>%
        dplyr::select(base_model, combined) %>%
        tidyr::pivot_wider(names_from = base_model, values_from = combined) %>%
        as.data.frame()
      names(dm_wide_df) <- paste0("DM__", names(dm_wide_df))
      dm_row <- dm_wide_df
    }
  }
  if (ncol(dm_row) > 0) cbind(cfg, mcs_row, dm_row) else cbind(cfg, mcs_row)
}

make_latex_table <- function(data, prefix, caption, label, cols, model_display, fmt, footnote_text) {
  cols_needed <- c("Lags", "Window", "Horizon")
  tab <- data %>% dplyr::select(dplyr::all_of(cols_needed), dplyr::starts_with(prefix))
  
  for (m in cols) {
    col_code    <- paste0(prefix, m)
    display     <- model_display[[m]]
    if (col_code %in% names(tab)) {
      val          <- tab[[col_code]]
      tab[[display]] <- if (prefix == "MCS__" && is.numeric(val)) fmt(val) else ifelse(is.na(val), "--", val)
    } else {
      tab[[display]] <- "--"
    }
  }
  
  final_cols <- c("Lags", "Window", "Horizon", unname(model_display[cols]))
  tab        <- tab %>% dplyr::select(dplyr::all_of(final_cols))
  align_vec  <- c("c", "c", "c", rep("c", length(cols)))
  options(knitr.table.format = "latex")
  
  k <- kableExtra::kbl(tab, format = "latex", booktabs = TRUE, row.names = FALSE,
                       linesep = "", caption = caption, label = label, align = align_vec) %>%
    kableExtra::kable_styling(latex_options = c("hold_position")) %>%
    kableExtra::add_header_above(c("Configuration" = 3, "Models" = length(cols))) %>%
    kableExtra::footnote(general = footnote_text, threeparttable = TRUE)
  
  k_str <- as.character(k)
  k_str <- gsub("\\\\begin\\{tabular\\}", "\\\\resizebox{\\\\textwidth}{!}{%\n\\\\begin{tabular}", k_str)
  k_str <- gsub("\\\\end\\{tabular\\}",   "\\\\end{tabular}}",  k_str)
  return(k_str)
}

print_metric_tables <- function(loss_key, display_name, outs, alpha, model_order, model_display, fmt) {
  df_list <- lapply(outs, extract_one_metric, loss_key = loss_key, alpha = alpha, model_order = model_order)
  df <- dplyr::bind_rows(df_list) %>%
    dplyr::mutate(Lags = as.integer(Lags), Window = as.integer(Window), Horizon = as.integer(Horizon)) %>%
    dplyr::arrange(Lags, Window, Horizon)
  
  cat(paste0("\n\n%%% MCS TABLE (", display_name, ") %%%\n"))
  mcs_foot <- "Values represent the proportion of stocks included in the 95\\% MCS."
  cat(make_latex_table(df, "MCS__", paste0("Model Confidence Set Inclusion Rates (", display_name, ")"),
                       paste0("tab:mcs_", tolower(loss_key)), model_order, model_display, fmt, mcs_foot))
  
  cat(paste0("\n\n%%% DM TABLE (", display_name, ") %%%\n"))
  dm_foot        <- "Win Rate / Loss Rate. Win = StarTime significantly better; Loss = StarTime significantly worse."
  dm_competitors <- setdiff(model_order, "predictions_startime")
  cat(make_latex_table(df, "DM__", paste0("Diebold-Mariano Win / Loss Rates: StarTime vs Competitors (", display_name, ")"),
                       paste0("tab:dm_", tolower(loss_key)), dm_competitors, model_display, fmt, dm_foot))
}


# EXECUTION LOGIC
# --- PHASE A: Run statistical tests and save results to disk ---
if (RUN_STATISTICAL_TESTS) {
  message("-> Phase A: Running DM/MCS tests and saving to disk...")
  
  if (REPLICATE_PAPER_RESULTS) {
    target_files <- list.files(in_folder, pattern = "^Financial_lag.*\\.RData$", full.names = TRUE)
  } else {
    grid_df <- expand.grid(lag = FIN_MAX_LAG, h = FIN_H, window = FIN_WINDOW_SIZE, stringsAsFactors = FALSE)
    fname   <- sprintf("Financial_lag%d_h%d_window%d.RData", grid_df$lag, grid_df$h, grid_df$window)
    target_files <- file.path(in_folder, fname)
  }
  
  target_files <- target_files[file.exists(target_files)]
  target_files <- sort(target_files)
  
  if (length(target_files) == 0) {
    message("-> SKIPPING TESTS: No result .RData files found in: ", in_folder)
  } else {
    run_batch(
      files             = target_files,
      save_dir          = tests_dir,
      model_names       = c("predictions_startime", "predictions_simple_startime",
                            "predictions_har", "predictions_ar1", "predictions_rw"),
      main_models       = "predictions_startime",
      summary_benchmark = "predictions_har",
      mcs_B             = 5000,
      mcs_alpha         = 0.05,
      mcs_statistic     = "TR",
      dm_alternative    = "two.sided"
    )
    message("-> Tests saved to: ", tests_dir)
  }
}

# --- PHASE B: Load saved test results and print LaTeX tables ---
if (GENERATE_TABLES) {
  message("-> Phase B: Printing LaTeX tables from saved test results...")
  
  test_files <- list.files(tests_dir, pattern = "_tests\\.RData$", full.names = TRUE)
  test_files <- sort(test_files)
  
  if (length(test_files) == 0) {
    message("-> SKIPPING TABLE PRINT: No test files found in: ", tests_dir)
    message("   Set RUN_STATISTICAL_TESTS <- TRUE and re-run to generate them.")
  } else {
    outs <- lapply(test_files, function(f) { env <- new.env(); load(f, envir = env); env$out })
    outs <- Filter(Negate(is.null), outs)
    
    if (length(outs) > 0) {
      model_order   <- c("predictions_rw", "predictions_ar1", "predictions_har",
                         "predictions_startime", "predictions_simple_startime")
      model_display <- c(
        predictions_rw             = "Random Walk",
        predictions_ar1            = "AR1",
        predictions_har            = "HAR",
        predictions_startime       = "Post StarTime",
        predictions_simple_startime = "Simple StarTime"
      )
      fmt <- function(x, digits = 3) ifelse(is.na(x), "--", sprintf(paste0("%.", digits, "f"), x))
      
      print_metric_tables("MSE",   "MSFE",  outs, 0.05, model_order, model_display, fmt)
      print_metric_tables("QLIKE", "QLIKE", outs, 0.05, model_order, model_display, fmt)
    }
  }
}

message("Financial Table Generation Complete!")
message("=========================================================")