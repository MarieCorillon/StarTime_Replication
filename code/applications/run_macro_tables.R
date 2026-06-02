message("=========================================================")
message("Starting Macro Table Generation...")

in_folder  <- file.path(RESULTS_DIR, "applications/macro")
tests_dir  <- file.path(RESULTS_DIR, "applications/macro/tests")
tables_dir <- file.path(FIGURES_DIR, "applications/macro/tables")

if (!dir.exists(tests_dir))  dir.create(tests_dir,  recursive = TRUE)
if (!dir.exists(tables_dir)) dir.create(tables_dir, recursive = TRUE)

# 1. DATA EXTRACTION FUNCTION
total_obs      <- 131
global_is_peak <- rep(FALSE, total_obs)
global_is_peak[111:118] <- TRUE   # COVID-19 shock indices

extract_metrics <- function(filepath, is_forecast, window_size) {
  env <- new.env()
  load(filepath, envir = env)
  res <- env$res
  
  y_true     <- if (!is.null(res$penalized$evaluation$y_true)) res$penalized$evaluation$y_true else res$random_walk$evaluation$y_true
  pred_post  <- res$penalized$evaluation$pred_star_post
  pred_simple <- res$penalized$results$simple_predictions
  pred_midas  <- res$penalized$results$midasml_predictions
  
  if (is_forecast) {
    pred_ar1 <- res$ar1$evaluation$pred
    pred_rw  <- res$random_walk$evaluation$pred
    y_true     <- y_true[-length(y_true)]
    pred_post   <- pred_post[-length(pred_post)]
    pred_simple <- pred_simple[-length(pred_simple)]
    pred_midas  <- pred_midas[-length(pred_midas)]
    pred_ar1    <- pred_ar1[-length(pred_ar1)]
    pred_rw     <- pred_rw[-length(pred_rw)]
  } else {
    pred_ar1 <- rep(NA, length(y_true))
    pred_rw  <- rep(NA, length(y_true))
  }
  
  global_idx <- if (window_size == 66) 67:131 else 106:131
  is_peak    <- global_is_peak[global_idx]
  
  calc_msfe <- function(p, y) {
    if (all(is.na(p))) return(c(TOTAL = NA, OFF_PEAK = NA, PEAK = NA))
    c(TOTAL    = mean((p - y)^2,               na.rm = TRUE) * 1000,
      OFF_PEAK = mean((p[!is_peak] - y[!is_peak])^2, na.rm = TRUE) * 1000,
      PEAK     = mean((p[is_peak]  - y[is_peak])^2,  na.rm = TRUE) * 1000)
  }
  
  calc_se <- function(p, y) {
    if (all(is.na(p))) return(rep(NA, length(y)))
    ((p - y)^2) * 100
  }
  
  list(
    msfe = list(
      Post   = calc_msfe(pred_post,   y_true), Simple = calc_msfe(pred_simple, y_true),
      MIDAS  = calc_msfe(pred_midas,  y_true), AR1    = calc_msfe(pred_ar1,    y_true),
      RW     = calc_msfe(pred_rw,     y_true)
    ),
    se = data.frame(
      Global_Index = global_idx, Is_Peak = is_peak,
      Post   = calc_se(pred_post,   y_true), Simple = calc_se(pred_simple, y_true),
      MIDAS  = calc_se(pred_midas,  y_true), AR1    = calc_se(pred_ar1,    y_true),
      RW     = calc_se(pred_rw,     y_true)
    )
  )
}


# 2. MCS COMPUTATION FUNCTION (returns p-values only — no printing)
compute_mcs_pvals <- function(d, off_peak = FALSE, seed = 271) {
  b_mat <- function(s, suffix) {
    df  <- if (off_peak) s$se[!s$se$Is_Peak, ] else s$se
    mat <- df %>% dplyr::select(Post, Simple, MIDAS, AR1, RW)
    names(mat) <- paste0(suffix, names(mat))
    mat
  }
  
  mat_full_h0 <- b_mat(d$full_h0, "FULL_H0_") %>% dplyr::select(-FULL_H0_AR1, -FULL_H0_RW)
  mat_full_h1 <- b_mat(d$full_h1, "FULL_H1_")
  mat_red_h0  <- b_mat(d$red_h0,  "RED_H0_")  %>% dplyr::select(-RED_H0_AR1,  -RED_H0_RW)
  mat_red_h1  <- b_mat(d$red_h1,  "RED_H1_")  %>% dplyr::select(-RED_H1_AR1,  -RED_H1_RW)
  
  mcs_matrix <- as.matrix(cbind(mat_full_h0, mat_full_h1, mat_red_h0, mat_red_h1))
  
  set.seed(seed)
  
  capture.output(
    suppressWarnings(
      mcs_res <- MCSprocedure(mcs_matrix, alpha = 0.05, verbose = FALSE)
    )
  )
  
  res_df <- if (.hasSlot(mcs_res, "show")) mcs_res@show else as.data.frame(mcs_res)
  pvals  <- setNames(
    if ("MCS_M" %in% colnames(res_df)) res_df[, "MCS_M"] else res_df[, 1],
    rownames(res_df)
  )
  pvals  # named numeric vector: names are column names, values are p-values
}


# 3. LATEX PRINTING FUNCTIONS
print_msfe_latex <- function(w, d) {
  fmt <- function(x) sprintf("%.3f", x)
  cat(sprintf("\\begin{table}[h]
\\centering
\\caption{MSFE for window $W = %d$}
\\begin{threeparttable}
\\begin{tabular}{lcccccccc}
\\toprule
& \\multicolumn{2}{c}{Post StarTime}
& \\multicolumn{2}{c}{Simple StarTime}
& \\multicolumn{2}{c}{MIDAS-ML}
& AR(1) & Random walk \\\\
\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}
& $h=0$ & $h=1$ & $h=0$ & $h=1$ & $h=0$ & $h=1$ &  &  \\\\
\\midrule
\\addlinespace[0.3em]
\\multicolumn{9}{l}{\\textit{Reduced Set}} \\\\
", w))
  
  for (metric in c("TOTAL", "OFF_PEAK", "PEAK")) {
    label <- switch(metric,
                    TOTAL    = "\\hspace{1em}Total    ",
                    OFF_PEAK = "\\hspace{1em}Non-peak ",
                    PEAK     = "\\hspace{1em}COVID-Peak")
    cat(sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n",
                label,
                fmt(d$red_h0$msfe$Post[[metric]]),   fmt(d$red_h1$msfe$Post[[metric]]),
                fmt(d$red_h0$msfe$Simple[[metric]]),  fmt(d$red_h1$msfe$Simple[[metric]]),
                fmt(d$red_h0$msfe$MIDAS[[metric]]),   fmt(d$red_h1$msfe$MIDAS[[metric]]),
                fmt(d$full_h1$msfe$AR1[[metric]]),    fmt(d$full_h1$msfe$RW[[metric]])))
  }
  
  cat("\n\\addlinespace[0.5em]\n\\multicolumn{9}{l}{\\textit{Full Set}} \\\\\n")
  
  for (metric in c("TOTAL", "OFF_PEAK", "PEAK")) {
    label <- switch(metric,
                    TOTAL    = "\\hspace{1em}Total    ",
                    OFF_PEAK = "\\hspace{1em}Non-peak ",
                    PEAK     = "\\hspace{1em}COVID-Peak")
    cat(sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s \\\\\n",
                label,
                fmt(d$full_h0$msfe$Post[[metric]]),  fmt(d$full_h1$msfe$Post[[metric]]),
                fmt(d$full_h0$msfe$Simple[[metric]]), fmt(d$full_h1$msfe$Simple[[metric]]),
                fmt(d$full_h0$msfe$MIDAS[[metric]]),  fmt(d$full_h1$msfe$MIDAS[[metric]]),
                fmt(d$full_h1$msfe$AR1[[metric]]),    fmt(d$full_h1$msfe$RW[[metric]])))
  }
  
  cat(sprintf("\\bottomrule
\\end{tabular}
\\begin{tablenotes}
\\item \\textit{Note: MSFE values are multiplied by $1000$ due to the small scale of GDP growth.}
\\end{tablenotes}
\\end{threeparttable}
\\label{tab:macro_w%d}
\\end{table}\n\n", w))
}

print_mcs_latex <- function(w, pvals, off_peak = FALSE) {
  fmt        <- function(x) sprintf("%.3f", x)
  get_p      <- function(k) if (k %in% names(pvals)) fmt(as.numeric(pvals[k])) else "-"
  label_type <- if (off_peak) "Non-Peak" else "Total"
  label_key  <- if (off_peak) "nonpeak"  else "total"
  
  cat(sprintf("\\begin{table}[!h]
\\centering
\\caption{\\label{tab:mcs_macro_w%d_%s} Model Confidence Set $p$-values: $T_{max}$ (%s, $W=%d$)}
\\resizebox{\\textwidth}{!}{
\\begin{threeparttable}
\\begin{tabular}[t]{lccccccc}
\\toprule
\\multicolumn{3}{c}{Configuration} & \\multicolumn{5}{c}{Models} \\\\
\\cmidrule(l{3pt}r{3pt}){1-3} \\cmidrule(l{3pt}r{3pt}){4-8}
Reduction & Window & Horizon & Post StarTime & Simple StarTime & MIDAS-ML & AR(1) & Random Walk\\\\
\\midrule
\\addlinespace[0.3em]
\\multicolumn{8}{l}{}\\\\
\\hspace{1em}Full    & %d & 0 & %s & %s & %s & & \\\\
\\hspace{1em}Full    & %d & 1 & %s & %s & %s & & \\\\
\\hspace{1em}Reduced & %d & 0 & %s & %s & %s & & \\\\
\\hspace{1em}Reduced & %d & 1 & %s & %s & %s & & \\\\
\\hspace{0em}\\multicolumn{2}{l}{Univariate Benchmark} & & & & & %s & %s \\\\
\\bottomrule
\\end{tabular}
\\begin{tablenotes}
\\item \\textit{Note: Values represent $p$-values for the MCS using the $T_{max}$ statistic. $p$-values $\\ge 0.05$ indicate inclusion in the MCS.}
\\end{tablenotes}
\\end{threeparttable}}
\\end{table}\n\n",
              w, label_key, label_type, w,
              w, get_p("FULL_H0_Post"), get_p("FULL_H0_Simple"), get_p("FULL_H0_MIDAS"),
              w, get_p("FULL_H1_Post"), get_p("FULL_H1_Simple"), get_p("FULL_H1_MIDAS"),
              w, get_p("RED_H0_Post"),  get_p("RED_H0_Simple"),  get_p("RED_H0_MIDAS"),
              w, get_p("RED_H1_Post"),  get_p("RED_H1_Simple"),  get_p("RED_H1_MIDAS"),
              get_p("FULL_H1_AR1"),     get_p("FULL_H1_RW")))
}


# 4. EXECUTION LOGIC
data_store_file <- file.path(tests_dir, "macro_data_store.RData")

# --- PHASE A: Compute data_store + all MCS p-values, save everything to disk ---
if (RUN_STATISTICAL_TESTS) {
  message("-> Phase A: Extracting metrics and computing MCS p-values...")
  
  expected_files <- c(
    "nowcast_full_w66.RData",    "nowcast_reduced_w66.RData",
    "forecast_full_w66.RData",   "forecast_reduced_w66.RData",
    "nowcast_full_w105.RData",   "nowcast_reduced_w105.RData",
    "forecast_full_w105.RData",  "forecast_reduced_w105.RData"
  )
  full_paths  <- file.path(in_folder, expected_files)
  files_exist <- sapply(full_paths, file.exists)
  
  if (!all(files_exist)) {
    message("-> SKIPPING TESTS: Not all required macro result files are present.")
    message("   Missing: ", paste(expected_files[!files_exist], collapse = ", "))
  } else {
    # Build data store
    data_store <- list()
    for (w in c(66, 105)) {
      data_store[[paste0("w", w)]] <- list(
        full_h0 = extract_metrics(file.path(in_folder, sprintf("nowcast_full_w%d.RData",    w)), FALSE, w),
        full_h1 = extract_metrics(file.path(in_folder, sprintf("forecast_full_w%d.RData",   w)), TRUE,  w),
        red_h0  = extract_metrics(file.path(in_folder, sprintf("nowcast_reduced_w%d.RData",  w)), FALSE, w),
        red_h1  = extract_metrics(file.path(in_folder, sprintf("forecast_reduced_w%d.RData", w)), TRUE,  w)
      )
    }
    
    # Pre-compute all MCS p-values and store in a named list
    mcs_pvals <- list(
      w66_total    = compute_mcs_pvals(data_store$w66,  off_peak = FALSE, seed = 271),
      w66_nonpeak  = compute_mcs_pvals(data_store$w66,  off_peak = TRUE,  seed = 272),
      w105_total   = compute_mcs_pvals(data_store$w105, off_peak = FALSE, seed = 273),
      w105_nonpeak = compute_mcs_pvals(data_store$w105, off_peak = TRUE,  seed = 274)
    )
    
    save(data_store, mcs_pvals, file = data_store_file)
    message("-> All test results saved to: ", data_store_file)
  }
}

# --- PHASE B: Load saved results and print LaTeX tables (zero recomputation) ---
if (GENERATE_TABLES) {
  message("-> Phase B: Printing macro LaTeX tables from saved results...")
  
  if (!file.exists(data_store_file)) {
    message("-> SKIPPING TABLE PRINT: No saved data store found at: ", data_store_file)
    message("   Set RUN_STATISTICAL_TESTS <- TRUE and re-run to generate it.")
  } else {
    env <- new.env()
    load(data_store_file, envir = env)
    data_store <- env$data_store
    mcs_pvals  <- env$mcs_pvals
    
    cat("\n%%% ================= MSFE LATEX TABLES ================= %%%\n\n")
    print_msfe_latex(66,  data_store$w66)
    print_msfe_latex(105, data_store$w105)
    
    cat("\n%%% ================= MCS LATEX TABLES ================== %%%\n\n")
    print_mcs_latex(66,  mcs_pvals$w66_total,    off_peak = FALSE)
    print_mcs_latex(66,  mcs_pvals$w66_nonpeak,  off_peak = TRUE)
    print_mcs_latex(105, mcs_pvals$w105_total,   off_peak = FALSE)
    print_mcs_latex(105, mcs_pvals$w105_nonpeak, off_peak = TRUE)
  }
}

message("Macro Table Generation Complete!")
message("=========================================================")