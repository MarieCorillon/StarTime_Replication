
load("data/aligned_macro_data.RData")

# --- Tree Structure Parameters (unchanged) ---
daily_tree <- c(60,12,3,1)
weekly_tree <- c(12,3,1)
monthly_tree <- c(3,1)
quarterly_tree <- c(1)

daily_sizes <- c(5,4,3,1)
weekly_sizes <- c(4,3,1)
monthly_sizes <- c(3,1)
quarterly_sizes <- c(1)

# ====================================================================
# Add quarter/year columns to the RAW data
# ====================================================================

aligned_macro_data$daily_final$quarters <- quarters(aligned_macro_data$daily_final$date)
aligned_macro_data$daily_final$years <- lubridate::year(aligned_macro_data$daily_final$date)

week_start <- as.Date(aligned_macro_data$weekly_final$week)
week_end <- ceiling_date(week_start, unit = "week")
aligned_macro_data$weekly_final$quarters <- quarters(week_end)
aligned_macro_data$weekly_final$years <- lubridate::year(aligned_macro_data$weekly_final$week)

aligned_macro_data$monthly_final$quarters <- quarters(aligned_macro_data$monthly_final$date)
aligned_macro_data$monthly_final$years <- lubridate::year(aligned_macro_data$monthly_final$date)

aligned_macro_data$quarterly_final$quarters <- quarters(aligned_macro_data$quarterly_final$date)
aligned_macro_data$quarterly_final$years <- lubridate::year(aligned_macro_data$quarterly_final$date)

# --- Create y data_frame ---
y_base <- aligned_macro_data$quarterly_final %>%
  arrange(date) %>%
  transmute(
    date_q   = as.Date(date),
    year_q   = years,
    quarter_q= quarters,
    GDP 
  )

add_q_index <- function(df) {
  df %>%
    mutate(
      q_index = (year_q - min(year_q)) * 4L +
        as.integer(substr(quarter_q, 2, 2))
    )
}

y_now <- y_base %>% add_q_index()

y_lag <- y_base %>%
  add_q_index() %>%
  transmute(
    q_index_tm1 = q_index + 1L,
    year_tm1    = year_q,
    quarter_tm1 = quarter_q
  )

y_df <- y_now %>%
  left_join(y_lag, by = c("q_index" = "q_index_tm1")) %>%
  select(-q_index)

# ====================================================================
# DAILY LAGS (from RAW data)
# ====================================================================
d_daily <- aligned_macro_data$daily_final

vars_to_lag <- d_daily %>%
  select(-date, -quarters, -years) %>%
  select(where(is.numeric)) %>%
  names()

L <- daily_tree[1] - 1

for (var in vars_to_lag) {
  for (k in 0:L) {
    lag_name <- paste0(var, "_lag", k)
    d_daily[[lag_name]] <- dplyr::lag(d_daily[[var]], k)  
  }
}

daily_q_lags <- d_daily %>%
  group_by(years, quarters) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup()

y_with_daily <- y_df %>%
  left_join(daily_q_lags, by = c("year_tm1" = "years", "quarter_tm1" = "quarters"))

# ====================================================================
# WEEKLY LAGS (from RAW data)
# ====================================================================
w_weekly <- aligned_macro_data$weekly_final

vars_w_to_lag <- w_weekly %>%
  select(-week, -quarters, -years) %>%
  select(where(is.numeric)) %>%
  names()

Lw <- weekly_tree[1] - 1

for (var in vars_w_to_lag) {
  for (k in 0:Lw) {
    lag_name <- paste0(var, "_wlag", k)
    w_weekly[[lag_name]] <- dplyr::lag(w_weekly[[var]], k)  
  }
}

weekly_q_lags <- w_weekly %>%
  group_by(years, quarters) %>%
  slice_max(week, n = 1, with_ties = FALSE) %>%
  ungroup()

y_with_daily_weekly <- y_with_daily %>%
  left_join(weekly_q_lags, by = c("year_tm1" = "years", "quarter_tm1" = "quarters"))

# ====================================================================
# MONTHLY LAGS 
# ====================================================================
m_monthly <- aligned_macro_data$monthly_final

vars_m_to_lag <- m_monthly %>%
  select(-date, -quarters, -years) %>%
  select(where(is.numeric)) %>%
  names()

Lm <- monthly_tree[1] - 1

for (var in vars_m_to_lag) {
  for (k in 0:Lm) {
    lag_name <- paste0(var, "_mlag", k)
    m_monthly[[lag_name]] <- dplyr::lag(m_monthly[[var]], k) 
  }
}

monthly_q_lags <- m_monthly %>%
  group_by(years, quarters) %>%
  slice_max(date, n = 1, with_ties = FALSE) %>%
  ungroup()

y_with_daily_weekly_monthly <- y_with_daily_weekly %>%
  left_join(monthly_q_lags, by = c("year_tm1" = "years", "quarter_tm1" = "quarters"))

# ====================================================================
# QUARTERLY LAGS (from RAW data)
# ====================================================================
q_quarterly <- aligned_macro_data$quarterly_final

vars_q_to_lag <- q_quarterly %>%
  select(-date, -quarters, -years) %>%
  select(where(is.numeric)) %>%
  names()

Lq <- quarterly_tree[1] - 1

for (var in vars_q_to_lag) {
  for (k in 0:Lq) {
    lag_name <- paste0(var, "_qlag", k)
    q_quarterly[[lag_name]] <- dplyr::lag(q_quarterly[[var]], k) 
  }
}

quarterly_q_lags <- q_quarterly %>%
  select(years, quarters, matches("_qlag[0-9]+$"))

y_full <- y_with_daily_weekly_monthly %>%
  left_join(quarterly_q_lags, by = c("year_tm1" = "years", "quarter_tm1" = "quarters"))

# ====================================================================
# Final X and y (RAW, UNSCALED)
# ====================================================================
clean_df <- y_full %>%
  tidyr::drop_na()

y_raw <- as.matrix(clean_df$GDP)  

X_raw_df <- clean_df %>%
  select(
    matches("_lag[0-9]+$"),
    matches("_wlag[0-9]+$"),
    matches("_mlag[0-9]+$"),
    matches("_qlag[0-9]+$")
  )

X_raw <- as.matrix(X_raw_df)  

# ====================================================================
# Create Reduced Dataset (if needed)
# ====================================================================
create_reduced_dataset <- function(X_full, selected_regressors) {
  lag_patterns <- list(
    weekly = c("wlag0", "wlag1", "wlag2", "wlag3", "wlag4", "wlag5", 
               "wlag6", "wlag7", "wlag8", "wlag9", "wlag10", "wlag11"),
    monthly = c("mlag0", "mlag1", "mlag2"),
    quarterly = c("qlag0")
  )
  
  code_to_pattern <- list(
    NFCI = "NFCI", ICSA = "ICSA", AMTMNO = "AMTMNO",
    CPIAUCSL = "CPIAUCSL", RSAFS = "RSAFS",
    HOUST = "HOUST", PERMIT = "PERMIT", UNRATE = "UNRATE",
    PAYEMS = "PAYEMS", INDPRO = "INDPRO", GDP = "GDP"
  )
  
  keep_cols <- character(0)
  
  for(regressor_code in selected_regressors) {
    pattern <- code_to_pattern[[regressor_code]]
    if(is.null(pattern)) next
    
    if(regressor_code == "GDP") {
      lags <- lag_patterns$quarterly
    } else if(regressor_code %in% c("NFCI", "ICSA")) {
      lags <- lag_patterns$weekly
    } else {
      lags <- lag_patterns$monthly
    }
    
    for(lag_suffix in lags) {
      col_pattern <- paste0("^", pattern, "_", lag_suffix, "$")
      matching_cols <- grep(col_pattern, colnames(X_full), value = TRUE)
      keep_cols <- c(keep_cols, matching_cols)
    }
  }
  
  X_full[, keep_cols, drop = FALSE]
}

selected_codes <- c("NFCI", "ICSA", "AMTMNO", "CPIAUCSL", "RSAFS", 
                    "HOUST", "PERMIT", "UNRATE", "PAYEMS", "INDPRO", "GDP")

X_raw_reduced <- create_reduced_dataset(X_raw, selected_codes)


# ====================================================================
# NOWCAST DATASET (h = 0): predict the SAME y_raw (GDP_t) as forecasts
# - X_raw/y_raw remain untouched (master)
# - X_nowcast is aligned to the SAME target quarters as clean_df
# - Missing nowcast features are forward-filled (no look-ahead)
# ====================================================================

# Master target index
master_index <- clean_df %>%
  dplyr::mutate(.row_id = dplyr::row_number()) %>%
  dplyr::select(.row_id, date_q, year_q, quarter_q, year_tm1, quarter_tm1)

# Build contemporaneous HF quarterly snapshots (quarter t)
# daily_q_lags / weekly_q_lags / monthly_q_lags already exist in your script
daily_q_now  <- daily_q_lags  %>% dplyr::rename(year_q = years, quarter_q = quarters)
weekly_q_now <- weekly_q_lags %>% dplyr::rename(year_q = years, quarter_q = quarters)
monthly_q_now<- monthly_q_lags%>% dplyr::rename(year_q = years, quarter_q = quarters)

nowcast_full <- y_df %>%
  dplyr::left_join(daily_q_now,   by = c("year_q", "quarter_q")) %>%
  dplyr::left_join(weekly_q_now,  by = c("year_q", "quarter_q")) %>%
  dplyr::left_join(monthly_q_now, by = c("year_q", "quarter_q")) %>%
  dplyr::left_join(quarterly_q_lags, by = c("year_tm1" = "years", "quarter_tm1" = "quarters"))


# Force perfect target alignment: keep exactly the same quarters/rows as clean_df
cleannowdf <- master_index %>%
  dplyr::left_join(
    nowcast_full,
    by = c("date_q", "year_q", "quarter_q", "year_tm1", "quarter_tm1")
  ) %>%
  dplyr::arrange(.row_id)

# Define nowcast target
y_nowcast <- y_raw

# Build X_nowcast with EXACT SAME columns/order as X_raw 
required_cols <- colnames(X_raw)
missing_required <- setdiff(required_cols, colnames(cleannowdf))
if (length(missing_required) > 0) {
  cat("Nowcast: adding missing columns (will be imputed):\n")
  print(missing_required)
  for (cc in missing_required) cleannowdf[[cc]] <- NA_real_
}

# Keep only the forecasting feature columns, in the same order
X_nowcast_df <- cleannowdf %>%
  dplyr::select(dplyr::all_of(required_cols))

# Identify which rows have at least one NA
na_row_indices <- which(rowSums(is.na(X_nowcast_df)) > 0)

# Print a diagnostic message to see exactly which rows are affected
cat("Rows with NA values after forward-fill:", paste(na_row_indices, collapse = ", "), "\n")

X_nowcast_df <- X_nowcast_df[-na_row_indices,]

# Check if there are any NAs remaining
which(rowSums(is.na(X_nowcast_df)) > 0)

X_nowcast <- as.matrix(X_nowcast_df)
cat("Final X_nowcast dimensions:", dim(X_nowcast)[1], "x", dim(X_nowcast)[2], "\n")

X_nowcast <- as.matrix(X_nowcast_df)

y_nowcast <- y_nowcast[-na_row_indices]

# Reduced nowcast dataset
if (exists("selected_codes")) {
  X_nowcast_reduced <- create_reduced_dataset(X_nowcast, selected_codes)
}



save(X_raw, file = "4_applications/Macro/Data/X_macro_raw.RData")
save(X_raw_reduced, file = "4_applications/Macro/Data/X_macro_raw_reduced.RData")
save(y_raw, file = "4_applications/Macro/Data/y_macro_raw.RData")

save(X_nowcast, file = "4_applications/Macro/Data/X_macro_nowcast.RData")
save(X_nowcast_reduced, file = "4_applications/Macro/Data/X_macro_nowcast_reduced.RData")
save(y_nowcast, file = "4_applications/Macro/Data/y_macro_nowcast.RData")
