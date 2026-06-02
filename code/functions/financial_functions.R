prepare_window_data <- function(log_rv, t_now, window_size, max_lag, h) {
  # t_now: The index of the most recent observation we have.
  # "Window Size" = 1000 means we look at the last 1000 raw observations.
  effective_window <- window_size+max_lag
  # 1. Identify Raw Data Range
  # Start: t_now - 1000 + 1. End: t_now.
  raw_start <- t_now - effective_window + 1
  raw_end   <- t_now
  
  # Check for start of series
  if (raw_start < 1) stop("Not enough data for this window size")
  
  # Extract Raw Data 
  raw_block <- as.numeric(log_rv[raw_start:raw_end])
  
  # 2. Standardize
  mu <- mean(raw_block)
  sigma <- sd(raw_block)
  if (sigma == 0) sigma <- 1
  scaled_block <- (raw_block - mu) / sigma
  
  # 3. Create Lags
  # If raw_block has 1000 points, embed creates (1000 - max_lag) rows.
  # This is the standard "Effective Sample Size" reduction.
  mat <- stats::embed(scaled_block, max_lag + 1)
  
  # INITIAL (Default) split: Gap is 1
  y_all <- mat[, 1]
  X_all <- mat[, -1, drop=FALSE]
  
  if (h > 1) {
    n <- length(y_all)
    
    # Starts at index h
    y_train <- y_all[h:n]
    
    # Row i of X_all is Time T. Row i of y_all is Time T+1.
    # We want Target T+h. That is Row i + (h-1) of y_all.
    
    X_startime_train <- X_all[1:(n - h + 1), , drop=FALSE]
  } else {
    y_train <- y_all
    X_startime_train <- X_all
  }
  
  # 4. Test Features
  # The features for predicting t_now + h are the lags ending at t_now.
  # These are simply the LAST 'max_lag' values of our scaled_block.
  scaled_test_features <- tail(scaled_block, max_lag)
  
  X_test_vector <- rev(scaled_test_features)
  X_test_matrix <- matrix(X_test_vector, nrow = 1)
  
  # 5. True Target (Evaluation)
  if ((t_now + h) <= length(log_rv)) {
    raw_true_target <- log_rv[t_now + h]
    y_true_scaled <- (raw_true_target - mu) / sigma
  } else {
    raw_true_target <- NA; y_true_scaled <- NA
  }
  
  # 6. HAR Features
  if (max_lag == 20 || max_lag == 40) {
    rv_lag1 <- X_startime_train[, 1]
    rv_week <- rowMeans(X_startime_train[, 1:5, drop=FALSE])
    rv_month <- rowMeans(X_startime_train[, 1:20, drop=FALSE])
    X_har_train <- cbind(rv_lag1, rv_week, rv_month)
    
    rv_lag1_test <- X_test_matrix[, 1]
    rv_week_test <- mean(X_test_matrix[, 1:5])
    rv_month_test <- mean(X_test_matrix[, 1:20])
    X_har_test <- matrix(c(rv_lag1_test, rv_week_test, rv_month_test), nrow = 1)
    colnames(X_har_test) <- c("rv_lag1", "rv_week", "rv_month")
  } else {
    X_har_train <- NULL; X_har_test <- NULL
  }
  
  list(y_train=y_train, X_startime_train=X_startime_train, X_har_train=X_har_train,
       X_test_startime=X_test_matrix, X_test_har=X_har_test,
       y_true_scaled=y_true_scaled, raw_true_target = raw_true_target,
       mu=mu, sigma=sigma)
}


ols_rolling_predict <- function(log_rv, window_size, h, max_lag, type = "HAR") {
  
  # --- Index Definition (Matching Startime) ---
  start_t <- window_size+max_lag
  end_t   <- length(log_rv) - h
  
  n_steps <- end_t - start_t + 1
  
  # Initialize vectors
  predictions <- vector("numeric", n_steps)         # De-standardized (Log-RV scale)
  predictions_scaled <- vector("numeric", n_steps)  # Standardized scale (for MSE)
  true_vals_scaled <- vector("numeric", n_steps)    # Standardized truth (for MSE)
  true_vals <- vector("numeric", n_steps)
  betas <- vector("list", n_steps)
  
  # Loop counter
  j <- 1
  
  for (t_now in start_t:end_t) {
    
    # 1. Prepare Data using t_now
    data <- prepare_window_data(log_rv, t_now, window_size, max_lag, h)
    
    # 2. Select Features based on Model Type
    if (type == "HAR") {
      X_train <- data$X_har_train
      X_test  <- data$X_test_har
    } else if (type == "AR1") {
      # AR1 uses only the first column of HAR (Lag 1)
      X_train <- data$X_har_train[, 1, drop = FALSE]
      X_test  <- data$X_test_har[, 1, drop = FALSE]
    } else {
      stop("Unknown type for ols_rolling_predict. Use 'HAR' or 'AR1'.")
    }
    
    # 3. Estimate (OLS on Standardized Data)
    beta_hat <- ols(X = X_train, y = data$y_train)
    
    # 4. Predict (Result is Scaled)
    pred_scaled_val <- as.numeric(X_test %*% beta_hat)
    
    # 5. Store Results
    # A. Scaled Prediction (for MSE comparison with StarTime)
    predictions_scaled[j] <- pred_scaled_val
    
    # B. De-standardized Prediction (Log-RV scale) - For final output/QLIKE input
    predictions[j] <- pred_scaled_val * data$sigma + data$mu
    
    # C. Scaled Truth (for MSE)
    true_vals_scaled[j] <- data$y_true_scaled
    
    # D. Truth
    true_vals[j] <- data$raw_true_target
    
    betas[[j]] <- beta_hat
    
    j <- j + 1
  }
  
  return(list(
    predictions = predictions,           # Log-RV scale
    predictions_scaled = predictions_scaled, # Standardized scale
    true_vals_scaled = true_vals_scaled, # Standardized Truth, 
    true_vals = true_vals,
    betas = betas
  ))
}




startime_rolling_predict <- function(log_rv, max_lag, window_size, h, up,
                                     P_list, Ki_list, rho,
                                     K_in, stop_thresh, ic_type, thresh,
                                     nlambda) {
  
  A <- create_matrix_A(P_list, Ki_list)
  sp_Astar <- create_matrix_A_star(P_list, Ki_list)
  
  # --- Index Definition ---
  # First t_now = window_size
  # Last t_now  = length(log_rv) - h (because we need target at t+h)
  
  start_t <- window_size+max_lag
  end_t   <- length(log_rv) - h
  
  n_steps <- end_t - start_t + 1
  
  # Initialize vectors
  predictions <- vector("numeric", n_steps)
  simple_predictions <- vector("numeric", n_steps)
  predictions_scaled <- vector("numeric", n_steps)
  simple_predictions_scaled <- vector("numeric", n_steps)
  true_vals_scaled <- vector("numeric", n_steps) 
  true_vals <- vector("numeric", n_steps)
  
  betas <- vector("list", n_steps)
  simple_betas <- vector("list", n_steps)
  lambda1_vals <- vector("numeric", n_steps)
  lambda2_vals <- vector("numeric", n_steps)
  simple_lambda1_vals <- vector("numeric", n_steps)
  simple_lambda2_vals <- vector("numeric", n_steps)
  
  last_lambda1 <- NA; last_lambda2 <- NA
  simple_last_lambda1 <- NA; simple_last_lambda2 <- NA 
  
  # Loop counter j for storage
  j <- 1
  
  for (t_now in start_t:end_t) {
    # 1. Prepare Data using t_now (Time-based index)
    data <- prepare_window_data(log_rv, t_now, window_size, max_lag, h)
    
    X_train <- data$X_startime_train
    y_train <- data$y_train
    X_test  <- data$X_test_startime
    
    # 2. Hyperparameter Tuning (every 'up' steps)
    if ((j - 1) %% up == 0) {
      l1_grid <- lambda_grid(X_train, y_train, nlambda, "node sparsity", A, sp_Astar, NULL)
      l2_grid <- lambda_grid(X_train, y_train, nlambda, "leaf sparsity", A, sp_Astar, NULL)
      
      # Full
      est <- startime_ic(X_train, y_train, P_list, Ki_list, post=TRUE,
                         lambda1_grid = l1_grid, lambda2_grid = l2_grid,
                         rho = rho, K_in = K_in, stop_thresh = stop_thresh, thresh = thresh)
      last_lambda1 <- est$lambda1_best
      last_lambda2 <- est$lambda2_best
      
      # Simple
      est_simp <- startime_ic(X_train, y_train, P_list, Ki_list, post=FALSE,
                              lambda1_grid = l1_grid, lambda2_grid = l2_grid,
                              rho = rho, K_in = K_in, stop_thresh = stop_thresh, thresh = thresh)
      simple_last_lambda1 <- est_simp$lambda1_best
      simple_last_lambda2 <- est_simp$lambda2_best
    }
    
    # 3. Estimation
    est_final <- startime(X_train, y_train, P_list, Ki_list, post=TRUE,
                          lambda1 = last_lambda1, lambda2 = last_lambda2,
                          rho = rho, K_in = K_in, stop_thresh = stop_thresh)
    
    est_simp_final <- startime(X_train, y_train, P_list, Ki_list, post=FALSE,
                               lambda1 = simple_last_lambda1, lambda2 = simple_last_lambda2,
                               rho = rho, K_in = K_in, stop_thresh = stop_thresh)
    
    # 4. Predict (Result is Scaled)
    pred_scaled <- X_test %*% est_final$beta
    simp_pred_scaled <- X_test %*% est_simp_final$beta
    
    # 5. De-standardize and Store
    predictions_scaled[j] <- pred_scaled
    simple_predictions_scaled[j] <- simp_pred_scaled
    predictions[j] <- pred_scaled * data$sigma + data$mu
    simple_predictions[j] <- simp_pred_scaled * data$sigma + data$mu
    true_vals_scaled[j] <- data$y_true_scaled 
    true_vals[j] <- data$raw_true_target
    
    # Save Metadata
    betas[[j]] <- est_final$beta
    simple_betas[[j]] <- est_simp_final$beta
    lambda1_vals[j] <- last_lambda1
    lambda2_vals[j] <- last_lambda2
    simple_lambda1_vals[j] <- simple_last_lambda1
    simple_lambda2_vals[j] <- simple_last_lambda2
    
    j <- j + 1
  }
  
  return(list(
    predictions = predictions, 
    simple_predictions = simple_predictions,
    predictions_scaled = predictions_scaled, 
    simple_predictions_scaled = simple_predictions_scaled,
    betas = betas, 
    simple_betas = simple_betas,
    lambda1 = lambda1_vals, 
    lambda2 = lambda2_vals,
    simple_lambda1 = simple_lambda1_vals, 
    simple_lambda2 = simple_lambda2_vals,
    true_vals_scaled = true_vals_scaled,
    true_vals = true_vals
  ))
}



random_walk_rolling_predict <- function(log_rv, window_size, h, max_lag) {
  
  # Index Definition
  start_t <- window_size+max_lag
  end_t   <- length(log_rv) - h
  
  n_steps <- end_t - start_t + 1
  
  # Initialize vectors
  predictions <- vector("numeric", n_steps)        # Log-RV scale
  predictions_scaled <- vector("numeric", n_steps) # Standardized scale (for MSE)
  true_vals_scaled <- vector("numeric", n_steps)   # Standardized Truth
  true_vals <- vector("numeric", n_steps)   # Truth
  
  j <- 1
  
  for (t_now in start_t:end_t) {
    
    # 1. Get Window Stats (mu, sigma) to standardize
    data <- prepare_window_data(log_rv, t_now, window_size, max_lag, h)
    
    # 2. Random Walk Logic
    # The forecast for t+h is simply the value at t_now.
    rw_forecast_raw <- log_rv[t_now]
    
    # 3. Standardize the Forecast
    # (RW_raw - mu) / sigma
    rw_forecast_scaled <- (rw_forecast_raw - data$mu) / data$sigma
    
    # 4. Store Results
    predictions[j] <- rw_forecast_raw
    predictions_scaled[j] <- rw_forecast_scaled
    true_vals_scaled[j] <- data$y_true_scaled
    true_vals[j] <- data$raw_true_target
    
    j <- j + 1
  }
  
  return(list(
    predictions = predictions,
    predictions_scaled = predictions_scaled,
    true_vals_scaled = true_vals_scaled, 
    true_vals = true_vals
  ))
}



qlike <- function(actual, predicted) {
  mean((exp(actual) / exp(predicted)) - (actual - predicted) - 1)
}




parallel_run_all_stocks <- function(variances_df, max_lag, window_size, h, up,
                                    P_list, Ki_list, rho, K_in, 
                                    stop_thresh = 1e-05, ic_type = "BIC", 
                                    thresh = 1, nlambda = 10) {
  
  stock_names <- colnames(variances_df)[-1]
  
  results <- foreach(
    stock = stock_names,
    .packages = c(
      "Rcpp",
      "RcppArmadillo",
      "glmnet",
      "doParallel",
      "doRNG",
      "foreach",
      "ggplot2",
      "dplyr",
      "tidyr",
      "StarTime"
    ),
    .export = c(
      "prepare_window_data",          
      "startime_rolling_predict", 
      "ols_rolling_predict",
      "random_walk_rolling_predict",
      "qlike",
      "compute_mse"
    )
  ) %dopar% {
    
    # 1. Get log-transformed realized volatility
    log_rv <- log(variances_df[[stock]])
    
    # 2. Run StarTime (with both post and simple versions)
    res_startime <- startime_rolling_predict(
      log_rv = log_rv,
      max_lag = max_lag,
      window_size = window_size,
      h = h,
      up = up,
      P_list = P_list,
      Ki_list = Ki_list,
      rho = rho,
      K_in = K_in,
      stop_thresh = stop_thresh,
      ic_type = ic_type,
      thresh = thresh,
      nlambda = nlambda
    )
    
    # 3. Run HAR (OLS with HAR features)
    res_har <- ols_rolling_predict(
      log_rv = log_rv,
      window_size = window_size,
      h = h,
      max_lag = max_lag,
      type = "HAR"
    )
    
    # 4. Run AR(1)
    res_ar1 <- ols_rolling_predict(
      log_rv = log_rv,
      window_size = window_size,
      h = h,
      max_lag = max_lag,
      type = "AR1"
    )
    
    # 5. Run Random Walk
    res_rw <- random_walk_rolling_predict(
      log_rv = log_rv,
      window_size = window_size,
      h = h,
      max_lag = max_lag
    )
    
    # 6. Extract True Values (all models return the same true_vals)
    true_values <- res_startime$true_vals  # Raw Log-RV scale
    
    # 7. Calculate MSE (on Log-RV scale)
    mse_startime <- compute_mse(true_values, res_startime$predictions)
    mse_simple_startime <- compute_mse(true_values, res_startime$simple_predictions)
    mse_HAR <- compute_mse(true_values, res_har$predictions)
    mse_ar1 <- compute_mse(true_values, res_ar1$predictions)
    mse_rw <- compute_mse(true_values, res_rw$predictions)
    
    # 8. Calculate QLIKE (on Realized Variance scale)
    qlike_startime <- qlike(exp(true_values), exp(res_startime$predictions))
    qlike_simple_startime <- qlike(exp(true_values), exp(res_startime$simple_predictions))
    qlike_HAR <- qlike(exp(true_values), exp(res_har$predictions))
    qlike_ar1 <- qlike(exp(true_values), exp(res_ar1$predictions))
    qlike_rw <- qlike(exp(true_values), exp(res_rw$predictions))
    
    # 9. Return Results
    list(
      # Predictions (Log-RV scale)
      predictions_startime = res_startime$predictions,
      predictions_simple_startime = res_startime$simple_predictions,
      predictions_har = res_har$predictions,
      predictions_ar1 = res_ar1$predictions,
      predictions_rw = res_rw$predictions,
      
      # True Values
      true_values = true_values,
      
      # MSE (Log-RV scale)
      mse_startime = mse_startime,
      mse_simple_startime = mse_simple_startime,
      mse_HAR = mse_HAR,
      mse_ar1 = mse_ar1,
      mse_rw = mse_rw,
      
      # QLIKE (RV scale)
      qlike_startime = qlike_startime,
      qlike_simple_startime = qlike_simple_startime,
      qlike_HAR = qlike_HAR,
      qlike_ar1 = qlike_ar1,
      qlike_rw = qlike_rw,
      
      # Coefficients (for interpretation/diagnostics)
      startime_coefficients = res_startime$betas,
      simple_startime_coefficients = res_startime$simple_betas,
      har_coefficients = res_har$betas
    )
  }
  
  names(results) <- stock_names
  return(results)
}