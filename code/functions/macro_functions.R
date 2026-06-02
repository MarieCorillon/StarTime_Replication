valid_lambdas <- function(ic_values, lambda1_grid, lambda2_grid){
  lambda1 <- rep(lambda1_grid, length(lambda1_grid))
  lambda2 <- rep(lambda2_grid, each=length(lambda2_grid))
  bic_df <- data.frame(cbind(ic_values, lambda1, lambda2))
  
  ## Find the smallest lambda1 lambda2
  # Convert all Inf values to NA so ggplot can handle them
  bic_df$ic_values[is.infinite(bic_df$ic_values)] <- NA
  # Filter for valid rows (where BIC is not NA)
  valid_df <- bic_df[!is.na(bic_df$ic_values), ]
  
  # Find the row with the smallest simultaneous lambdas that is valid
  min_corner <- valid_df[which.min(valid_df$lambda1 * valid_df$lambda2), ]
  valid_lambda1 <- min_corner$lambda1
  valid_lambda2 <- min_corner$lambda2
  
  return(list(valid_lambda1 = valid_lambda1, valid_lambda2=valid_lambda2))
}


penalized_rolling_predict <- function(X, y, window_size, h, up,
                                      P_list, Ki_list, rho,
                                      K_in, stop_thresh, ic_type, thresh,
                                      nlambda, b2_final) {
  
  A <- create_matrix_A(P_list, Ki_list)
  sp_Astar <- create_matrix_A_star(P_list, Ki_list)
  
  n <- nrow(X)
  n_steps<- n - window_size - h
  
  # Initialize vectors
  predictions <- vector("numeric", n_steps)
  betas <- vector("list", n_steps)
  lambda1_vals <- vector("numeric", n_steps)
  lambda2_vals <- vector("numeric", n_steps)
  time_star_post <- vector("numeric", n_steps)
  
  simple_predictions <- vector("numeric", n_steps)
  simple_betas <- vector("list", n_steps)
  simple_lambda1_vals <- vector("numeric", n_steps)
  simple_lambda2_vals <- vector("numeric", n_steps)
  time_simple <- vector("numeric", n_steps)
  
  midasml_predictions <- vector("numeric", n_steps)
  midasml_betas <- vector("list", n_steps)
  midasml_lambda_vals <- vector("numeric", n_steps)
  time_midasml <- vector("numeric", n_steps)
  
  last_lambda1 <- NA; last_lambda2 <- NA
  simple_last_lambda1 <- NA; simple_last_lambda2 <- NA
  midasml_last_lambda <- NA
  
  for (j in 1:n_steps) {
    # Forecast using point i (X_i)
    i <- window_size + j
    train_start_X <- j
    train_end_X <- j + window_size - h - 1
    
    X_train_raw <- X[train_start_X:train_end_X, , drop=FALSE]
    
    y_train_raw <- y[(train_start_X + h):(train_end_X + h)]
    
    X_test_raw <- X[i, , drop=FALSE]
    
    # STANDARDIZE INTERNALLY
    X_train_means <- colMeans(X_train_raw)
    X_train_sds   <- apply(X_train_raw, 2, sd)
    X_train <- sweep(X_train_raw, 2, X_train_means, "-")
    X_train <- sweep(X_train, 2, X_train_sds, "/")
    X_test  <- sweep(X_test_raw, 2, X_train_means, "-")
    X_test  <- sweep(X_test, 2, X_train_sds, "/")
    
    y_mean <- mean(y_train_raw)
    y_sd   <- sd(y_train_raw)
    if (y_sd == 0) y_sd <- 1
    y_train <- (y_train_raw - y_mean) / y_sd
    
    # TUNING & ESTIMATION
    if ((j - 1) %% up == 0) {
      start_time_star_post <- proc.time()
      l1_grid <- lambda_grid(X_train, y_train, nlambda, "node sparsity", A, sp_Astar, NULL)
      l2_grid <- lambda_grid(X_train, y_train, nlambda, "leaf sparsity", A, sp_Astar, NULL)
      
      startime_estimate <- startime_ic(X_train, y_train, P_list, Ki_list, post=TRUE,
                                       lambda1_grid = l1_grid, lambda2_grid = l2_grid,
                                       rho = rho, K_in = K_in, stop_thresh = stop_thresh,
                                       thresh = thresh, b2_final=FALSE)
      
      valid <- valid_lambdas(startime_estimate$ic_values, l1_grid, l2_grid)
      l1_grid <- lambda_grid(X_train, y_train, nlambda, "node sparsity", A, sp_Astar, lambda_min = valid$valid_lambda1)
      l2_grid <- lambda_grid(X_train, y_train, nlambda, "leaf sparsity", A, sp_Astar, lambda_min = valid$valid_lambda2)
      
      startime_estimate <- startime_ic(X_train, y_train, P_list, Ki_list, post=TRUE,
                                       lambda1_grid = l1_grid, lambda2_grid = l2_grid,
                                       rho = rho, K_in = K_in, stop_thresh = stop_thresh,
                                       thresh = thresh, b2_final=FALSE)
      
      last_lambda1 <- startime_estimate$lambda1_best
      last_lambda2 <- startime_estimate$lambda2_best
      time_star_post[j] <- (proc.time() - start_time_star_post)[["elapsed"]]
      
      start_time_simple <- proc.time()
      simple_est <- startime_ic(X_train, y_train, P_list, Ki_list, post=FALSE,
                                lambda1_grid = l1_grid, lambda2_grid = l2_grid,
                                rho = rho, K_in = K_in, stop_thresh = stop_thresh,
                                thresh = thresh, b2_final=TRUE)
      simple_last_lambda1 <- simple_est$lambda1_best
      simple_last_lambda2 <- simple_est$lambda2_best
      time_simple[j] <- (proc.time() - start_time_simple)[["elapsed"]]
      
      # MIDAS
      start_time_midasml <- proc.time()
      sg_lasso <- sglfit(x = X_train, y = y_train, gamma = 1)
      sg_lasso_mat <- as.matrix(sg_lasso$beta)
      ic_idx <- which.min(apply(sg_lasso_mat, 2, function(b) {
        ic_deriv(X_train, y_train, b, b, b, ic_type, thresh)
      }))
      midasml_last_lambda <- sg_lasso$lambda[as.numeric(ic_idx)]
      time_midasml[j] <- (proc.time() - start_time_midasml)[["elapsed"]]
    }
    
    # Estimation (Fixed Lambda)
    start_time_star_post <- proc.time()
    est_final <- startime(X_train, y_train, P_list, Ki_list, post=TRUE,
                          lambda1 = last_lambda1, lambda2 = last_lambda2,
                          rho = rho, K_in = K_in, stop_thresh = stop_thresh)
    time_star_post[j] <- time_star_post[j] + (proc.time() - start_time_star_post)[["elapsed"]]
    
    start_time_simple <- proc.time()
    est_simp_final <- startime(X_train, y_train, P_list, Ki_list, post=FALSE,
                               lambda1 = simple_last_lambda1, lambda2 = simple_last_lambda2,
                               rho = rho, K_in = K_in, stop_thresh = stop_thresh, b2_final=TRUE)
    time_simple[j] <- time_simple[j] + (proc.time() - start_time_simple)[["elapsed"]]
    
    start_time_midasml <- proc.time()
    sg_lasso <- sglfit(x = X_train, y = y_train, gamma = 1, lambda = midasml_last_lambda)
    time_midasml[j] <- time_midasml[j] + (proc.time() - start_time_midasml)[["elapsed"]]
    sg_lasso_matrix <- as.matrix(sg_lasso$beta)
    
    # Prediction
    pred_scaled <- X_test %*% est_final$beta
    predictions[j] <- pred_scaled * y_sd + y_mean
    betas[[j]] <- est_final$beta
    rownames(betas[[j]]) <- colnames(X)
    lambda1_vals[j] <- last_lambda1
    lambda2_vals[j] <- last_lambda2
    
    simp_pred_scaled <- X_test %*% est_simp_final$beta
    simple_predictions[j] <- simp_pred_scaled * y_sd + y_mean
    simple_betas[[j]] <- est_simp_final$beta
    rownames(simple_betas[[j]]) <- colnames(X)
    simple_lambda1_vals[j] <- simple_last_lambda1
    simple_lambda2_vals[j] <- simple_last_lambda2
    
    mid_pred_scaled <- X_test %*% sg_lasso_matrix
    midasml_predictions[j] <- mid_pred_scaled * y_sd + y_mean
    midasml_betas[[j]] <- sg_lasso_matrix
    midasml_lambda_vals[j] <- midasml_last_lambda
    
  }

  return(list(predictions = predictions, betas = betas, lambda1 = lambda1_vals, lambda2 = lambda2_vals,
              simple_predictions = simple_predictions, simple_betas = simple_betas, 
              simple_lambda1 = simple_lambda1_vals, simple_lambda2 = simple_lambda2_vals,
              midasml_predictions = midasml_predictions, midasml_betas = midasml_betas, 
              midasml_lambda_vals = midasml_lambda_vals,
              time_star_post = time_star_post, time_simple = time_simple, time_midasml = time_midasml))
}


ar1_rolling_predict <- function(X, y, window_size, h) {
  X <- as.matrix(X)
  n <- nrow(X)
  
  n_steps<- n - window_size - h
  
  predictions <- vector("numeric", n_steps)
  betas <- vector("list", n_steps)
  
  for (j in 1:n_steps) {
    i <- window_size + j
    
    train_start_X <- j
    train_end_X <- j + window_size - h - 1
    
    X_train_raw <- X[train_start_X:train_end_X, , drop=FALSE]
    
    y_train_raw <- y[(train_start_X + h):(train_end_X + h)]
    
    X_test_raw <- X[i, , drop=FALSE]
    
    # Standardize
    X_mean <- mean(X_train_raw); X_sd <- sd(X_train_raw); if(X_sd==0) X_sd<-1
    y_mean <- mean(y_train_raw); y_sd <- sd(y_train_raw); if(y_sd==0) y_sd<-1
    
    X_train <- (X_train_raw - X_mean)/X_sd
    X_test  <- (X_test_raw - X_mean)/X_sd
    y_train <- (y_train_raw - y_mean)/y_sd
    
    beta_hat <- ols(X = X_train, y = y_train)
    pred_scaled <- X_test %*% beta_hat
    predictions[j] <- pred_scaled * y_sd + y_mean
    betas[[j]] <- beta_hat
  }
  return(list(predictions = predictions, betas = betas))
}


random_walk_rolling_predict <- function(X,y, window_size, h) {
  n <- nrow(X)
  n_steps<- n - window_size - h
  predictions <- vector("numeric", n_steps)

  for (j in 1:n_steps) {
    i <- window_size + j
    predictions[j] <- y[(i-1)]
  }

  return(list(predictions=predictions))
}


evaluate_rolling_forecasts <- function(rolling_results, y, window_size, h) {
  
  n <- length(y)
  
  n_steps<- n - window_size - h
  
  first_target_idx <- window_size + h + 1
  
  test_indices <- seq(from = first_target_idx, by = 1, length.out = n_steps)

  if (tail(test_indices, 1) > length(y)) {
    warning("Evaluation indices exceed data length. Truncating.")
    valid_len <- sum(test_indices <= length(y))
    test_indices <- test_indices[1:valid_len]
  }
  
  y_true <- y[test_indices]
  
  if ("simple_predictions" %in% names(rolling_results)) {
    pred_star_post    <- unlist(rolling_results$predictions)
    pred_star_simple  <- unlist(rolling_results$simple_predictions)
    pred_midasml      <- unlist(rolling_results$midasml_predictions)
    
    mse_star_post    <- compute_mse(y_true, pred_star_post)
    mse_star_simple  <- compute_mse(y_true, pred_star_simple)
    mse_midasml      <- compute_mse(y_true, pred_midasml)
    
    return(list(mse_star_post=mse_star_post, mse_star_simple=mse_star_simple, mse_midasml=mse_midasml,
                y_true=y_true, pred_star_post=pred_star_post))
  } else {
    pred <- unlist(rolling_results$predictions)
    mse <- compute_mse(y_true, pred)
    return(list(mse=mse, y_true=y_true, pred=pred))
  }
}


macro_run <- function(X, y, window_size, h, up,
                      P_list, Ki_list, rho,
                      K_in, stop_thresh, ic_type, thresh,
                      nlambda, b2_final, dataset) {
  
  # Run Penalized (Standardized internally)
  penalized_results <- penalized_rolling_predict(X, y, window_size, h, up,
                                                 P_list, Ki_list, rho,
                                                 K_in, stop_thresh, ic_type, thresh,
                                                 nlambda, b2_final)
  
  penalized_eval <- evaluate_rolling_forecasts(penalized_results, y, window_size, h)
  
  rw_results <- NA
  ar1_results <- NA
  
  rw_eval        <- NA
  rw_eval$mse <- NA
  ar1_eval       <- NA
  ar1_eval$mse       <- NA
  
  if (dataset == "forecast"){
    rw_results <- random_walk_rolling_predict(X,y, window_size, h)
    ar1_results <- ar1_rolling_predict(X[, "GDP_qlag0", drop=FALSE], y, window_size, h)
    
    rw_eval        <- evaluate_rolling_forecasts(rw_results, y, window_size, h)
    ar1_eval       <- evaluate_rolling_forecasts(ar1_results, y, window_size, h)
  }
  
  return(list(
    penalized   = list(results = penalized_results, evaluation = penalized_eval),
    random_walk = list(results = rw_results,        evaluation = rw_eval),
    ar1         = list(results = ar1_results,       evaluation = ar1_eval)
  ))
}


get_Xy <- function(dataset, reduction) {
  # Map "forecast" user API to the underlying "raw" data files
  if (dataset == "forecast" || dataset == "raw") {
    if (reduction == "full") {
      load("data/macro/X_macro_raw.RData")
      load("data/macro/y_macro_raw.RData")
      return(list(X = X_raw, y = y_raw))
    } else if (reduction == "reduced") {
      load("data/macro/X_macro_raw_reduced.RData")
      load("data/macro/y_macro_raw.RData") 
      return(list(X = X_raw_reduced, y = y_raw))
    }
  } else if (dataset == "nowcast") {
    if (reduction == "full") {
      load("data/macro/X_macro_nowcast.RData")
      load("data/macro/y_macro_nowcast.RData")
      return(list(X = X_nowcast, y = y_nowcast))
    } else if (reduction == "reduced") {
      load("data/macro/X_macro_nowcast_reduced.RData")
      load("data/macro/y_macro_nowcast.RData") 
      return(list(X = X_nowcast_reduced, y = y_nowcast))
    }
  }
  stop(paste("Unknown dataset/reduction:", dataset, reduction))
}

build_groups <- function(dataset = c("nowcast", "forecast"),
                         reduction = c("full", "reduced")) {
  
  # Base Tree Definitions
  daily_tree    <- c(60, 12, 3, 1)
  weekly_tree   <- c(12, 3, 1)
  monthly_tree  <- c(3, 1)
  quarterly_tree<- c(1)
  
  daily_sizes    <- c(5, 4, 3, 1)
  weekly_sizes   <- c(4, 3, 1)
  monthly_sizes  <- c(3, 1)
  quarterly_sizes<- c(1)
  
  if (dataset == "forecast" || dataset == "raw"){
    if (reduction == "full") {
      no_daily_trees    <- 5   
      no_weekly_trees   <- 3   
      no_monthly_trees  <- 20  
      no_quarterly_trees<- 2   
    } else { 
      no_daily_trees    <- 0
      no_weekly_trees   <- 2   
      no_monthly_trees  <- 8   
      no_quarterly_trees<- 1   
    }
  } else if (dataset == "nowcast"){
    if (reduction == "full") {
      no_daily_trees    <- 5   
      no_weekly_trees   <- 3   
      no_monthly_trees  <- 20  
      no_quarterly_trees<- 2   
    } else { 
      no_daily_trees    <- 0
      no_weekly_trees   <- 2   
      no_monthly_trees  <- 8   
      no_quarterly_trees<- 1   
    }
  }
  
  P_list <- c(
    rep(list(daily_tree),     no_daily_trees),
    rep(list(weekly_tree),    no_weekly_trees),
    rep(list(monthly_tree),   no_monthly_trees),
    rep(list(quarterly_tree), no_quarterly_trees)
  )
  
  Ki_list <- c(
    rep(list(daily_sizes),     no_daily_trees),
    rep(list(weekly_sizes),    no_weekly_trees),
    rep(list(monthly_sizes),   no_monthly_trees),
    rep(list(quarterly_sizes), no_quarterly_trees)
  )
  
  list(P_list = P_list, Ki_list = Ki_list)
}