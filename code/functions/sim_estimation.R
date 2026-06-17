#' Robust MIDAS Model Estimation with Randomized Restarts
#'
#' Attempts to fit a MIDAS regression model, retrying with randomized starting 
#' values until convergence is achieved.
#'
#' @param formula MIDAS regression formula (formula)
#' @param start List of starting values for optimization (list)
#'
#' @return Fitted MIDAS model object (from \code{midas_r})
#'
#' @details
#' The function repeatedly calls \code{midas_r} with the provided formula and starting values.
#' If an error occurs (e.g., non-convergence), it perturbs the starting values and retries.
#' This continues until a successful fit is obtained.
#'
#' @seealso \code{\link[midasr]{midas_r}}
#'
safe_midas_r <- function(formula, start, max_attempts = 100) {
  initial_start <- start
  for (attempt in seq_len(max_attempts)) {  # Capped so it can never hang
    result <- tryCatch({
      midasr::midas_r(formula, start = start)
    }, error = function(e) {
      message("Attempt ", attempt, " failed: ", e$message)
      return(NULL)
    })

    if (!is.null(result)) {
      message("Success on attempt ", attempt)
      return(result)  # If successful, return the result
    }

    # Perturb the starting values and retry. Each element must be perturbed
    # individually. 
    start <- lapply(initial_start, function(s) {
      runif(length(s), min = 0.5, max = 1.5) * s
    })
  }
  stop("safe_midas_r: MIDAS model did not converge after ", max_attempts,
       " attempts.")
}



#' Configure Lambda Grids for Penalized Estimation
#'
#' Helper for grid search setup for each penalty type.
#'
#' @param penalty_name Name of penalty type. (character)
#' @param X Design matrix. (numeric matrix)
#' @param y Response vector. (numeric vector)
#' @param penalties Penalty specification. (list)
#' @param grid_length Number of lambda grid points. (integer)
#' @param A Tree encoding matrix. (numeric matrix)
#' @param sp_Astar Tree structure matrix for StarTime. (numeric matrix)
#'
#' @return List containing:
#' \itemize{
#'   \item l1: Sparsity penalty grid (numeric vector)
#'   \item l2: Aggregation penalty grid (numeric vector)
#' }
#' 
#' @export
configure_lambda_grids <- function(penalty_name, X, y, penalties, grid_length,
                                   A, P_list, sp_Astar) {
  list(
    l1 = if(penalty_name != "lasso") {
      lambda_grid(
        X, y, grid_length, penalties$node, A, sp_Astar,
        lambda_min_ratio = NULL
      )
    } else 0,
    
    l2 = 
      lambda_grid(
        X, y, grid_length, penalties$leaf, A, sp_Astar,
        lambda_min_ratio = NULL
      )
  )
}



#' Process Penalized Estimation Methods for AR and MIDAS Simulation
#'
#' Runs grid search and post-estimation for each penalty, returning metrics and estimates.
#'
#' @param X Design matrix. (numeric matrix)
#' @param y Response vector. (numeric vector)
#' @param penalties List of penalty specifications. (list)
#' @param true_beta True coefficient vector. (numeric vector)
#' @param grid_length Number of lambda grid points. (integer)
#' @param rho ADMM penalty parameter. (numeric)
#' @param K_in Number of inner ADMM iterations. (integer)
#' @param P_list List of tree node counts per level. (list of integer vectors)
#' @param Ki_list List of tree group sizes per level. (list of integer vectors)
#' @param A Tree encoding matrix. (numeric matrix)
#' @param sp_Astar Tree structure matrix for StarTime. (numeric matrix)
#' @param b_init Initial beta values. (numeric vector)
#' @param g_init Initial gamma values. (numeric vector)
#' @param ic_type Information criterion type. (character)
#' @param thresh Threshold for IC calculation. (numeric)
#' @param stop_thresh
#'
#' @return A list containing:
#' \itemize{
#'   \item oracle_simple: Oracle-selected metrics from simple sparse 
#'   tree-aggregation (list)
#'   \item ic_simple: IC-selected metrics (list)
#'   \item oracle_simple_betas: Oracle-selected coefficients (vector)
#'   \item ic_simple_betas: IC-selected coefficients (vector)
#'   \item oracle_post: Oracle-selected metrics from post-estimated sparse
#'    tree-aggregation (list)
#'   \item ic_post: IC-selected metrics (list)
#'   \item oracle_post_betas: Oracle-selected coefficients (vector)
#'   \item ic_post_betas: IC-selected coefficients (vector)
#' }
#'

process_penalties <- function(X, y, penalties, true_beta, grid_length, rho, K_in,
                              stop_thresh,
                              P_list, Ki_list, A, sp_Astar, ic_type, 
                              b_init, g_init,
                              thresh) {
  
  oracle_simple_results <- list()
  oracle_post_results   <- list()
  ic_simple_results     <- list()
  ic_post_results       <- list()
  
  oracle_simple_betas   <- list()
  ic_simple_betas       <- list()
  oracle_post_betas     <- list()
  ic_post_betas         <- list()
  
  oracle_simple_gamma1  <- list()
  ic_simple_gamma1      <- list()
  oracle_post_gamma1    <- list()
  ic_post_gamma1        <- list()
  
  oracle_simple_beta2   <- list()
  ic_simple_beta2       <- list()
  oracle_post_beta2     <- list()
  ic_post_beta2         <- list()
  
  for (penalty_name in names(penalties)) {
    lambda_config <- configure_lambda_grids(
      penalty_name, X, y, penalties[[penalty_name]], 
      grid_length, A, sp_Astar
    )
    
    beta_grid <- admm_grid(
      X, y, lambda_config$l1, lambda_config$l2,
      rho, K_in, stop_thresh, P_list, Ki_list, b_init, g_init
    )
    
    ## (Simple) Oracle
    oracle_simple <- mse_rand(beta_grid, true_beta)
    oracle_staggr <- perform_staggr(
      X, y,
      oracle_simple$best_lambda1, oracle_simple$best_lambda2,
      rho, K_in, stop_thresh,
      P_list, Ki_list, b_init, g_init, true_beta, A
    )
    
    oracle_simple_results[[penalty_name]] <- oracle_staggr[c("betadiff", "ari", "f1")]
    oracle_simple_betas[[penalty_name]]   <- oracle_staggr$beta
    oracle_simple_gamma1[[penalty_name]]  <- oracle_staggr$gamma1
    oracle_simple_beta2[[penalty_name]]   <- oracle_staggr$beta2
    
    ## (Simple) IC
    ic_simple <- ic_select_apply(X, y, A, beta_grid, ic_type, thresh)
    ic_staggr <- perform_staggr(
      X, y,
      ic_simple$best_lambda1, ic_simple$best_lambda2,
      rho, K_in, stop_thresh,
      P_list, Ki_list, b_init, g_init, true_beta, A
    )
    
    ic_simple_results[[penalty_name]] <- ic_staggr[c("betadiff", "ari", "f1")]
    ic_simple_betas[[penalty_name]]   <- ic_staggr$beta
    ic_simple_gamma1[[penalty_name]]  <- ic_staggr$gamma1
    ic_simple_beta2[[penalty_name]]   <- ic_staggr$beta2
    
    ## Post-processing
    post_beta <- post_estimation(X, y, beta_grid, A, P_list)
    
    # (Post) Oracle
    oracle_post <- mse_rand(post_beta, true_beta)
    oracle_post_beta <- post_beta$beta_matrix[, oracle_post$best_row]
    oracle_post_metrics <- compute_metrics_post(X, y, oracle_post_beta, true_beta)
    
    oracle_post_results[[penalty_name]] <- oracle_post_metrics
    oracle_post_betas[[penalty_name]]   <- oracle_post_beta
    # No gamma1/beta2 here – post_estimation already gives final betas only
    
    # (Post) IC
    ic_post <- ic_select_post(X, y, post_beta, ic_type, thresh)
    ic_post_beta <- post_beta$beta_matrix[, ic_post$best_row]
    ic_post_metrics <- compute_metrics_post(X, y, ic_post_beta, true_beta)
    
    ic_post_results[[penalty_name]] <- ic_post_metrics
    ic_post_betas[[penalty_name]]   <- ic_post_beta
  }
  
  list(
    oracle_simple        = oracle_simple_results,
    ic_simple            = ic_simple_results,
    oracle_post          = oracle_post_results,
    ic_post              = ic_post_results,
    
    oracle_simple_betas  = oracle_simple_betas,
    ic_simple_betas      = ic_simple_betas,
    oracle_post_betas    = oracle_post_betas,
    ic_post_betas        = ic_post_betas,
    
    oracle_simple_gamma1 = oracle_simple_gamma1,
    ic_simple_gamma1     = ic_simple_gamma1,
    oracle_simple_beta2  = oracle_simple_beta2,
    ic_simple_beta2      = ic_simple_beta2
  )
}



#' Process Ridge Regression for AR Simulation
#'
#' Fits a ridge regression model (using \code{glmnet::glmnet} with \code{alpha = 0})
#' and selects the best solution by both oracle (lowest MSE to true coefficients)
#' and information criterion (IC). Computes performance metrics for both.
#'
#' @param X Design matrix. (numeric matrix)
#' @param y Response vector. (numeric vector)
#' @param true_beta True coefficient vector. (numeric vector)
#' @param ic_type Information criterion type. (character)
#' @param thresh Threshold for IC calculation. (numeric)
#'
#' @return List containing:
#' \itemize{
#'   \item oracle: Oracle-selected metrics (list)
#'   \item oracle_beta: Oracle-selected coefficients (numeric vector)
#'   \item ic: IC-selected metrics (list)
#'   \item ic_beta: IC-selected coefficients (numeric vector)
#' }
#'
#' @export
process_ridge <- function(X, y, true_beta, ic_type, thresh) {
  # Fit ridge regression path
  ridge_fit <- glmnet::glmnet(X, y, alpha = 0)
  ridge_betas <- as.matrix(ridge_fit$beta)
  
  # Oracle: select beta with lowest MSE to true_beta
  oracle_idx <- which.min(apply(ridge_betas, 2, compute_mse, true_beta))
  # IC: select beta with lowest information criterion
  ic_idx <- which.min(apply(ridge_betas, 2, function(b) {
    ic_deriv(X, y, b, b, b, ic_type, thresh)
  }))
  
  list(
    oracle = compute_metrics_post(X, y, ridge_betas[, oracle_idx], true_beta),
    oracle_beta = ridge_betas[, oracle_idx],
    ic = compute_metrics_post(X, y, ridge_betas[, ic_idx], true_beta),
    ic_beta = ridge_betas[, ic_idx]
  )
}



#' Process MIDAS-ML Estimation
#'
#' @param X Design matrix (numeric matrix)
#' @param y Response vector (numeric vector)
#' @param true_beta True coefficient vector (numeric vector)
#' @param ic_type Information criterion type (character)
#' @param thresh Threshold for IC calculation (numeric)
#' @return List of MIDAS-ML results
#' @seealso \code{\link[sglr]{sglfit}}
process_midas_ml <- function(X, y, true_beta, ic_type, thresh) {
  # Fit MIDAS-ML model
  sg_lasso <- sglfit(x = X, y = y, gamma = 1)
  sg_lasso_matrix <- as.matrix(sg_lasso$beta)
  
  # Oracle selection
  oracle_idx <- which.min(apply(sg_lasso_matrix, 2, compute_mse, true_beta))
  oracle_beta <- sg_lasso_matrix[, oracle_idx]
  
  # IC selection
  ic_idx <- which.min(apply(sg_lasso_matrix, 2, function(b) {
    ic_deriv(X, y, b, b, b, ic_type, thresh)
  }))
  ic_beta <- sg_lasso_matrix[, ic_idx]
  
  list(
    oracle = compute_metrics_post(X, y, oracle_beta, true_beta),
    oracle_beta = oracle_beta,
    ic = compute_metrics_post(X, y, ic_beta, true_beta),
    ic_beta = ic_beta
  )
}



#' Run AR Simulation with Penalized and Benchmark Methods
#'
#' Runs a full simulation for an AR(p) DGP, including penalized estimation,
#' OLS, and ridge regression, and computes performance metrics.
#'
#' @param sim_id Identifier for the simulation run. (integer)
#'
#' @return A list with:
#' \itemize{
#'   \item data: List with simulated y (numeric vector) and X (numeric matrix)
#'   \item metrics: Data frame of performance metrics for all methods
#'   \item beta_estimates: List of coefficient estimates for all methods
#' }
#'
run_AR_simulation <- function(sim_id) {
  # Validate required parameters exist in environment
  required_objects <- c("n", "P_list", "Ki_list", "burnin", "sd", "beta",
                        "b_init", "g_init", "A", "sp_Astar")
  
  
  stopifnot(
    is.numeric(sim_id),
    length(sim_id) == 1,
    all(sapply(required_objects, exists)),
    exists("config") && is.list(config),
    all(c("grid_length", "rho", "K_in", "stop_thresh", "ic_type", "thresh") %in% names(config))
  )
  
  # 1. Data Generation 
  data <- generate_data(
    n = n, 
    coef = beta, 
    sd = sd, 
    burnin = burnin
  )
  X <- data$X
  y <- data$y
  
  # 2. Penalized Estimation 
  penalty_results <- process_penalties(
    X = X,
    y = y,
    penalties = config$penalties,
    true_beta = beta,
    grid_length = config$grid_length,
    rho = config$rho,
    K_in = config$K_in,
    stop_thresh = config$stop_thresh,
    P_list = P_list,
    Ki_list = Ki_list,
    A = A,
    sp_Astar = sp_Astar,
    b_init = b_init,
    g_init = g_init,
    ic_type = config$ic_type,
    thresh = config$thresh
  )
  
  # 3. Basic OLS 
  ols_est <- ols(X, y)
  ols_results <- list(
    basic = compute_metrics_post(X, y, ols_est, beta)
  )
  
  # 4. Ridge Regression 
  ridge_results <- process_ridge(
    X = X, 
    y = y, 
    true_beta = beta, 
    ic_type = config$ic_type, 
    thresh = config$thresh)
  
  # 5. Compile Results 
  list(
    data = data,
    metrics = data.frame(
      sim_id = sim_id,
      sapply(penalty_results$oracle_simple, unlist),
      sapply(penalty_results$ic_simple, unlist),
      sapply(penalty_results$oracle_post, unlist),
      sapply(penalty_results$ic_post, unlist),
      ols_basic    = unlist(ols_results$basic),
      ridge_oracle = unlist(ridge_results$oracle),
      ridge_ic     = unlist(ridge_results$ic)
    ),
    beta_estimates = list(
      ss_oracle_simple      = penalty_results$oracle_simple_betas$ss,
      ss_ic_simple          = penalty_results$ic_simple_betas$ss,
      ss_oracle_post        = penalty_results$oracle_post_betas$ss,
      ss_ic_post            = penalty_results$ic_post_betas$ss,
      lasso_oracle_simple   = penalty_results$oracle_simple_betas$lasso,
      lasso_ic_simple       = penalty_results$ic_simple_betas$lasso,
      lasso_oracle_post     = penalty_results$oracle_post_betas$lasso,
      lasso_ic_post         = penalty_results$ic_post_betas$lasso,
      ols_basic             = ols_est,
      ridge_oracle          = ridge_results$oracle_beta,
      ridge_ic              = ridge_results$ic_beta,
      
      ss_oracle_simple_gamma1 = penalty_results$oracle_simple_gamma1$ss,
      ss_ic_simple_gamma1     = penalty_results$ic_simple_gamma1$ss,
      ss_oracle_simple_beta2  = penalty_results$oracle_simple_beta2$ss,
      ss_ic_simple_beta2      = penalty_results$ic_simple_beta2$ss
    )
  )
}




#' Run Mixed Simulation with Penalized, Benchmark, and MIDAS Methods
#'
#' @param sim_id Simulation identifier (integer)
#' @return List containing simulation results
run_Mixed_simulation <- function(sim_id) {
  # Validate required parameters exist in environment
  required_objects <- c("n", "m", "burnin", "X_ar", "sd", "weight_params", "r",
                        "guess", "config", "P_list", "Ki_list", "A", "sp_Astar", 
                        "b_init", "g_init", "omega")
  stopifnot(
    is.numeric(sim_id),
    length(sim_id) == 1,
    all(sapply(required_objects, exists)),
    is.list(config),
    all(c("penalties", "grid_length", "rho", "K_in", "stop_thresh", "ic_type", "thresh") %in% names(config))
  )
  
  # 1. Data Generation 
  data <- simulate_ardl_midas(n, m, burnin, X_ar, sd, weight_params, r, guess, omega)
  working_data <- construct_midas_design(data$y, data$X, m, r)
  X <- working_data$X
  y <- working_data$y
  beta <- c(r, unlist(omega))
  
  # 2. Penalized Estimation 
  penalty_results <- process_penalties(
    X = X,
    y = y,
    penalties = config$penalties,
    true_beta = beta,
    grid_length = config$grid_length,
    rho = config$rho,
    K_in = config$K_in,
    stop_thresh = config$stop_thresh,
    P_list = P_list,
    Ki_list = Ki_list,
    A = A,
    sp_Astar = sp_Astar,
    b_init = b_init,
    g_init = g_init,
    ic_type = config$ic_type,
    thresh = config$thresh
  )
  
  # 3. OLS Benchmark 
  if (dim(X)[1] > dim(X)[2]){
    ols_est <- ols(X, y)
  } else {
    ols_est <- NA
  }
  ols_results <- compute_metrics_post(X, y, ols_est, beta)
  
  # 4. MIDAS-ML Estimation 
  midas_ml_results <- process_midas_ml(X, y, beta, config$ic_type, config$thresh)
  
  # 5. MIDAS Model 
  midas_fit <- data$midas_model
  midas_beta <- midas_fit$midas_coefficients
  midas_results <- compute_metrics_post(X, y, midas_beta, beta)
  
  # 6. Compile Results 
  list(
    data = working_data,
    metrics = data.frame(
      sim_id = sim_id,
      sapply(penalty_results$oracle_simple, unlist),
      sapply(penalty_results$ic_simple, unlist),
      sapply(penalty_results$oracle_post, unlist),
      sapply(penalty_results$ic_post, unlist),
      ols_basic        = unlist(ols_results),
      midas_ml_oracle  = unlist(midas_ml_results$oracle),
      midas_ml_ic      = unlist(midas_ml_results$ic),
      midas            = unlist(midas_results)
    ),
    beta_estimates = list(
      penalty_oracle_simple = penalty_results$oracle_simple_betas,
      penalty_ic_simple     = penalty_results$ic_simple_betas,
      penalty_oracle_post   = penalty_results$oracle_post_betas,
      penalty_ic_post       = penalty_results$ic_post_betas,
      ols_basic             = ols_est,
      midas_ml_oracle       = midas_ml_results$oracle_beta,
      midas_ml_ic           = midas_ml_results$ic_beta,
      midas_beta            = midas_beta,
      
      penalty_oracle_simple_gamma1 = penalty_results$oracle_simple_gamma1,
      penalty_ic_simple_gamma1     = penalty_results$ic_simple_gamma1,
      penalty_oracle_simple_beta2  = penalty_results$oracle_simple_beta2,
      penalty_ic_simple_beta2      = penalty_results$ic_simple_beta2
    )
  )
}
