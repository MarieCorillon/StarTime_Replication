

#' Generate AR(p) Time Series Data
#' 
#' Simulates an autoregressive process of order p with Gaussian noise.
#' 
#' @param n Number of observations to generate (excluding burn-in period). (integer)
#' @param coef Numeric vector of AR coefficients (length p). (numeric vector)
#' @param sd Standard deviation of Gaussian noise. (numeric, default = 1)
#' @param burnin Number of initial observations to discard. (integer, default = 100)
#' 
#' @return A list containing:
#' \itemize{
#'   \item y: Numeric vector of length n (response series)
#'   \item X: Numeric matrix of dimension n x p (lagged predictors)
#' }
#' 
#' @details
#' The AR(p) process follows:
#' \deqn{y_t = \sum_{j=1}^p \phi_j y_{t-j} + \varepsilon_t}
#' where \eqn{\varepsilon_t \sim N(0, \sigma^2)}. 
#' 
generate_data <- function(n, coef, sd = 1, burnin = 100) {
  # Input validation
  stopifnot(
    n > 0,
    length(coef) >= 1,
    sd > 0,
    burnin >= 0
  )
  
  p <- length(coef)
  ntot <- n + burnin
  
  # Initialize series with noise
  series <- numeric(ntot)
  series[1:p] <- rnorm(p, mean = 0, sd = sd)  # Initial values
  
  # Generate AR process
  for (t in (p + 1):ntot) {
    ar_component <- sum(coef * series[(t-1):(t-p)])
    series[t] <- ar_component + rnorm(1, mean = 0, sd = sd)
  }
  
  # Discard burn-in and create lag matrix
  y_burnt <- series[(burnin + 1):ntot]
  y_embedded <- stats::embed(y_burnt, dimension = p + 1)
  
  return(list(
    y = y_embedded[, 1],
    X = y_embedded[, -1, drop = FALSE]
  ))
}




#' Generate High-Frequency AR(1) Process
#'
#' Simulates a stationary AR(1) process at high frequency with configurable parameters.
#'
#' @param n_low Number of low-frequency periods (post-burnin) (positive integer)
#' @param m Frequency ratio (high-freq observations per low-freq period) (positive integer)
#' @param burnin Number of initial low-frequency periods to discard (non-negative integer)
#' @param ar_coef AR(1) coefficient (numeric, between -1 and 1 for stationarity)
#' @param sd_innov Standard deviation of innovations (positive numeric, default = 1)
#' 
#' @return Numeric vector of high-frequency AR(1) series
#' 
#' @details
#' The AR(1) process follows:
#' \deqn{x_t = \phi x_{t-1} + \varepsilon_t}
#' where \eqn{\varepsilon_t \sim N(0, \sigma^2)}. The series is initialized from the stationary distribution.
#'
#' @examples
#' # Generate 10 low-frequency periods (after burning 5) with m=3
#' hf_data <- generate_hf_data(n_low = 10, m = 3, burnin = 5, ar_coef = 0.5)
generate_hf_data <- function(n_low, m, burnin, ar_coef, sd_innov = 1) {
  # Input validation
  stopifnot(
    n_low > 0,
    m > 0,
    burnin >= 0,
    abs(ar_coef) < 1,  # Stationarity condition
    sd_innov > 0
  )
  
  # Total low-frequency periods (including burnin)
  total_low <- n_low + burnin + 1
  
  # High-frequency series length
  n_hf <- total_low * m
  
  # Initialize from stationary distribution
  x <- numeric(n_hf)
  x[1] <- rnorm(1, mean = 0, sd = sd_innov / sqrt(1 - ar_coef^2))
  
  # Generate innovations
  eps <- rnorm(n_hf - 1, mean = 0, sd = sd_innov)
  
  # Simulate AR(1)
  for (t in 2:n_hf) {
    x[t] <- ar_coef * x[t - 1] + eps[t - 1]
  }
  
  x
}



#' Simulate ARDL-MIDAS Data and Fit MIDAS Model
#'
#' @param n Number of low-frequency periods (integer)
#' @param m Frequency ratio (integer)
#' @param burnin Burn-in period (integer)
#' @param X_ar List of AR(1) coefficients for exogenous variables (length 3 or 10) (list of numeric)
#' @param sd Standard deviation of noise (numeric)
#' @param weight_params List of MIDAS weight parameters (list)
#' @param r AR coefficients for y (numeric vector)
#' @param guess Starting values for MIDAS optimization (numeric vector)

#' @return List with y, X, omega, formula, midas_model

#' @details
#' Simulates MIDAS data for 3 or 10 exogenous variables, generates high-frequency data,
#' computes MIDAS weights, simulates the response variable with ARDL structure, aligns
#' time series, and fits a MIDAS regression model.
#'
simulate_ardl_midas <- function(n, m, burnin, X_ar, sd, weight_params, r, guess, omega) {
  stopifnot(length(X_ar) %in% c(3, 10))
  
  exVar <- length(X_ar)
  n_total <- n + burnin + m
  y_noburn <- rnorm(n_total)
  
  # Generate high-frequency X variables
  x_hf <- lapply(seq_len(exVar), function(k) {
    ts(generate_hf_data(n + m - 1, m, burnin, X_ar[[k]], sd), frequency = m)
  })
  
  # Build lagged X matrices
  X_lagged <- lapply(x_hf, function(x) mls(x, 0:(m-1), m))
  
  # Simulate ARDL-MIDAS process
  for (t in (length(r) + 1):n_total) {
    y_lags <- sum(r * y_noburn[(t - 1):(t - length(r))])
    x_effects <- sum(sapply(seq_len(exVar), function(k) {
      X_lagged[[k]][t, ] %*% omega[[k]]  # Last 7 will be 0 when exVar=10
    }))
    y_noburn[t] <- x_effects + y_lags + rnorm(1)
  }
  
  # Finalize series
  y <- ts(y_noburn[(burnin + m + 1):n_total], start = (burnin + m + 1), frequency = 1)
  
  # Align X variables
  X_aligned <- lapply(x_hf, function(x) window(x, start = start(y)))
  X_mat <- do.call(cbind, X_aligned)
  colnames(X_mat) <- paste0("X", seq_len(exVar))
  
  # Needed for midas_r to know what the Xs are
  for (i in seq_len(ncol(X_mat))) {
    assign(paste0("X", i), X_mat[, i], envir = environment())
  }
  assign("y", y, envir = environment())
  
  # Build formula string
  formula_terms <- paste0("mls(X", seq_len(ncol(X_mat)), ", 0:(m-1), m, nealmon)")
  formula_str <- paste(
    "y ~ mls(y, 1:length(r), 1)",      # AR term
    "+",                               # Plus sign
    paste(formula_terms, collapse = " + "),  # All X terms
    "- 1"                              # No intercept
  )
  formula <- as.formula(formula_str)
  
  
  # Prepare start list
  start_list <- setNames(replicate(ncol(X_mat), guess, simplify = FALSE), paste0("X", seq_len(ncol(X_mat))))
  
  # Now call safe_midas_r (no data argument)
  midas_model <- safe_midas_r(formula, start = start_list)
  
  
  list(
    y = y,
    X = X_mat,
    omega = omega,
    formula = formula,
    midas_model = midas_model
  )
}



#' Construct MIDAS Design Matrix
#'
#' @param y Response vector (numeric)
#' @param X Exogenous variables (matrix or data.frame)
#' @param m Frequency ratio (integer)
#' @param r AR coefficients for y (numeric vector)
#' @return List with X (design matrix) and y (aligned response)
construct_midas_design <- function(y, X, m, r) {
  n_exo <- if (is.matrix(X)) ncol(X) else 1
  # Build lagged exogenous matrices
  X_comp <- lapply(seq_len(n_exo), function(j) mls(X[, j], 0:(m - 1), m))
  names(X_comp) <- if (is.matrix(X)) colnames(X) else "X"
  
  # Add AR lags
  y_comp <- mls(y, 1:length(r), 1)
  X_mat <- do.call(cbind, c(list(y_comp), X_comp))
  # Remove initial rows lost to lagging
  rows_to_remove <- seq_len(length(r))
  list(X = X_mat[-rows_to_remove, ], y = y[-rows_to_remove])
}