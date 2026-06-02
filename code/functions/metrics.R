# Defines and documents performance metrics for estimation methods,
#          including MSE, Adjusted Rand Index (ARI), and F1 score, as well as
#          utility functions for metric computation and model evaluation.


#' @description
#' Provides functions to compute key performance metrics for regression and clustering
#' models, including mean squared error (MSE), Adjusted Rand Index (ARI), F1 score,
#' and associated utility functions for model evaluation in simulation studies.
#'
#' @details
#' This file is intended to be sourced at the start of your simulation scripts.
#' All functions are self-contained and can be used for both AR and MIDAS pipelines.


#' Adjusted Rand Index (ARI)
#' 
#' Computes the similarity measure between two clusterings, correcting for chance.
#' 
#' @param labels_true True cluster labels (numeric/character vector)
#' @param labels_pred Predicted cluster labels (numeric/character vector)
#' @return Numeric value between -1 (no agreement) and 1 (perfect agreement)
#' 
#' @details
#' Implements the Hubert-Arabie formulation of ARI:
#' \deqn{ARI = \frac{\text{Index} - \text{Expected Index}}{\text{Max Index} - \text{Expected Index}}}
#' Handles edge cases where denominator approaches zero by returning 0.
#' 
adjusted_rand_index <- function(labels_true, labels_pred) {
  # Input validation
  stopifnot(
    length(labels_true) == length(labels_pred),
    length(labels_true) > 0,
    is.atomic(labels_true),
    is.atomic(labels_pred)
  )
  
  # Convert inputs to vectors if needed
  labels_true <- as.vector(labels_true)
  labels_pred <- as.vector(labels_pred)
  
  # Create contingency table
  cont_table <- table(labels_true, labels_pred)
  
  # Sum over all combinations in contingency table
  sum_comb_c <- sum(choose(cont_table, 2))
  
  # Sum combinations for row and column marginals
  sum_comb_rows <- sum(choose(rowSums(cont_table), 2))
  sum_comb_cols <- sum(choose(colSums(cont_table), 2))
  
  # Calculate total number of pairs
  n <- length(labels_true)
  total_comb <- choose(n, 2)
  
  # Calculate expected and maximum index
  expected_index <- (sum_comb_rows * sum_comb_cols) / total_comb
  max_index <- (sum_comb_rows + sum_comb_cols) / 2
  
  # Handle edge case where expected == max
  if (abs(max_index - expected_index) < .Machine$double.eps) {
    return(0)
  }
  
  # Calculate ARI
  ARI <- (sum_comb_c - expected_index) / (max_index - expected_index)
  return(ARI)
}


#' F1 Score for Coefficient Recovery
#' 
#' Computes the F1 score between true and estimated non-zero coefficients.
#' 
#' @param labels_true True coefficient vector (numeric)
#' @param labels_pred Predicted coefficient vector (numeric)
#' @return F1 score between 0 (no recovery) and 1 (perfect recovery)
#' 
#' @details
#' Calculates the harmonic mean of precision and recall for non-zero coefficient identification:
#' \deqn{F1 = \frac{2 \cdot TP}{2 \cdot TP + FP + FN}}
#' Where:
#' - TP: True positives (correct non-zero coefficients)
#' - FP: False positives (incorrect non-zero coefficients)
#' - FN: False negatives (missed non-zero coefficients)
#'
F1_score <- function(labels_true,labels_pred){ 
  # Convert inputs to vectors if needed
  labels_true <- as.vector(labels_true)
  labels_pred <- as.vector(labels_pred)
  
  # Input validation
  stopifnot(
    length(labels_true) == length(labels_pred),
    is.numeric(labels_true),
    is.numeric(labels_pred)
  )
  
  # Convert to binary (0/1) indicators
  true_binary <- as.integer(labels_true != 0)
  pred_binary <- as.integer(labels_pred != 0)
  
  # Calculate confusion matrix components
  tp <- sum(true_binary & pred_binary)
  fp <- sum(!true_binary & pred_binary)
  fn <- sum(true_binary & !pred_binary)
  
  # Handle division by zero
  denominator <- (2*tp + fp + fn)
  if (denominator == 0) return(0)
  
  # Return F1 score
  return((2 * tp) / denominator)
}


#' Compute Tree Aggregation Metrics
#' 
#' Calculates performance metrics for sparse tree-aggregated regression models.
#' 
#' @param X Design matrix (numeric matrix)
#' @param y Response vector (numeric vector)
#' @param taggr Model object from la_admm() (from the StarTime package)
#' @param true_beta True coefficient vector (numeric)
#' @param A Tree-encoding matrix (numeric matrix)
#' @return List containing:
#' \itemize{
#'   \item betadiff: MSE between estimated and true coefficients
#'   \item ari: Adjusted Rand Index for grouping similarity
#'   \item f1: F1 score for coefficient recovery
#' }
#' @export
compute_metrics <- function(X, y, taggr, true_beta, A) {
  beta_final <- taggr$beta
  gamma1_final <- taggr$gamma1
  beta2_final <- taggr$beta2
  
  # Check if beta_final is all NA or a single NA
  if (all(is.na(beta_final))) {
    return(list(
      betadiff = NA_real_,
      ari = NA_real_,
      f1 = NA_real_
    ))
  } else {
    return(list(
      betadiff = compute_mse(true_beta, beta_final),
      ari = adjusted_rand_index(round(true_beta,4),round(as.vector(A %*% gamma1_final),4)),
      f1 = F1_score(true_beta, beta2_final)
    ))
  }
}


#' Compute "Post" Model Metrics
#' 
#' Calculates performance metrics for post estimates.
#' 
#' @param X Design matrix (numeric matrix)
#' @param y Response vector (numeric vector)
#' @param estimated_beta Estimated coefficient vector (numeric)
#' @param true_beta True coefficient vector (numeric)
#' @return List containing:
#' \itemize{
#'   \item betadiff: MSE between estimated and true coefficients
#'   \item ari: Adjusted Rand Index for grouping similarity
#'   \item f1: F1 score for coefficient recovery
#' }
#' 
#' @export
compute_metrics_post <- function(X, y, estimated_beta, true_beta) {
  # Check if estimated_beta is all NA or single NA
  if (all(is.na(estimated_beta))) {
    return(list(
      betadiff = NA_real_,
      ari = NA_real_,
      f1 = NA_real_
    ))
  } else {
    return(list(
      betadiff = compute_mse(estimated_beta, true_beta),
      ari = adjusted_rand_index(round(as.vector(estimated_beta),4), round(as.vector(true_beta),4)),
      f1 = F1_score(true_beta, estimated_beta)
    ))
  }
}

#' Run Tree-Aggregated Regression
#' 
#' Wrapper function for sparse tree-aggregation estimation with automatic 
#' metric calculation.
#' 
#' @param X Design matrix (numeric matrix)
#' @param y Response vector (numeric vector)
#' @param lambda1 L1 sparsity penalty parameter (numeric)
#' @param lambda2 Tree-aggregation penalty parameter (numeric)
#' @param rho ADMM increment parameter (integer)
#' @param K_in Number of inner LA-ADMM loops (integer)
#' @param stop_thresh
#' @param P_list List containing the tree-specific vectors with the number of 
#' nodes on each level of the trees (list of vectors)
#' @param Ki_list List containing the tree-specific vectors with the group size 
#' on each level of the trees (list of vectors)
#' @param leaf_penalty Penalty type on the betas (string)
#' @param node_penalty Penalty type on the gammas (string)
#' @param b_init Initial beta values (numeric vector)
#' @param g_init Initial gamma values (numeric vector)
#' @param true_beta True coefficient vector (numeric)
#' @param A Tree-encoding matrix (numeric matrix)
#' @return List containing:
#' \itemize{
#'   \item beta: Final coefficient estimates
#'   \item gamma1: First gamma copy in the ADMM algorithm (related to grouping)
#'   \item beta2: Second beta copy in the ADMM algoithm (related to sparsity)
#'   \item metrics: Performance metrics from compute_metrics()
#' }
#' 
#' @export
perform_staggr <- function(X, y, lambda1, lambda2, rho, K_in, stop_thresh, 
                           P_list, Ki_list, b_init, g_init,
                           true_beta, A) {
  taggr <- admm_solver(
    X, y, lambda1 = lambda1, lambda2 = lambda2, rho = rho, K_in = K_in,
    stop_thresh = stop_thresh, P_list = P_list, Ki_list = Ki_list,
    b_init = b_init, g_init = g_init
  )
  metrics <- compute_metrics(X, y, taggr, true_beta, A)
  
  # Return beta, gamma1, beta2 plus metrics
  c(
    list(
      beta   = taggr$beta,
      gamma1 = taggr$gamma1,
      beta2  = taggr$beta2
    ),
    metrics
  )
}
