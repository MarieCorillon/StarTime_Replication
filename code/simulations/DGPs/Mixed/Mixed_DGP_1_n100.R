# Mixed Simulation Parameters (DGP1)
# --- Core Simulation Parameters ---
n <- 100            # Low-frequency sample size
m <- 12             # Frequency ratio (e.g., monthly to annual)
burnin <- 200       # Burn-in period for initialization
sd <- 1             # Standard deviation of noise

# --- AR and MIDAS Parameters ---
r <- c(0.3, 0.01, 0, 0)    # AR coefficients for y
X_ar <- rep(0.2, 3)          # AR(1) coefficients for each X

# --- MIDAS Weight Parameters ---
weight_params <- list(
  c(1,1,3),   # For X1
  c(1,2,3),   # For X2
  c(1,2,2)    # For X3
)

# Compute MIDAS weights: first 3 variables use weights, others get zero
omega <- lapply(seq_len(length(X_ar)), function(k) {
    if (k <= 3) nbeta(p = weight_params[[k]], d = m) 
    else rep(0, m)
})

# --- Penalization Configuration ---
config <- list(
  penalties = list(
    ss = list(node = "node sparsity", leaf = "leaf sparsity")
  ),
  grid_length = 10,
  rho = 1,
  K_in = 1000,
  stop_thresh = 1e-05,
  ic_type = "BIC",
  thresh = 0.3
)

# --- Group Structure for Penalized Methods ---
P_list = list(c(4,1),c(12,4,1), c(12,4,1),c(12,4,1))
Ki_list = list(c(4,1),c(3,4,1), c(3,4,1), c(3,4,1))

# --- Derived Parameters (do not modify) ---
A <- create_matrix_A(P_list, Ki_list)
sp_Astar <- create_matrix_A_star(P_list, Ki_list)
b_init <- rep(0, sum_first_elements(P_list))
g_init <- rep(0, sum_all_elements(P_list))

# --- MIDAS Model Fitting Parameters ---
guess <- c(1, -0.5)    # Starting values for MIDAS weights