# AR Simulation Parameters (DGP2)
n <- 100
P_list <- list(c(20, 4, 1))
Ki_list <- list(c(5, 4, 1))
burnin <- 200
sd <- 1
beta <- c(0.5, rep(-0.1, 4), rep(0, 15))

# Penalization Configuration
config <- list(
  penalties = list(
    ss = list(node = "node sparsity", leaf = "leaf sparsity"),
    lasso = list(node = "node sparsity", leaf = "leaf sparsity")
  ),
  grid_length = 10,
  rho = 1,
  K_in = 1000,
  stop_thresh = 1e-05,
  ic_type = "BIC",
  thresh = 1
)

# Derived Parameters (do not modify)
b_init <- rep(0, sum_first_elements(P_list))
g_init <- rep(0, sum_all_elements(P_list))
A <- create_matrix_A(P_list, Ki_list)
sp_Astar <- create_matrix_A_star(P_list, Ki_list)


