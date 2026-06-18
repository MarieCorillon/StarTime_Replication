# ==============================================================================
# 1. INITIALIZATION
# ==============================================================================

# Automatically set the working directory to the project root.
# Works when sourced interactively (source("user_run_file.R")) or run from
# the command line (Rscript user_run_file.R). No external packages required.
local({
  script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) {
      flag <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      if (length(flag)) dirname(normalizePath(sub("--file=", "", flag))) else NULL
    }
  )
  if (!is.null(script_dir) && nzchar(script_dir)) setwd(script_dir)
})

# Set to TRUE to automatically download and launch the Rtools installer on
# Windows if it is not already present. Rtools is required to compile StarTime.
# On macOS/Linux, clear manual instructions will be shown instead.
INSTALL_RTOOLS <- TRUE

source("setup/install_dependencies.R")

# Load required libraries
packages <- c("glmnet", "doParallel", "doRNG", "foreach", "ggplot2", 
              "dplyr", "tidyr", "patchwork", "readr", "StarTime",
              "RColorBrewer", "kableExtra", "MCS", "midasml", "stringr",
              "knitr", "forecast", "midasr")
invisible(lapply(packages, library, character.only = TRUE, quietly = TRUE))

# ==============================================================================
# 2. USER CONTROLS (THE DASHBOARD)
# ==============================================================================

# --- 3A. Paper Replication ---
REPLICATE_PAPER_RESULTS <- TRUE

# Global Parallel Toggle for estimations
USE_PARALLEL <- FALSE
NUM_CORES <- max(1, parallel::detectCores() - 1)

# --- 3B. Global Output Settings (What to generate) ---
# Set these to TRUE to execute specific phases of the pipeline.
RUN_ESTIMATIONS       <- FALSE
RUN_STATISTICAL_TESTS <- FALSE
GENERATE_PLOTS        <- FALSE
GENERATE_TABLES       <- FALSE    # Prints the tables right after saving the tests

# --- 3C. Simulation Controls ---
RUN_SIMULATIONS   <- FALSE
M_SIM             <- 500       

RUN_SIM_AR        <- FALSE
SELECTED_DGPS_AR  <- c("AR_DGP_1_n100", "AR_DGP_2_n100", "AR_DGP_2_n100", 
                       "AR_DGP_1_n200", "AR_DGP_2_n200", "AR_DGP_3_n200") 

RUN_SIM_MIXED     <- FALSE
SELECTED_DGPS_MIXED  <- c("Mixed_DGP_1_n100","Mixed_DGP_2_n100","Mixed_DGP_3_n100",
                          "Mixed_DGP_1_n200","Mixed_DGP_2_n200","Mixed_DGP_3_n200") 

# --- 3D. Application Controls (Which modules to target) ---
APP_RUN_INTRO_GRAPH <- FALSE
APP_RUN_FINANCIAL   <- FALSE
# Flexibly accept single values or multiple values using c()
FIN_MAX_LAG         <- c(20,40)         
FIN_WINDOW_SIZE     <- c(1000,250,125)       
FIN_H               <- c(1, 5, 20)         

APP_RUN_MACRO       <- FALSE
MACRO_DATASET       <- c("nowcast","forecast")          
MACRO_REDUCTION     <- c("reduced","full")          
MACRO_WINDOW        <- c(105,66)                

# --- 3E. Paper vs User Results Directory ---
USE_PAPER_RESULTS <- TRUE             

if (USE_PAPER_RESULTS && !REPLICATE_PAPER_RESULTS) {
  RESULTS_DIR <- "results/data/paper"
  
  # Safety check: Prevent running heavy computations if we are pointing to the paper's folder
  if (RUN_ESTIMATIONS || RUN_SIMULATIONS || RUN_STATISTICAL_TESTS) {
    stop("SAFETY ABORT: You have computation toggles = TRUE while USE_PAPER_RESULTS = TRUE. 
          This would overwrite the original paper data/tests. 
          Set execution toggles to FALSE to just print tables/plots, or USE_PAPER_RESULTS <- FALSE to compute new data.")
  }
} else {
  RESULTS_DIR <- "results/data/user"
}

# Figures ALWAYS go to the user folder
FIGURES_DIR <- "results/figures/user" 

# ==============================================================================
# 3. EXECUTION LOGIC (Main Driver)
# ==============================================================================

if (REPLICATE_PAPER_RESULTS) {
  message("=========================================================")
  message("FULL PAPER REPLICATION MODE ACTIVATED")
  message("Overrides active: Running all grids, applications, and outputs.")
  message("=========================================================")
  
  # Force Output Flags
  GENERATE_PLOTS  <- TRUE
  GENERATE_TABLES <- TRUE
  
  # Run Simulations
  message("Running All Simulations...")
  source("code/simulations/run_simulations.R")
  
  # Run Applications
  message("Running Introductory Graph...")
  source("code/applications/run_intro_graph.R")
  
  message("Running Full Financial Application...")
  source("code/applications/run_financial_app.R")
  
  message("Running Full Macro Application...")
  source("code/applications/run_macro_app.R")
  source("code/applications/run_macro_coefficient_summaries.R")
  
  # Generate All Outputs
  message("Generating Paper Plots and Tables...")
  source("code/applications/run_financial_plots.R")
  source("code/applications/run_financial_tables.R")
  source("code/applications/run_macro_plots.R")
  source("code/applications/run_macro_tables.R")
  
} else {
  message("=========================================================")
  message("CUSTOM EXECUTION MODE")
  message("=========================================================")
  
  # ---------------------------------------------------------
  # Phase 1: Model Estimations & Data Generation
  # ---------------------------------------------------------
  if (RUN_ESTIMATIONS) {
    message("--- Running Estimations ---")
    if (APP_RUN_FINANCIAL) source("code/applications/run_financial_app.R")
    if (APP_RUN_MACRO) {
      source("code/applications/run_macro_app.R")
      source("code/applications/run_macro_coefficient_summaries.R")
    }
  }
  
  if (RUN_SIMULATIONS) {
    message("--- Running Simulations ---")
    source("code/simulations/run_simulations.R")
  }
  
  # ---------------------------------------------------------
  # Phase 2: Output Generation (Plots)
  # ---------------------------------------------------------
  if (GENERATE_PLOTS) {
    message("--- Generating Plots ---")
    if (APP_RUN_INTRO_GRAPH) source("code/applications/run_intro_graph.R")  
    if (APP_RUN_FINANCIAL)   source("code/applications/run_financial_plots.R")
    if (APP_RUN_MACRO)       source("code/applications/run_macro_plots.R")
    if (RUN_SIM_AR || RUN_SIM_MIXED) source("code/simulations/run_simulation_plots.R")
  }
  
  # ---------------------------------------------------------
  # Phase 3: Output Generation (Tables)
  # ---------------------------------------------------------
  if (GENERATE_TABLES) {
    message("--- Generating Tables ---")
    if (APP_RUN_FINANCIAL) {
      source("code/applications/run_financial_tables.R")
    }
    if (APP_RUN_MACRO) {
      source("code/applications/run_macro_tables.R")
    }
  }
}

message("=========================================================")
message("PIPELINE COMPLETED SUCCESSFULLY")
message("=========================================================")