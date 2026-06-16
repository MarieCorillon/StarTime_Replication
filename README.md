# Replication Package: Sparse Tree-Based Aggregation for Time Series Regressions

[![arXiv](https://img.shields.io/badge/arXiv-2606.03665-b31b1b.svg)](https://arxiv.org/abs/2606.03665)

> This repository contains the complete replication code for *Sparse Tree-Based Aggregation for Time Series Regressions*.  
> It includes all Monte Carlo simulations (AR and Mixed-Frequency) and both empirical applications (Financial and Macro).  
> **Runtime:** Full replication takes approximately 13.1 hours on a workstation equipped with an AMD Ryzen Threadripper PRO 5965WX, using 16 cores, with 128 GB of RAM.

---

## Table of Contents
1. [Getting Started](#getting-started)
2. [Repository Structure](#repository-structure)
3. [Installation](#installation)
4. [The Master Dashboard (`user_run_file.R`)](#the-master-dashboard-user_run_filer)
5. [Workflow Recipes](#workflow-recipes)
6. [Output Guide](#output-guide)
7. [Reproducibility Notes](#reproducibility-notes)
8. [Data](#data)
9. [FAQ & Troubleshooting](#faq--troubleshooting)

---

## Getting Started

### Prerequisites
You need **R** installed on your system.  
*(This code was developed and tested on R `= 4.5.2`. To check your installed version, open R and run `R.version.string`.)*

### Obtaining the Code
If you are familiar with Git, clone the repository:

```bash
git clone https://github.com/MarieCorillon/StarTime_Replication.git
cd StarTime_Replication
```

If you do not use Git, click the green **<> Code** button at the top of this page, select **Download ZIP**, and extract the archive to your local machine.

---

## Repository Structure

An overview of the key directories:

```
├── setup/                            # Environment setup
│   ├── install_dependencies.R        # Installs exact package versions via pak
│   └── StarTime_1.0.tar.gz          # Local package (not on CRAN)
├── code/
│   ├── functions/                    # Shared utilities
│   │   ├── sim_results_functions.R  # Plotting and aggregation helpers
│   │   ├── metrics.R                # MSE, ARI, F1 computation
│   │   ├── macro_functions.R        # Macro application engine
│   │   └── financial_functions.R    # Financial application engine
│   ├── simulations/
│   │   ├── DGPs/                    # Data Generating Process configurations
│   │   │   ├── AR/                  # AR(1)–AR(3) DGP scripts
│   │   │   └── Mixed/               # Mixed(1)-Mixed(3) DGP scripts
│   │   ├── AR_run.R                 # AR simulation runner
│   │   ├── Mixed_run.R              # Mixed-frequency simulation runner
│   │   ├── run_simulations.R        # Simulation controller
│   │   └── run_simulation_plots.R   # Generates simulation bar plots
│   └── applications/
│       ├── run_financial_app.R      # Financial estimation
│       ├── run_financial_plots.R    # Financial figures
│       ├── run_financial_tables.R   # Financial tables (MCS, DM)
│       ├── run_macro_app.R          # Macro estimation (nowcast/forecast)
│       ├── run_macro_plots.R        # Macro figures
│       ├── run_macro_tables.R       # Macro tables
│       └── run_macro_coefficient_summaries.R  # StarTime coefficient analysis
├── data/
│   ├── financial/
│   │   └── Variances_10min.csv      # Realized volatility input data
│   └── macro/                       # Macro input data and date indices
└── results/
    ├── data/
    │   ├── paper/                   # ⭐ PRE-COMPUTED results shipped with repo
    │   │   ├── simulations/
    │   │   │   ├── AR/
    │   │   │   └── Mixed/
    │   │   └── applications/
    │   │       ├── financial/
    │   │       └── macro/
    │   └── user/                    # Your newly computed results
    │       ├── simulations/
    │       │   ├── AR/
    │       │   └── Mixed/
    │       └── applications/
    │           ├── financial/
    │           └── macro/
    └── figures/
        ├── paper/                   # ⭐ PRE-COMPUTED figures shipped with repo
        └── user/                    # Your newly generated figures
```

### Important Notes on the Dual-Folder System
- **`results/data/paper/`** and **`results/figures/paper/`** contain the exact pre-computed results and figures from the paper.  
- **`results/data/user/`** and **`results/figures/user/`** are where **your** outputs go when running custom code.  
- The toggle `USE_PAPER_RESULTS` controls which `data` directory is read from/written to.  
- **Figures are always written to `results/figures/user/`** (to prevent accidental overwriting of published figures), but you can read pre-computed figures from `results/figures/paper/`.

---

## Installation

### Step 1 — Run the master file

Open `user_run_file.R` in your editor and source it:

```r
source("user_run_file.R")
```

Or from the terminal, navigate to the project folder and run:

```bash
Rscript user_run_file.R
```

The working directory is set automatically — you do not need to call `setwd()` yourself.

The script will then:
1. Check that a C++ compiler is available (required to build `StarTime`)
2. Install all R packages at the exact versions used in the paper
3. Launch the full replication pipeline

---

### Step 2 — Compiler requirement (Rtools / Xcode / build tools)

`StarTime` is a custom C++ package that must be compiled on your machine the first time it is installed. The script detects whether the necessary tools are present and guides you through the process automatically.

> **Note:** `StarTime` is also available as a standalone package at [github.com/MarieCorillon/StarTime](https://github.com/MarieCorillon/StarTime) if you want to use it independently from this replication.

#### Windows

The required toolchain is called **Rtools**. The `INSTALL_RTOOLS` toggle at the top of `user_run_file.R` controls what happens if it is missing (set to `TRUE` by default).

| Your situation | What happens |
|---|---|
| Rtools already installed | Nothing extra — compilation proceeds silently |
| Rtools not installed · `INSTALL_RTOOLS <- TRUE` (default) | The script downloads the correct Rtools version and launches an installer window. Click through the prompts (~2 minutes), then **restart R once** and run `user_run_file.R` again |
| Rtools not installed · `INSTALL_RTOOLS <- FALSE` | The script stops immediately and prints a download link and manual instructions |

#### macOS

Two components are required:

1. **Xcode Command Line Tools** — provides the C/C++ compiler. If not found, the script stops and instructs you to run the following in your Terminal, then restart R:

   ```bash
   xcode-select --install
   ```

2. **GNU Fortran** — required to compile `StarTime`'s Fortran sources. It is *not* included with Xcode Command Line Tools and must be installed separately. Download the installer matching your R version and CPU architecture (Apple Silicon vs. Intel) from the official CRAN tools page:

   <https://mac.r-project.org/tools/>

   After installation, restart R before re-running `user_run_file.R`.

#### Linux

The `r-base-dev` package (or equivalent) is required. If not found, the script stops and prints the appropriate command:

```bash
# Ubuntu / Debian
sudo apt-get install r-base-dev

# Fedora / CentOS
sudo yum install R-devel
```

---

### How many times will I need to run the file?

| Situation | Runs needed |
|---|---|
| Compiler already installed | **Once** — installs packages and runs the pipeline |
| Windows without Rtools (`INSTALL_RTOOLS <- TRUE`) | **Twice** — first run installs Rtools and asks for a restart; second run completes everything |
| macOS / Linux without build tools | **Twice** — first run prints the terminal command; after running it and restarting R, the second run completes everything |

After the first successful run, all packages are installed and subsequent runs go straight to the pipeline with no setup overhead.

---

## The Master Dashboard (`user_run_file.R`)

All execution flows through a single file: **`user_run_file.R`**. You do not need to edit any other script.

### Core Concepts

| Toggle | Purpose |
|--------|---------|
| `REPLICATE_PAPER_RESULTS` | Global override. If `TRUE`, runs the exact full paper grid and ignores custom settings. |
| `USE_PAPER_RESULTS` | Storage/read toggle. If `TRUE`, points to `results/data/paper/`. If `FALSE`, points to `results/data/user/`. |
| `USE_PARALLEL` | Enable/disable parallel computing (`FALSE` by default). The number of cores is set via `NUM_CORES`, defaulting to `detectCores() - 1`. |

### Safety Mechanism
If `USE_PAPER_RESULTS = TRUE` while any computation toggle (`RUN_ESTIMATIONS`, `RUN_SIMULATIONS`, or `RUN_STATISTICAL_TESTS`) is also `TRUE`, the script **aborts immediately** to protect the shipped paper results from accidental overwrite. To run computations, you must set `USE_PAPER_RESULTS <- FALSE`.

### Three Execution Modes

#### Mode A: Full Paper Replication
Set `REPLICATE_PAPER_RESULTS <- TRUE`. The script automatically:
- Runs all simulation DGPs (AR 1–3, Mixed 1–3, both `n=100` and `n=200`).
- Runs all application configurations (Financial full lag/window/h grid; Macro full nowcast/forecast/full/reduced/window grid).
- Generates all plots and tables.
- Results are saved to `results/data/user/` and `results/figures/user/`.

#### Mode B: Custom Execution
Set `REPLICATE_PAPER_RESULTS <- FALSE` and configure the granular toggles in Section 3 of the dashboard:

**Simulations**
- `RUN_SIMULATIONS <- TRUE`
- `RUN_SIM_AR` / `RUN_SIM_MIXED`
- `SELECTED_DGPS_AR` / `SELECTED_DGPS_MIXED` (e.g., `c("AR_DGP_2_n100")`)
- `M_SIM` (default 500)

**Applications**
- `APP_RUN_FINANCIAL <- TRUE` + `FIN_MAX_LAG`, `FIN_WINDOW_SIZE`, `FIN_H`
  - These are fed into `expand.grid()`, so `c(20,40) × c(1000,250,125) × c(1,5,20)` runs all combinations automatically.
- `APP_RUN_MACRO <- TRUE` + `MACRO_DATASET`, `MACRO_REDUCTION`, `MACRO_WINDOW`
  - Also uses `expand.grid()`

#### Mode C: Generate Plots & Tables from Existing Results
To produce figures or LaTeX tables **without re-running any estimation**:

```r
REPLICATE_PAPER_RESULTS <- FALSE
RUN_ESTIMATIONS         <- FALSE
RUN_SIMULATIONS         <- FALSE
GENERATE_PLOTS          <- TRUE   # and/or GENERATE_TABLES <- TRUE
USE_PAPER_RESULTS       <- TRUE   # Read from shipped paper results
```

### Critical Warning: Statistical Tests Require the Full Grid
The **Diebold-Mariano (DM)** and **Model Confidence Set (MCS)** tables for each application require **all configurations** of that application to be available simultaneously (they are computed across the full cross-section of models).

- If you want these tables, either:
  1. Set `USE_PAPER_RESULTS <- TRUE` (to use the shipped full results), **or**
  2. Ensure you have run the **entire** application grid (all lags/windows/horizons for Financial; all dataset/reduction/window combinations for Macro) with `USE_PAPER_RESULTS <- FALSE`.

Running a subset of configurations and then requesting tables will produce errors or incomplete output.

---

## Workflow Recipes

Below are four copy-paste ready configurations. Place them at the top of `user_run_file.R` (Section 2) and source the file.

### Recipe 1: Full Paper Replication (Zero Configuration)
```r
REPLICATE_PAPER_RESULTS <- TRUE # It overrides the rest of the toggles.
```

### Recipe 2: Run Only Mixed-Frequency Simulations for DGP 2
```r
REPLICATE_PAPER_RESULTS <- FALSE
USE_PAPER_RESULTS       <- FALSE

RUN_SIMULATIONS         <- TRUE
RUN_SIM_AR              <- FALSE
RUN_SIM_MIXED           <- TRUE
SELECTED_DGPS_MIXED     <- c("Mixed_DGP_2_n100", "Mixed_DGP_2_n200")
M_SIM                   <- 500

RUN_ESTIMATIONS         <- FALSE
GENERATE_PLOTS          <- TRUE   # Plots generated automatically
GENERATE_TABLES         <- FALSE
```

### Recipe 3: Generate Only Financial Plots from Pre-Computed Paper Results
```r
REPLICATE_PAPER_RESULTS <- FALSE
USE_PAPER_RESULTS       <- TRUE   # Read from paper/ folder

RUN_ESTIMATIONS         <- FALSE
RUN_SIMULATIONS         <- FALSE
GENERATE_PLOTS          <- TRUE
GENERATE_TABLES         <- FALSE

APP_RUN_FINANCIAL       <- TRUE
APP_RUN_MACRO           <- FALSE
APP_RUN_INTRO_GRAPH     <- FALSE
```

### Recipe 4: Custom Macro Application (Nowcast, Reduced, Window = 66) + Plots + Tables
```r
REPLICATE_PAPER_RESULTS <- FALSE
USE_PAPER_RESULTS       <- FALSE

RUN_ESTIMATIONS         <- TRUE
RUN_SIMULATIONS         <- FALSE
GENERATE_PLOTS          <- TRUE
GENERATE_TABLES         <- FALSE

APP_RUN_MACRO           <- TRUE
MACRO_DATASET           <- c("nowcast")
MACRO_REDUCTION         <- c("reduced")
MACRO_WINDOW            <- c(66)

APP_RUN_FINANCIAL       <- FALSE
APP_RUN_INTRO_GRAPH     <- FALSE
```

---

## Output Guide

After execution, your outputs will be in the following locations:

| Output Type | Location |
|-------------|----------|
| **Simulation data** (`.RData`) | `results/data/user/simulations/{AR,Mixed}/` |
| **Simulation bar plots** (`.pdf`) | `results/figures/user/simulations/{AR,Mixed}/` |
| **Financial app data** | `results/data/user/applications/financial/` |
| **Financial line & coefficient plots** | `results/figures/user/applications/financial/{lines,coefficients}/` |
| **Macro app data** | `results/data/user/applications/macro/` |
| **Macro plots** | `results/figures/user/applications/macro/` |
| **Macro coefficient summaries** | `results/data/user/applications/macro/coefficient_summaries/` |
| **Tables** | Printed to the R console as raw LaTeX |

Pre-computed paper results are preserved in `results/data/paper/` and `results/figures/paper/` and are never overwritten by user runs.

---

## Reproducibility Notes

- **Parallel random number generation:** Each Monte Carlo replication `i` is seeded deterministically inside the worker via `set.seed(i + 10000)`. This ensures identical results regardless of the number of parallel cores (`USE_PARALLEL`) or the machine architecture.
- **Exact dependency versions:** The `setup/install_dependencies.R` script forces specific CRAN package versions (e.g., `glmnet 4.1-8`, `ggplot2 4.0.1`, etc.) using `pak::pkg_install()`, removing reliance on lockfile restoration.
- **Adaptive parallelization:** `NUM_CORES` defaults to `parallel::detectCores() - 1`, but you may cap it on shared systems.

---

## Data

All data required to reproduce the empirical results are included in this repository under the **`data/`** directory, so the full pipeline runs without any additional downloads. The source, access details, and original publication (where applicable) of each dataset are documented below.

### Simulation Data
The Monte Carlo experiments rely on no external data. All series are generated at runtime from the data-generating processes (DGPs) defined in **`code/simulations/DGPs/`** (AR and Mixed-Frequency), and are fully reproducible from the deterministic seeds described in the [Reproducibility Notes](#reproducibility-notes).

### Financial Data
The financial application uses the 10-minute realized variances of the 30 Dow Jones Industrial Average constituents, provided in **`data/financial/Variances_10min.csv`**. These data form part of the replication materials accompanying:

> Hecq, A., Margaritella, L., and Smeekes, S. (2023). Granger Causality Testing in High-Dimensional VARs: A Post-Double-Selection Procedure. *Journal of Financial Econometrics*, 21(3), 915–958. https://doi.org/10.1093/jjfinec/nbab023

They were retrieved from the authors' public data repository (<https://drive.google.com/drive/folders/1Yr8CSZGqV9g_bK2yK--V21SmtKQmVVKf>), accessed on October 17, 2025.

### Macroeconomic Data
The macroeconomic application draws on publicly available series from the following sources, all accessed on November 7, 2025:

- **FRED** (Federal Reserve Bank of St. Louis);
- the **Economic Policy Uncertainty** index (<https://www.policyuncertainty.com>);
- **Yahoo Finance** (S&P 500 and Dow Jones Industrial Average indices).

The transformed and aligned dataset is provided in **`data/macro/aligned_macro_data.RData`**, exactly as described in the paper. The corresponding model-ready inputs, the predictor matrices (**`X_macro_*.RData`**) and target vectors (**`y_macro_*.RData`**), are already aligned and can be passed directly to the StarTime functions.

---

## FAQ & Troubleshooting

**Q: Where did my results go?**  
A: Check the `USE_PAPER_RESULTS` toggle. If it was `TRUE`, results are read from `results/data/paper/`. If `FALSE`, they are written to `results/data/user/`.

**Q: Can I run a single DGP?**  
A: Yes. Set `REPLICATE_PAPER_RESULTS <- FALSE`, `RUN_SIMULATIONS <- TRUE`, and fill `SELECTED_DGPS_AR` or `SELECTED_DGPS_MIXED` with your chosen DGP name(s).

**Q: Why are the DM or MCS tables empty/failing?**  
A: These tests require the *full* cross-section of models. Either (1) set `USE_PAPER_RESULTS <- TRUE` to use the shipped full results, or (2) ensure you have estimated all configurations for that application (e.g., all financial lag/window/horizon combinations).

---
