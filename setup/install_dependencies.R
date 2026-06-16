# ==============================================================================
# Dependency installer
# Installs all required packages at the exact versions used in the paper.
# Uses pak for fast, reliable installation from CRAN archives.
# ==============================================================================

# --- Step 1: Check for compilation tools (required to build StarTime) --------
#
# StarTime is distributed as a C++ source package (.tar.gz) and must be
# compiled on the user's machine. This requires:
#   Windows : Rtools  (https://cran.r-project.org/bin/windows/Rtools/)
#   macOS   : Xcode Command Line Tools  (run: xcode-select --install)
#   Linux   : r-base-dev  (run: sudo apt-get install r-base-dev)
#
# Set INSTALL_RTOOLS <- TRUE in user_run_file.R to attempt an automatic
# Rtools installation on Windows (launches an installer, ~2 min, one restart).

if (!requireNamespace("pkgbuild", quietly = TRUE)) install.packages("pkgbuild")

if (!pkgbuild::has_build_tools(debug = FALSE)) {

  if (.Platform$OS.type == "windows") {

    if (exists("INSTALL_RTOOLS") && isTRUE(INSTALL_RTOOLS)) {
      message("Rtools not found. Downloading and launching the Rtools installer...")
      message("Please follow the installer prompts (click Next/Install).")
      message("Once complete, RESTART R and run user_run_file.R again.")
      if (!requireNamespace("installr", quietly = TRUE)) install.packages("installr")
      installr::install.Rtools()
      # Attempt to register the newly installed Rtools without a restart
      pkgbuild::find_rtools(TRUE)
      if (!pkgbuild::has_build_tools(debug = FALSE)) {
        stop(
          "Rtools installation may not be complete yet.\n",
          "Please restart R and run user_run_file.R again."
        )
      }
    } else {
      stop(
        "Rtools is required to compile StarTime but was not found.\n\n",
        "Option 1 (automatic): set INSTALL_RTOOLS <- TRUE in user_run_file.R\n",
        "Option 2 (manual)   : download from https://cran.r-project.org/bin/windows/Rtools/\n",
        "After installing, restart R and run user_run_file.R again."
      )
    }

  } else if (Sys.info()[["sysname"]] == "Darwin") {
    stop(
      "Xcode Command Line Tools are required to compile StarTime but were not found.\n\n",
      "Run the following command in your Terminal, then restart R:\n",
      "  xcode-select --install"
    )

  } else {
    stop(
      "C++ compilation tools are required to compile StarTime but were not found.\n\n",
      "Install them from your terminal, then restart R:\n",
      "  Ubuntu/Debian : sudo apt-get install r-base-dev\n",
      "  Fedora/CentOS : sudo yum install R-devel"
    )
  }
}

# --- Step 2: Ensure pak is available -----------------------------------------

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

# --- Step 3: Package list and pinned versions --------------------------------
#
# Versions match the environment used to produce the paper results.
# pak installs from the CRAN archive when the pinned version differs from
# what is currently available on CRAN.

packages <- c("glmnet", "doParallel", "doRNG", "foreach", "ggplot2",
              "dplyr", "tidyr", "patchwork", "readr", "StarTime",
              "RColorBrewer", "kableExtra", "MCS", "midasml", "stringr",
              "knitr", "forecast", "midasr")

versions <- c(
  "glmnet"       = "4.1-8",
  "doParallel"   = "1.0.17",
  "doRNG"        = "1.8.6.2",
  "foreach"      = "1.5.2",
  "ggplot2"      = "4.0.1",
  "dplyr"        = "1.1.4",
  "tidyr"        = "1.3.1",
  "patchwork"    = "1.3.2",
  "readr"        = "2.1.5",
  "RColorBrewer" = "1.1-3",
  "kableExtra"   = "1.4.0",
  "MCS"          = "0.1.3",
  "midasml"      = "0.1.10",
  "stringr"      = "1.5.1",
  "knitr"        = "1.51",
  "forecast"     = "8.24.0",
  "midasr"       = "0.8"
)

# --- Step 4: Install / upgrade to pinned versions ----------------------------

pkgs_not_available <- c()

for (pkg in packages) {

  if (pkg == "StarTime") {
    if (!requireNamespace("StarTime", quietly = TRUE)) {
      cat("Installing StarTime from local tarball...\n")
      tryCatch(
        install.packages("setup/StarTime_1.0.tar.gz",
                         repos = NULL, type = "source"),
        error = function(e) {
          cat("StarTime install FAILED:", conditionMessage(e), "\n")
          pkgs_not_available <<- c(pkgs_not_available, "StarTime")
        }
      )
    }

  } else {
    # For CRAN packages: install if missing OR if the installed version
    # does not match the pinned version
    current_ver <- tryCatch(
      as.character(packageVersion(pkg)),
      error = function(e) NA_character_
    )
    target_ver  <- versions[[pkg]]
    needs_install <- is.na(current_ver) ||
                     (!is.na(target_ver) && current_ver != target_ver)

    if (needs_install) {
      pkg_spec <- if (!is.na(target_ver)) paste0(pkg, "@", target_ver) else pkg
      cat("Installing", pkg_spec, "\n")
      tryCatch(
        pak::pkg_install(pkg_spec, ask = FALSE),
        error = function(e) {
          cat(pkg, "install FAILED:", e$message, "\n")
          pkgs_not_available <<- c(pkgs_not_available, pkg)
        }
      )
    }
  }
}

# --- Step 5: Report any failures ---------------------------------------------

if (length(pkgs_not_available) > 0) {
  warning(
    "The following packages could not be installed: ",
    paste(pkgs_not_available, collapse = ", ")
  )
}

# Verify every required package is loadable
for (p in setdiff(packages, pkgs_not_available)) {
  if (!requireNamespace(p, quietly = TRUE))
    cat("WARNING: package not available after install:", p, "\n")
}
