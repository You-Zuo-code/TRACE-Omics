local({
  if (!file.exists("DESCRIPTION") || !file.exists("tools/install_helpers.R"))
    stop("Run install.R from the TRACE-Omics directory.")
  source("tools/install_helpers.R", local = TRUE)
  mode <- trace_install_mode(commandArgs(TRUE), Sys.getenv("TRACE_INSTALL_MODE", "current"))
  lib <- file.path(".library", paste(R.version$major, R.version$minor, R.version$platform, sep = "-"))
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  lib <- normalizePath(lib)
  .libPaths(c(lib, .libPaths()))
  options(repos = c(CRAN = "https://cloud.r-project.org"),
          BioC_mirror = Sys.getenv("TRACE_BIOC_MIRROR", "https://bioconductor.posit.co"),
          timeout = max(1200, getOption("timeout", 60)))
  message("Installing for R ", getRversion(), " (", mode, " mode) into ", lib)
  if (mode == "locked") {
    if (!file.exists("renv.lock")) stop("renv.lock is missing.")
    if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", lib = lib)
    locked <- renv::lockfile_read("renv.lock")$Packages
    renv::restore(project = getwd(), library = lib, lockfile = "renv.lock",
                  repos = c(CRAN = "https://packagemanager.posit.co/cran/2026-02-21",
                            CRAN_fallback = "https://cloud.r-project.org"), prompt = FALSE)
    installed <- installed.packages(lib.loc = lib)
    stock <- installed.packages(priority = c("base", "recommended"), lib.loc = .Library)
    installed <- rbind(installed, stock[setdiff(rownames(stock), rownames(installed)), , drop = FALSE])
    wrong <- names(locked)[vapply(names(locked), function(p)
      !p %in% rownames(installed) || installed[p, "Version"] != locked[[p]]$Version, logical(1))]
    if (length(wrong)) stop("Lockfile verification failed: ", paste(wrong, collapse = ", "))
  } else {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib)
    message("Bioconductor release: ", BiocManager::version())
    BiocManager::install(trace_dependencies(), lib = lib, ask = FALSE, update = FALSE)
  }
  missing <- trace_dependencies()[!vapply(trace_dependencies(), requireNamespace,
                                         logical(1), quietly = TRUE)]
  if (length(missing)) stop("Packages could not be loaded: ", paste(missing, collapse = ", "),
                             ". See the installation messages above.")
  packages <- installed.packages(lib.loc = .libPaths())
  write.csv(packages[, c("Package", "Version", "LibPath"), drop = FALSE],
            file.path(lib, "installed-packages.csv"), row.names = FALSE)
  message("Installation complete. Start with: Rscript run.R")
  message("Check this environment with: Rscript tests/run.R")
})
