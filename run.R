local({
  frames <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  file <- if (length(frames)) tail(frames, 1)[[1]] else if (length(arg)) sub("^--file=", "", arg[1]) else "run.R"
  root <- dirname(normalizePath(file, mustWork = TRUE))
  lib <- file.path(root, ".library", paste(R.version$major, R.version$minor, R.version$platform, sep = "-"))
  if (dir.exists(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
  if (!requireNamespace("shiny", quietly = TRUE)) stop("Run Rscript install.R first.")
  shiny::runApp(root, host = "127.0.0.1", port = as.integer(Sys.getenv("TRACE_PORT", "3838")),
                launch.browser = interactive())
})
