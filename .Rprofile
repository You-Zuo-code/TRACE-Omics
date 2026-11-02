local({
  lib <- file.path(".library", paste(R.version$major, R.version$minor, R.version$platform, sep = "-"))
  if (dir.exists(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
})
