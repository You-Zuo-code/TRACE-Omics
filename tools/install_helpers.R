trace_install_mode <- function(args = character(), mode = "current", r = getRversion()) {
  if ("--locked" %in% args) mode <- "locked"
  if (!mode %in% c("current", "locked")) stop("Installation mode must be current or locked.")
  if (r < "4.5.0") stop("TRACE-Omics requires R 4.5 or later.")
  if (mode == "locked" && (r < "4.5.0" || r >= "4.6.0"))
    stop("The optional lockfile uses R 4.5.x / Bioconductor 3.21. Use the default installer with this R version.")
  mode
}

trace_dependencies <- function(file = "DESCRIPTION") {
  x <- read.dcf(file, fields = "Imports")[[1]]
  trimws(gsub("\\s*\\([^)]*\\)", "", strsplit(x, ",")[[1]]))
}
