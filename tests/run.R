if (!file.exists("global.R")) stop("Run tests/run.R from the TRACE-Omics directory.")
lib <- file.path(".library", paste(R.version$major, R.version$minor, R.version$platform, sep = "-"))
if (dir.exists(lib)) .libPaths(c(normalizePath(lib), .libPaths()))
library(testthat)
source("global.R", local = FALSE)
ui <- source("ui.R", local = new.env())$value
stopifnot(inherits(shiny::shinyAppDir("."), "shiny.appobj"))
html <- htmltools::renderTags(ui)
stopifnot(grepl("Longitudinal multi-omics analysis", html$html, fixed = TRUE),
          grepl("theme.css", html$head, fixed = TRUE))
testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
