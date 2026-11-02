test_that("current installation accepts newer R; only opt-in reproduction uses the old release", {
  root <- if (file.exists("DESCRIPTION")) "." else "../.."
  source(file.path(root, "tools", "install_helpers.R"), local = TRUE)
  for (v in c("4.5.0", "4.5.1", "4.6.0", "4.7.0", "5.0.0"))
    expect_identical(trace_install_mode(r = package_version(v)), "current")
  expect_identical(trace_install_mode("--locked", r = package_version("4.5.1")), "locked")
  expect_error(trace_install_mode("--locked", r = package_version("4.6.0")), "default installer")
  expect_error(trace_install_mode(mode = "invalid"), "mode")
  expect_error(trace_install_mode(r = package_version("4.3.0")), "4.5 or later")
  dependencies <- trace_dependencies(file.path(root, "DESCRIPTION"))
  expect_true(all(c("shiny", "systemfonts", "mixOmics", "org.Hs.eg.db") %in% dependencies))
  expect_false(anyDuplicated(dependencies) > 0)
})

test_that("cross-study Methods report both declared designs without NA placeholders", {
  p <- list(study_A = "TIV 2009", study_B = "TIV 2011", A_contrast = "D7 vs D0",
            B_contrast = "D7 vs D0", A_pairing = "paired within subject", B_pairing = "paired within subject",
            A_n_units = 28L, B_n_units = 71L, A_unit = "Subject", B_unit = "Subject",
            shared_features = 54715L, feature_rule = "exact IDs", fdr = .05,
            permutations = 0, outcome = "completed")
  text <- byod_methods_text(list(list(module = "cross-dataset differential concordance", params = p)))
  expect_match(text, "TIV 2009: D7 vs D0", fixed = TRUE)
  expect_match(text, "n=28 Subject", fixed = TRUE)
  expect_match(text, "n=71 Subject", fixed = TRUE)
  expect_match(text, "54715 shared features", fixed = TRUE)
  expect_false(grepl("on NA|using NA", text))
  extra <- byod_methods_text(list(list(module = "Module score", params = list(outcome = "completed"))))
  expect_false(grepl("on NA|using NA", extra))
  p$outcome <- "declined"; p$reason <- "insufficient units"
  declined <- byod_methods_text(list(list(module = "cross-dataset differential concordance", params = p)))
  expect_match(declined, "declined")
  expect_false(grepl("were fitted separately", declined, fixed = TRUE))
})

test_that("comparison archive keys are not counted as molecular layers", {
  inventory <- data.frame(dataset = c("study1", "study1", BUNDLE_PAIR_DATASET),
                          layer = c("RNA", "Protein", "study1--vs--study2"))
  expect_identical(bundle_recorded_layers(inventory), c("RNA", "Protein"))
  expect_length(bundle_recorded_layers(inventory[3, ]), 0L)
})

test_that("empty analysis panels render their message without trace warnings", {
  expect_warning(p <- plotly::plotly_build(emptyplot("No results yet")), NA)
  expect_match(p$x$layout$annotations[[1]]$text, "No results yet", fixed = TRUE)
})
