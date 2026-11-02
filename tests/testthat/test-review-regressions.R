test_that("binary pathway effects are independent of metadata order", {
  set.seed(1234)
  ids <- paste0("s", seq_len(12))
  meta <- data.frame(Subject = rep(ids, each = 2), Group = rep(rep(c("No", "Yes"), each = 6), each = 2),
                     row.names = paste0("sample", seq_len(24)))
  scores <- matrix(rnorm(4 * 24), 4, dimnames = list(paste0("path", 1:4), rownames(meta)))
  scores[, meta$Group == "Yes"] <- scores[, meta$Group == "Yes"] + 2
  d <- list(meta = meta, subject = "Subject")
  a <- byod_pathway_clinical(d, scores, "Group")
  d$meta <- d$meta[24:1, ]
  b <- byod_pathway_clinical(d, scores[, 24:1], "Group")
  expect_equal(a$tab[order(a$tab$pathway), ], b$tab[order(b$tab$pathway), ], ignore_attr = TRUE)
  expect_true(all(a$tab$tested == "Yes" & a$tab$reference == "No"))
  expect_true(all(a$tab$estimate > 0))
  d$clinical_contrasts <- list(Group = c("No", "Yes"))
  c <- byod_pathway_clinical(d, scores, "Group")
  expect_equal(c$tab$estimate, -a$tab$estimate)
  d$meta$Group[1] <- "No"
  bad <- byod_pathway_clinical(d, scores, "Group")
  expect_false(bad$ok)
})

test_that("between-subject pairing never guesses a metadata block", {
  set.seed(1234)
  meta <- expand.grid(Timepoint = c("t0", "t1"), Subject = paste0("s", 1:9), stringsAsFactors = FALSE)
  meta$Group <- rep(c(rep(c("A", "B"), 4), "A"), each = 2)
  meta$Pair <- rep(c(rep(paste0("pair", 1:4), each = 2), "incomplete"), each = 2)
  meta$Batch <- meta$Pair
  rownames(meta) <- paste0("sample", seq_len(nrow(meta)))
  d <- list(meta = meta, expr = matrix(rnorm(40 * nrow(meta)), 40,
    dimnames = list(paste0("gene", 1:40), rownames(meta))))
  refused <- byod_run_limma(d, "Group", "A", "B", "Subject", TRUE)
  expect_equal(nrow(refused), 0L)
  expect_match(attr(refused, "declined"), "explicitly declared")
  unpaired <- byod_run_limma(d, "Group", "A", "B", "Subject", FALSE)
  expect_equal(attr(unpaired, "n_units"), 9L)
  expect_equal(attr(unpaired, "n"), 18L)
  d$nest <- "Pair"
  paired <- byod_run_limma(d, "Group", "A", "B", "Subject", TRUE)
  expect_equal(attr(paired, "block"), "Pair")
  expect_equal(attr(paired, "n_units"), 4L)
  expect_equal(attr(paired, "n"), 16L)
  expect_equal(attr(paired, "excluded_profiles"), 2L)
  expect_equal(attr(paired, "excluded_subjects"), "s9")
  expect_equal(attr(byod_run_limma(d, "Group", "A", "B", NULL, TRUE), "declined"),
    "Paired analysis requires a declared subject ID for every selected sample. Declare the pairing or turn Paired off.")
})

test_that("constant-feature filtering is optional and recorded", {
  set.seed(1234)
  meta <- expand.grid(Timepoint = c("pre", "post", "other"), Subject = paste0("s", 1:8),
                      stringsAsFactors = FALSE)
  meta$Pair <- rep(rep(paste0("p", 1:4), each = 2), each = 3)
  rownames(meta) <- paste0("sample", seq_len(nrow(meta)))
  E <- matrix(rnorm(30 * nrow(meta)), 30,
              dimnames = list(paste0("g", 1:30), rownames(meta)))
  E[1, ] <- 0
  E[2, ] <- 5
  E[3, ] <- ifelse(meta$Timepoint == "other", seq_len(nrow(meta)), 0)
  E[4, ] <- seq_len(nrow(meta)) * 1e-10
  expect_false(byod_clean_defaults()$drop_constant)
  expect_identical(byod_clean_matrix(E)$M, E)
  cleaned <- byod_clean_matrix(E, list(drop_constant = TRUE))
  expect_equal(nrow(cleaned$M), 28L)
  expect_equal(cleaned$steps[[1]]$removed, 2L)
  expect_true("g4" %in% rownames(cleaned$M))
  for (mode in c("unpaired", "paired", "nested")) {
    d <- list(meta = meta, expr = E, tp_levels = c("pre", "post", "other"))
    if (mode == "nested") d$nest <- "Pair"
    run <- function(x, ...) suppressWarnings(byod_run_limma(x, "Timepoint", "pre", "post",
                         "Subject", paired = mode != "unpaired", ...))
    default <- run(d)
    expect_identical(default, run(d, filter_constant = FALSE))
    expect_equal(nrow(default), 30L)
    expect_false(attr(default, "filter_constant"))
    expect_length(attr(default, "removed_features"), 0L)
    filtered <- run(d, filter_constant = TRUE)
    expect_equal(nrow(filtered), 27L)
    expect_identical(attr(filtered, "removed_features"), c("g1", "g2", "g3"))
    manual <- d; manual$expr <- E[-(1:3), , drop = FALSE]
    reference <- run(manual)
    expect_equal(filtered$Feature, reference$Feature)
    expect_equal(filtered$P.Value, reference$P.Value)
    expect_equal(filtered$adj.P.Val, reference$adj.P.Val)
    expect_identical(d$expr, E)
    all_constant <- d; all_constant$expr[,] <- 0
    refused <- run(all_constant, filter_constant = TRUE)
    expect_match(attr(refused, "declined"), "Fewer than two features")
    expect_length(attr(refused, "removed_features"), 30L)
  }
  html <- htmltools::renderTags(source(file.path("..", "..", "ui.R"), local = new.env(), chdir = TRUE)$value)$html
  for (id in c("clean_const", "byod_de_filter_constant")) {
    tag <- regmatches(html, regexpr(paste0('<input[^>]*id="', id, '"[^>]*>'), html))
    expect_length(tag, 1L)
    expect_false(grepl("checked", tag, fixed = TRUE))
  }
})

test_that("the differential-page switch reaches the fit and result archive", {
  root <- normalizePath(if (file.exists("global.R")) "." else "../..")
  wd <- getwd(); setwd(root); on.exit(setwd(wd), add = TRUE)
  server <- source("server.R", local = TRUE)$value
  study <- STUDIES$nmp_hum
  study$layers <- study$layers["mRNA"]; study$layer_names <- "mRNA"
  d <- study$layers$mRNA
  set.seed(1234)
  d$expr <- matrix(rnorm(30 * ncol(d$expr)), 30,
    dimnames = list(paste0("test", 1:30), colnames(d$expr)))
  d$expr[1, ] <- 0; d$expr[2, ] <- 7
  study$layers$mRNA <- d
  shiny::testServer(server, {
    byod_study(study)
    session$setInputs(byod_layer = "mRNA", byod_de_var = "Timepoint",
      byod_de_g1 = "preNMP", byod_de_g2 = "6h", byod_de_paired = TRUE,
      byod_de_fdr = .05, byod_de_lfc = 0, byod_de_nlab = 0,
      byod_de_filter_constant = FALSE)
    suppressWarnings(session$setInputs(byod_de_run = 1))
    expect_equal(nrow(byod_de_res()), 30L)
    expect_false(byod_de_export()$result$contrast$filter_constant)
    session$setInputs(byod_de_filter_constant = TRUE)
    expect_false(byod_de_export()$result$contrast$filter_constant)
    session$setInputs(byod_de_run = 2)
    result <- byod_de_export()$result
    expect_equal(nrow(result$differential), 28L)
    expect_true(result$contrast$filter_constant)
    expect_equal(result$contrast$n_features_removed, 2L)
    expect_true(all(result$constant_features$removed))
    expect_identical(result$constant_features$Feature, c("test1", "test2"))
    expect_identical(byod_layer_d()$expr, d$expr)
  })
})

test_that("failed raster exports cannot lower another session's budget", {
  old <- options(byod.export.max_px = 1e7)
  on.exit(options(old), add = TRUE)
  scope <- new.env(parent = environment(byod_save_figure))
  save <- byod_save_figure
  environment(save) <- scope
  scope$ggsave <- function(...) stop("simulated device allocation failure")
  expect_error(save(ggplot2::ggplot(), tempfile(fileext = ".png"), w = 1, h = 1, dpi = 300),
               "simulated device allocation failure")
  expect_identical(getOption("byod.export.max_px"), 1e7)
  attempts <- 0L
  scope$ggsave <- function(...) { attempts <<- attempts + 1L; if (attempts == 1L) stop("allocation failure") }
  outcome <- save(ggplot2::ggplot(), tempfile(fileext = ".png"), w = 1, h = 1, dpi = 300)
  expect_equal(outcome$dpi, 150)
  expect_identical(getOption("byod.export.max_px"), 1e7)
  second <- save(ggplot2::ggplot(), tempfile(fileext = ".png"), w = 1, h = 1, dpi = 300)
  expect_equal(second$dpi, 300)
})

test_that("conditional panels stay visible on Shiny 1.14 and later", {
  css <- paste(readLines("../../www/theme.css", warn = FALSE), collapse = "\n")
  expect_match(css, ".shiny-panel-conditional.shiny-conditional--shown", fixed = TRUE)
  expect_match(css, "display: contents !important", fixed = TRUE)
})

test_that("interactive volcano and baseline views convert to plotly", {
  pig <- STUDIES$nmp_pig$layers$mRNA
  multi <- suppressWarnings(byod_multi_contrast(pig, "in_vivo", c("1h", "6h")))
  expect_s3_class(suppressWarnings(byod_multi_volcano_plot(multi)), "plotly")
  baseline <- suppressWarnings(byod_baseline_sensitivity(pig, "in_vivo", "preNMP", "6h"))
  expect_true(isTRUE(baseline$ok))
  expect_s3_class(suppressWarnings(byod_baseline_plot(baseline)), "plotly")
})
