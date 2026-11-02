library(testthat)

seed_test_root <- normalizePath(if (file.exists("global.R")) "." else "../..")
if (!exists("trace_seed")) {
  seed_test_wd <- getwd(); setwd(seed_test_root)
  suppressWarnings(suppressMessages(source("global.R")))
  setwd(seed_test_wd)
}

test_that("seed validation and local RNG scopes preserve the caller's state", {
  expect_identical(trace_seed(), 1234L)
  expect_identical(trace_validate_seed(0), 0L)
  expect_identical(trace_validate_seed(2147483646), 2147483646L)
  for (bad in list(NULL, NA, NaN, Inf, -1, 1.5, 2147483647, c(1, 2), "abc"))
    expect_error(trace_validate_seed(bad), "Random seed must")
  draw <- function(seed, fail = FALSE) {
    restore <- trace_rng_scope(seed); on.exit(restore())
    x <- runif(10)
    if (fail) stop("deliberate test error")
    x
  }
  original <- .Random.seed
  expect_identical(draw(1234), draw(1234))
  expect_false(identical(draw(1234), draw(4321)))
  expect_error(draw(1234, TRUE), "deliberate")
  expect_identical(.Random.seed, original)
})

test_that("UMAP uses the selected seed instead of cached legacy coordinates", {
  pca <- outer(seq_len(36), seq_len(6), function(i, j) sin(i * j) + cos(i / (j + 1)))
  rownames(pca) <- paste0("s", seq_len(nrow(pca)))
  d <- list(pca = pca, umap = data.frame(UMAP1 = rep(0, 36), UMAP2 = rep(0, 36)),
            species = "seed fixture", omic = "RNA")
  before <- .Random.seed
  a <- trace_umap(d, 1234)
  expect_identical(a, trace_umap(d, 1234))
  expect_false(isTRUE(all.equal(a, trace_umap(d, 4321))))
  expect_false(identical(a, d$umap))
  expect_identical(rownames(a), rownames(pca))
  expect_identical(.Random.seed, before)
})

test_that("grouped DIABLO CV repeats are reproducible without duplicating every fold", {
  s <- STUDIES$nmp_pig
  ids <- unique(s$layers$mRNA$meta$Pig_Number)[seq_len(12)]
  for (ln in s$layer_names) {
    d <- s$layers[[ln]]
    use <- d$meta$Pig_Number %in% ids & as.character(d$meta$Timepoint) %in% c("in_vivo", "6h")
    d$meta <- d$meta[use, , drop = FALSE]
    d$expr <- d$expr[seq_len(80), use, drop = FALSE]
    s$layers[[ln]] <- d
  }
  s$matched_samples <- Reduce(intersect, lapply(s$layers, function(d) colnames(d$expr)))
  cv <- function(seed = 1234L) suppressMessages(byod_diablo_cv(s, "mRNA", "Protein",
    nfeat = 40, keepx = 5, ncomp = 1, folds = 3, repeats = 2, seed = seed))
  a <- cv(); b <- cv(); c <- cv(4321L)
  expect_true(a$ok); expect_true(b$ok); expect_true(c$ok)
  expect_identical(a$seed, 1234L)
  expect_identical(a$pred, b$pred)
  prediction_key <- function(x) paste(x$repeat_i, x$sample)
  expect_false(identical(a$pred$fold, c$pred$fold[match(prediction_key(a$pred), prediction_key(c$pred))]))
  folds <- unique(a$pred[, c("repeat_i", "unit", "fold")])
  expect_equal(nrow(folds), 2L * a$n_units)
  f1 <- folds[folds$repeat_i == 1, ]; f2 <- folds[folds$repeat_i == 2, ]
  expect_false(identical(f1$fold, f2$fold[match(f1$unit, f2$unit)]))
})

test_that("the page has one seed control and applying it invalidates old results", {
  wd <- getwd(); setwd(seed_test_root); on.exit(setwd(wd), add = TRUE)
  ui <- source("ui.R", local = new.env())$value
  html <- htmltools::renderTags(ui)$html
  expect_equal(lengths(regmatches(html, gregexpr('id="analysis_seed"', html, fixed = TRUE))), 1)
  expect_match(html, 'value="1234"', fixed = TRUE)
  expect_false(grepl('id="xs_seed"', html, fixed = TRUE))
  expect_false(grepl("analysis_seed_status|Active seed:|After editing, click Apply seed", html))
  dom <- XML::htmlParse(html, asText = TRUE)
  panel <- XML::getNodeSet(dom, '//*[@id="analysis_seed"]/ancestor::div[@data-value][1]')
  expect_length(panel, 1L)
  expect_identical(XML::xmlGetAttr(panel[[1]], "data-value"), "Your data")
  server <- source("server.R", local = TRUE)$value
  shiny::testServer(server, {
    expect_identical(trace_seed(), 1234L)
    session$setInputs(analysis_seed = 4321)
    expect_identical(trace_seed(), 1234L)
    session$setInputs(apply_seed = 1)
    expect_identical(trace_seed(), 4321L)
    generation <- reactiveVal(0L)
    probe <- resettable(reactive(input$seed_probe_run), generation,
      function() list(value = runif(5), seed = trace_seed()))
    session$setInputs(seed_probe_run = 1)
    first <- probe()
    expect_identical(first$seed, 4321L)
    session$setInputs(analysis_seed = 1234, apply_seed = 2)
    expect_error(probe(), class = "shiny.silent.error")
    session$setInputs(seed_probe_run = 2)
    expect_identical(probe()$seed, 1234L)
    session$setInputs(analysis_seed = 1.5, apply_seed = 3)
    expect_identical(trace_seed(), 1234L)
    log_event("seed-test", engine = "stats")
    expect_identical(tail(analysis_log(), 1)[[1]]$params$random_seed, 1234L)
  })
  shiny::testServer(server, { expect_identical(trace_seed(), 1234L) })
})

test_that("GSEA and permutation tests reproduce their results with the same seed", {
  wd <- getwd(); setwd(seed_test_root); on.exit(setwd(wd), add = TRUE)
  d <- STUDIES$nmp_hum$layers$mRNA
  ids <- head(intersect(unique(unlist(annotation_sets(d, "Hallmark"))), rownames(d$expr)), 500)
  ranked <- setNames(sin(seq_along(ids) * 1.7), ids)
  before <- .Random.seed
  g1 <- run_gsea(ranked, "Hallmark", seed = 1234L, d = d)
  g2 <- run_gsea(ranked, "Hallmark", seed = 1234L, d = d)
  expect_gt(nrow(g1), 0)
  expect_identical(g1, g2)
  expect_identical(.Random.seed, before)

  d$expr <- d$expr[seq_len(80), , drop = FALSE]
  spec <- xs_spec("Timepoint", "6h", "preNMP", paired = TRUE)
  compare <- function(seed) xs_de_compare(d, d, spec, spec,
    n_perm = 8L, max_exact = 0L, seed = seed)
  a <- compare(1234L); b <- compare(1234L); c <- compare(4321L)
  expect_true(a$ok); expect_true(a$perm$ok)
  expect_identical(a$perm, b$perm)
  expect_false(identical(a$perm, c$perm))
  expect_identical(.Random.seed, before)
})
