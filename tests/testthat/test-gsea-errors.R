test_that("GSEA execution failures are not described as negative findings", {
  testthat::local_mocked_bindings(
    fgseaMultilevel = function(...) stop("test engine unavailable"),
    .package = "fgsea")
  d <- STUDIES$nmp_hum$layers$mRNA
  ids <- intersect(unique(unlist(annotation_sets(d, "Hallmark"))), rownames(d$expr))
  result <- run_gsea(setNames(seq_along(ids), ids), "Hallmark", d = d)
  expect_equal(nrow(result), 0L)
  expect_match(attr(result, "declined"), "test engine unavailable")
  expect_false(grepl("significance", attr(result, "declined")))
})

test_that("an empty pathway selection returns a message instead of an error", {
  for (pathway in list(NULL, character(), NA_character_, "", c("A", "B"))) {
    result <- byod_gsea_curve(numeric(), "Hallmark", pathway)
    expect_false(result$ok)
    expect_match(result$msg, "Select one pathway")
  }
})
