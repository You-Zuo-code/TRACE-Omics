test_that("reference studies retain their nesting and clinical variables", {
  pig <- STUDIES$nmp_pig
  human <- STUDIES$nmp_hum
  expect_equal(length(unique(pig$layers$mRNA$meta$Pig_Number)), 30)
  expect_equal(length(unique(human$layers$mRNA$meta$Kidney_ID)), 24)
  expect_gt(length(human$clinical_num), 10)
  expect_gt(length(pig$matched_samples), 50)
  expect_identical(pig$layers$mRNA$unit, "log2(TPM+1)")
})

test_that("porcine treatment comparison uses the animal as its unit", {
  pig <- STUDIES$nmp_pig
  result <- byod_run_limma(pig$layers$mRNA, "WI_status", "NWI", "WI",
                           subject_col = pig$subject, paired = TRUE)
  expect_gt(nrow(result), 1000)
  expect_identical(attr(result, "unit"), "Pig_Number")
  expect_equal(attr(result, "n_units"), 30)
  expect_true(all(c("logFC", "P.Value", "adj.P.Val") %in% names(result)))
})

test_that("human clinical grouping uses one observation per donor", {
  human <- STUDIES$nmp_hum
  d <- human$layers$mRNA
  groups <- levels(factor(d$meta[["Donor sex"]]))
  result <- byod_run_limma(d, "Donor sex", groups[1], groups[2],
                           subject_col = human$subject, paired = FALSE)
  expect_identical(attr(result, "unit"), "Kidney_ID")
  expect_lte(attr(result, "n_units"), 24)
  expect_lt(attr(result, "n_units"), ncol(d$expr))
})
