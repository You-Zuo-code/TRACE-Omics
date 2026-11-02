test_that("annotation guidance distinguishes supported resources and advanced choices", {
  expect_match(annotation_resource_hint("Human"), "no conversion table")
  expect_match(annotation_resource_hint("Human", "uniprot"), "org.Hs.eg.db")
  for (sp in c("Pig", "Mouse", "Rat"))
    expect_match(annotation_resource_hint(sp, "ensembl"), paste0(sp, " to human orthologs"))
  expect_match(annotation_resource_hint("Other"), "No built-in mapping")
  expect_match(annotation_resource_hint("Mouse", "uniprot"), "No built-in mapping")
  expect_match(annotation_resource_hint("Mouse", id_species = "human"), "no conversion table")
  expect_match(annotation_resource_hint("Other", mode = "none"), "no human mapping")
  expect_match(annotation_resource_hint("Human", mapping_choice = "custom"), "built-in mapping is not used")
})

test_that("the upload interface offers two own-cohort examples and explicit custom switches", {
  root <- normalizePath(if (file.exists("global.R")) "." else "../..")
  wd <- getwd(); setwd(root); on.exit(setwd(wd), add = TRUE)
  ui <- source("ui.R", local = new.env())$value
  html <- htmltools::renderTags(ui)$html
  expect_match(html, 'id="byod_example_human_rna"', fixed = TRUE)
  expect_match(html, 'id="byod_example_human_multi"', fixed = TRUE)
  expect_false(grepl("byod_demo|T-cell|GSE126805", html))
  expect_match(html, "Two views of the same 24-donor cohort", fixed = TRUE)
  for (i in 1:3) {
    expect_true(grepl(sprintf("input.byod_mapmode%d == &#39;custom&#39;", i), html, fixed = TRUE))
    expect_true(grepl(sprintf("input.byod_usegmt%d == true", i), html, fixed = TRUE))
  }
})

test_that("human examples retain curated values and work independently in workspace slots", {
  root <- normalizePath(if (file.exists("global.R")) "." else "../..")
  wd <- getwd(); setwd(root); on.exit(setwd(wd), add = TRUE)
  server <- source("server.R", local = TRUE)$value
  original <- STUDIES$nmp_hum
  shiny::testServer(server, {
    session$setInputs(byod_target = "ds1", byod_organism = "Mouse", byod_norm1 = "Log2",
      byod_mapmode1 = "custom", byod_usegmt1 = TRUE)
    session$setInputs(byod_example_human_rna = 1); session$flushReact()
    a <- byod_ws()$ds1$study
    expect_false(is.null(a)); expect_identical(a$layer_names, "mRNA")
    expect_identical(a$organism, "Human"); expect_identical(a$reference_key, "nmp_hum")
    expect_identical(a$subject, "Kidney_ID"); expect_length(a$tp_levels, 5L)
    expect_equal(length(unique(a$meta$Kidney_ID)), 24L)
    expect_identical(a$layers$mRNA$expr, original$layers$mRNA$expr)
    expect_identical(a$layers$mRNA$unit, "log2(TPM+1)")
    expect_identical(a$layers$mRNA$pca, original$layers$mRNA$pca)
    expect_length(a$matched_samples, 0L)
    expect_length(a$shared_genes, 0L)
    expect_true(a$layers$mRNA$genesets_apply)
    expect_true(is_uploaded()); expect_false(is.null(output$upload_built_summary))
    session$setInputs(byod_addslot = 1, byod_target = "ds2")
    session$setInputs(byod_example_human_multi = 1); session$flushReact()
    b <- byod_ws()$ds2$study
    expect_identical(b$layer_names, c("mRNA", "Protein"))
    expect_identical(b$meta, original$meta)
    expect_identical(b$matched_samples, original$matched_samples)
    expect_identical(b$shared_genes, original$shared_genes)
    expect_gt(length(b$matched_samples), 0L)
    for (ln in b$layer_names) {
      expect_identical(b$layers[[ln]]$expr, original$layers[[ln]]$expr)
      expect_identical(b$layers[[ln]]$unit, original$layers[[ln]]$unit)
      exported <- b$layers[[ln]]$input_expression
      expect_identical(exported$feature, rownames(original$layers[[ln]]$expr))
      expect_identical(names(exported)[-1], colnames(original$layers[[ln]]$expr))
      expect_identical(as.numeric(as.matrix(exported[, -1])), as.numeric(original$layers[[ln]]$expr))
    }
    expect_identical(byod_ws()$ds1$study$dataset_id, a$dataset_id)
    expect_false(identical(a$dataset_id, b$dataset_id))
    expect_identical(STUDIES$nmp_hum, original)
    expect_false(is.null(output$upload_built_summary))
    expect_match(output$analyse_source_card$html, "not a new cohort")
    expect_match(output$csr_upload_hint$html, "not independent studies")
    archive <- output$byod_annotation_dl
    entries <- unzip(archive, list = TRUE)$Name
    expect_true(all(c("annotation/layer_01/input_expression.csv.gz", "annotation/layer_02/input_expression.csv.gz",
      "annotation/input_metadata.csv") %in% entries))
  })
})
