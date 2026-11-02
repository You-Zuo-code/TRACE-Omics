test_that("custom species uploads export the original input, mapping and GMT", {
  root <- normalizePath(if (file.exists("global.R")) "." else "../..")
  wd <- getwd(); setwd(root); on.exit(setwd(wd), add = TRUE)
  tmp <- tempfile("species-upload-"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  ids <- head(unique(unlist(GENESETS$Hallmark)), 100)
  set.seed(3)
  meta <- expand.grid(Timepoint = c("pre", "post"), Subject = paste0("u", 1:8), stringsAsFactors = FALSE)
  meta$sample <- paste0("s", seq_len(nrow(meta)))
  x <- matrix(rnorm(100 * nrow(meta), 8, .4), 100, dimnames = list(NULL, meta$sample))
  x[1:20, meta$Timepoint == "post"] <- x[1:20, meta$Timepoint == "post"] + 2
  raw <- data.frame(Feature = paste0("custom_", 1:100), x, check.names = FALSE)
  hu <- raw; hu$Feature <- ids
  write.csv(raw, file.path(tmp, "native.csv"), row.names = FALSE)
  write.csv(hu, file.path(tmp, "human.csv"), row.names = FALSE)
  write.csv(meta[, c("sample", "Timepoint", "Subject")], file.path(tmp, "metadata.csv"), row.names = FALSE)
  write.csv(data.frame(source_id = raw$Feature, target_id = ids), file.path(tmp, "map.csv"), row.names = FALSE)
  writeLines(paste(c("custom_path", "fixture", raw$Feature[1:20]), collapse = "\t"), file.path(tmp, "sets.gmt"))
  writeLines("bad GMT", file.path(tmp, "bad.gmt"))
  fi <- function(f) { p <- file.path(tmp, f); data.frame(name = f, size = file.info(p)$size, type = "text/plain", datapath = p) }
  server <- source("server.R", local = TRUE)$value
  shiny::testServer(server, {
    session$setInputs(byod_target = "ds1", byod_expr1 = fi("native.csv"), byod_meta = fi("metadata.csv"),
      byod_type1 = "transcriptomics", byod_name1 = "RNA", byod_unit1 = "log2 intensity", byod_norm1 = "None",
      byod_idtype1 = "other", byod_idspace1 = "native", byod_annotation1 = "human",
      byod_mapping1 = fi("map.csv"), byod_mapsource1 = "test mapping v1", byod_gmt1 = fi("sets.gmt"),
      byod_mapmode1 = "custom", byod_usegmt1 = TRUE,
      byod_organism = "Other", byod_custom_organism = "Testus species (999999)",
      byod_tpcol = "Timepoint", byod_subj = "Subject", byod_studyname = "Native study")
    session$setInputs(byod_build = 1); session$flushReact()
    a <- byod_ws()$ds1$study; expect_false(is.null(a))
    expect_identical(a$organism, "Testus species (999999)")
    expect_identical(rownames(a$layers$RNA$expr), raw$Feature)
    expect_equal(as.numeric(a$layers$RNA$expr), as.numeric(as.matrix(raw[, -1])))
    expect_equal(a$layers$RNA$input_expression, raw)
    expect_equal(sum(a$layers$RNA$annotation$mapping$status == "mapped"), 100L)
    expect_true("Custom GMT" %in% annotation_collections(a$layers$RNA))
    downloaded <- output$byod_annotation_dl
    entries <- unzip(downloaded, list = TRUE)$Name
    expect_true(all(c("annotation/layer_01/mapping.csv", "annotation/layer_01/custom.gmt",
      "annotation/layer_01/input_expression.csv.gz", "annotation/input_metadata.csv", "annotation/resources.json") %in% entries))
    dest <- file.path(tmp, "extracted"); unzip(downloaded, exdir = dest)
    expect_equal(byod_read(file.path(dest, "annotation/layer_01/input_expression.csv.gz")), raw)
    expect_equal(nrow(read.csv(file.path(dest, "annotation/layer_01/mapping.csv"))), 100L)
    session$setInputs(byod_addslot = 1, byod_target = "ds2", byod_expr1 = fi("human.csv"),
      byod_mapmode1 = "builtin", byod_usegmt1 = FALSE,
      byod_idtype1 = "symbol", byod_organism = "Human", byod_studyname = "Human study")
    session$setInputs(byod_build = 2); session$flushReact()
    b <- byod_ws()$ds2$study; expect_false(is.null(b))
    expect_null(b$layers$RNA$annotation$config$mapping)
    expect_null(b$layers$RNA$annotation$gmt)
    expect_equal(sum(b$layers$RNA$annotation$mapping$status == "mapped"), 100L)
    session$setInputs(csr_A = "ds1", csr_B = "ds2", csr_layA = "RNA", csr_layB = "RNA",
      xs_varA = "Timepoint", xs_varB = "Timepoint", xs_testA = "post", xs_refA = "pre",
      xs_testB = "post", xs_refB = "pre", xs_paired = TRUE, xs_perm = FALSE, xs_fdr = .05)
    session$setInputs(xs_de_run = 1); session$flushReact()
    r <- xs_de_res(); expect_true(r$ok); expect_equal(r$n_shared, 100L); expect_equal(r$r, 1)
    expect_true(bundle_capture(session, "xs_de_sc_appr", xs_de_sc_export(), xs_de_scatter_gg(r),
      BUNDLE_PAIR_DATASET, bundle_pair_folder(a$dataset_id, b$dataset_id)))
    archive <- output$bundle_dl
    entries <- unzip(archive, list = TRUE)$Name
    expect_true("annotation/layer_01/mapping.csv" %in% entries)
    expect_true(any(grepl("^comparison_annotation/.+/layer_01/supplied_mapping.csv$", entries)))
    f <- entries[grepl("result_shared_features.csv$", entries)]
    expect_length(f, 1); expect_equal(nrow(read.csv(unz(archive, f))), 100L)
    old_id <- b$dataset_id
    session$setInputs(byod_usegmt1 = TRUE, byod_gmt1 = fi("bad.gmt"), byod_build = 3); session$flushReact()
    expect_identical(byod_ws()$ds2$study$dataset_id, old_id)
    expect_match(byod_msg()$txt, "GMT")
    session$setInputs(byod_gmt1 = NULL, byod_build = 4); session$flushReact()
    expect_identical(byod_ws()$ds2$study$dataset_id, old_id)
    expect_match(byod_msg()$txt, "upload a GMT")
    session$setInputs(byod_usegmt1 = FALSE, byod_mapmode1 = "custom", byod_mapping1 = NULL, byod_build = 5)
    session$flushReact()
    expect_identical(byod_ws()$ds2$study$dataset_id, old_id)
    expect_match(byod_msg()$txt, "upload the selected custom mapping")
    expect_identical(byod_ws()$ds1$study$dataset_id, a$dataset_id)
  })
})
