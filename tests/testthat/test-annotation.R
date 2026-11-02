annotation_fixture <- function(ids, organism, config = NULL, seed = 42L) {
  set.seed(seed)
  meta <- expand.grid(Timepoint = c("pre", "post"), Subject = paste0("u", 1:8), stringsAsFactors = FALSE)
  rownames(meta) <- paste0("sample", seq_len(nrow(meta)))
  meta$Timepoint <- factor(meta$Timepoint, levels = c("pre", "post"))
  x <- matrix(rnorm(length(ids) * nrow(meta)), length(ids), dimnames = list(ids, rownames(meta)))
  x[seq_len(min(20, nrow(x))), meta$Timepoint == "post"] <- x[seq_len(min(20, nrow(x))), meta$Timepoint == "post"] + 2
  list(expr = x, meta = meta, organism = organism, species = "test", subject = "Subject",
    unit = "log2 intensity", tp_levels = c("pre", "post"), annotation = annotation_build(ids, organism, config))
}

test_that("known human, pig, mouse and rat IDs map without rewriting input", {
  expect_null(annotation_organism(NA_character_)); expect_null(annotation_organism(list()))
  expect_error(annotation_organism(c("Pig", "Human")), "one source organism")
  for (id_type in c("symbol", "ensembl", "entrez")) {
    ids <- c(symbol = "Trp53", ensembl = "ENSMUSG00000059552.9", entrez = "22059")
    a <- annotation_build(ids[[id_type]], "Mouse", list(id_type = id_type))
    expect_true(any(a$mapping$status == "mapped" & a$mapping$target_id == "TP53"))
    expect_identical(a$organism, "Mouse")
  }
  expect_equal(annotation_build("ENSG00000141510.17", "Human", list(id_type = "ensembl"))$mapping$target_id, "TP53")
  expect_equal(annotation_build("7157", "Human", list(id_type = "entrez"))$mapping$target_id, "TP53")
  expect_equal(annotation_build("P04637", "Human", list(id_type = "uniprot"))$mapping$target_id, "TP53")
  for (sp in c("Pig", "Rat")) {
    z <- ANNOTATION_ORTHOLOGS
    z <- z[z$organism == sp & z$human_symbol == "TP53" & z$relationship == "one_to_one_in_HCOP_snapshot", ]
    expect_gt(nrow(z), 0)
    a <- annotation_build(z$ensembl[1], sp, list(id_type = "ensembl"))
    expect_true(any(a$mapping$target_id == "TP53" & a$mapping$status == "mapped"))
  }
  expect_true(all(annotation_build("22059", "Pig", list(id_type = "entrez"))$mapping$status != "mapped"))
  expect_true(all(annotation_build("P04637", "Mouse", list(id_type = "uniprot"))$mapping$status == "unmapped"))
  expect_error(byod_map_ids(matrix(1, 1, 1), "entrez_hs"), "explicit organism")
})

test_that("strict orthology never falls back and case is not an ortholog map", {
  p <- STUDIES$nmp_pig$layers$mRNA; h <- STUDIES$nmp_hum$layers$mRNA
  expect_equal(length(cs_feature_pairs(p, h)$features), 14287L)
  cls <- vapply(ORTHOLOG_MAP, function(z) z$class, character(1))
  good <- head(intersect(names(cls)[cls == "one_to_one"], intersect(rownames(p$expr), rownames(h$expr))), 20)
  bad <- head(intersect(names(cls)[cls != "one_to_one"], intersect(rownames(p$expr), rownames(h$expr))), 3)
  p$expr <- p$expr[c(good, bad), , drop = FALSE]; h$expr <- h$expr[c(good, bad), , drop = FALSE]
  expect_setequal(cs_feature_pairs(p, h)$features, good)
  fake <- paste0("not_a_gene_", 1:220)
  a <- annotation_fixture(fake, "Mouse"); b <- annotation_fixture(toupper(fake), "Human")
  expect_length(cs_feature_pairs(a, b)$features, 0)
  a <- annotation_fixture(c("Trp53", "Gapdh"), "Mouse")
  b <- annotation_fixture(c("TP53", "GAPDH"), "Human")
  expect_setequal(cs_feature_pairs(a, b)$features, "TP53")
  expect_identical(unname(cs_feature_pairs(a, b)$keyA), "Trp53")
  expect_true(all(annotation_build("Gapdh", "Mouse")$mapping$status == "ambiguous"))
  expect_setequal(cs_feature_pairs(a, a)$features, c("Trp53", "Gapdh"))
  a$organism <- "Arabidopsis thaliana"; a$annotation <- annotation_build(rownames(a$expr), a$organism)
  expect_length(cs_feature_pairs(a, b)$features, 0)
})

test_that("human-ID aliases cannot count the same gene twice", {
  a <- annotation_build(c("ENSG00000141510.1", "ENSG00000141510.2"), "Human", list(id_type = "ensembl"))
  expect_true(all(a$mapping$status == "duplicate_target"))
  d <- annotation_fixture(c("ENSG00000141510.1", "ENSG00000141510.2"), "Human", list(id_type = "ensembl"))
  expect_length(annotation_keys(d), 0)
})

test_that("the fallback reader preserves source identifiers and numeric metadata", {
  reader <- byod_read
  environment(reader) <- list2env(list(requireNamespace = function(...) FALSE), parent = environment(byod_read))
  p <- tempfile(fileext = ".dat"); on.exit(unlink(p), add = TRUE)
  writeLines(c("sample\tage\tgroup\tvalue", "001\t40\tA\t2.5", "002\t50\tB\t3.5"), p)
  x <- reader(p)
  expect_identical(x$sample, c("001", "002")); expect_true(is.numeric(x$age)); expect_true(is.numeric(x$value))
  gz <- tempfile(fileext = ".gz"); on.exit(unlink(gz), add = TRUE)
  con <- gzfile(gz, "wt"); writeLines(c("sample,age,group", "001,40,A", "002,50,B"), con); close(con)
  expect_equal(reader(gz)$sample, c("001", "002"))
})

test_that("GMT and custom mapping validate complete files and remain layer-local", {
  tmp <- tempfile("annotation-"); dir.create(tmp); on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  mp <- file.path(tmp, "map.csv"); gp <- file.path(tmp, "sets.gmt")
  write.csv(data.frame(source_id = c("p1", "p2", "p3", "not_uploaded"), target_id = c("TP53", "EGFR", "IL6", "IL6")), mp, row.names = FALSE)
  expect_error(annotation_custom_map(mp, ""), "source and version")
  m <- annotation_custom_map(mp, "manual fixture v1")
  a <- annotation_build(c("p1", "p2", "p3"), "Test species", list(mode = "human", mapping = m))
  expect_equal(sum(a$mapping$status == "mapped"), 2)
  expect_identical(a$mapping$status[a$mapping$source_id == "p3"], "ambiguous")
  writeLines(c("path_A\ttest\tp1\tp2\tp3", "path_B\ttest\tp4\tp5\tp6"), gp)
  g <- annotation_read_gmt(gp)
  d <- annotation_fixture(paste0("p", 1:20), "Test species", list(mode = "none", gmt = g))
  expect_identical(annotation_collections(d), "Custom GMT")
  expect_length(annotation_sets(d, "Hallmark"), 0)
  expect_identical(annotation_sets(d, "Custom GMT")$path_A, c("p1", "p2", "p3"))
  expect_gt(nrow(run_ora(c("p1", "p2", "p3"), "Custom GMT", rownames(d$expr), d)), 0)
  expect_equal(nrow(run_gsea(setNames(seq_len(20), rownames(d$expr)), "Hallmark", d = d)), 0)
  expect_true(annotation_compatible(d, d, "Custom GMT"))
  other <- d; other$annotation$gmt$hash <- "different"
  expect_false(annotation_compatible(d, other, "Custom GMT"))
  expect_false(xs_pathway_compare(d, other, xs_spec("Timepoint", "post", "pre"), xs_spec("Timepoint", "post", "pre"), "Custom GMT")$ok)
  writeLines(c("same_name\ttest\ta", "same name\ttest\tb"), gp)
  expect_error(annotation_read_gmt(gp), "unique")
  writeLines("not a GMT", gp); expect_error(annotation_read_gmt(gp), "tab")
})

test_that("mapped pathway tests equal the same matrix with human row names", {
  z <- ANNOTATION_ORTHOLOGS
  z <- z[z$organism == "Mouse" & z$relationship == "one_to_one_in_HCOP_snapshot" & z$support_n >= 3 & !is.na(z$symbol), ]
  ids <- unique(z$symbol)
  a <- annotation_fixture(ids, "Mouse")
  keys <- annotation_keys(a)
  keys <- keys[names(keys) %in% unique(unlist(GENESETS$Hallmark))]
  keys <- head(keys, 1500)
  a$expr <- a$expr[unname(keys), , drop = FALSE]
  a$annotation <- annotation_build(rownames(a$expr), "Mouse")
  b <- a; b$expr <- a$expr; rownames(b$expr) <- names(keys); b$organism <- "Human"
  b$annotation <- annotation_build(rownames(b$expr), "Human")
  stat <- seq_len(nrow(a$expr)) / nrow(a$expr) * 4 - 2
  ra <- run_gsea(setNames(stat, rownames(a$expr)), "Hallmark", d = a)
  rb <- run_gsea(setNames(stat, rownames(b$expr)), "Hallmark", d = b)
  expect_gt(nrow(ra), 5)
  expect_equal(ra, rb)
  oa <- run_ora(head(rownames(a$expr), 100), "Hallmark", rownames(a$expr), a)
  ob <- run_ora(head(rownames(b$expr), 100), "Hallmark", rownames(b$expr), b)
  expect_equal(oa, ob)
  ca <- byod_gsea_curve(setNames(stat, rownames(a$expr)), "Hallmark", ra$pathway[1], a)
  cb <- byod_gsea_curve(setNames(stat, rownames(b$expr)), "Hallmark", rb$pathway[1], b)
  expect_true(ca$ok); expect_equal(ca$curve, cb$curve)
  pa <- byod_pathway_activity(a, "Hallmark", "ssGSEA")
  pb <- byod_pathway_activity(b, "Hallmark", "ssGSEA")
  expect_true(pa$ok); expect_identical(pa$executed_method, "ssGSEA")
  expect_equal(dimnames(pa$mat), dimnames(pb$mat)); expect_equal(as.numeric(pa$mat), as.numeric(pb$mat))
  spec <- xs_spec("Timepoint", "post", "pre")
  dc <- xs_de_compare(a, b, spec, spec, n_perm = 0L)
  expect_true(dc$ok); expect_equal(dc$n_shared, nrow(a$expr)); expect_equal(dc$r, 1)
  pw <- xs_pathway_compare(a, b, spec, spec)
  expect_true(pw$ok); expect_equal(pw$r, 1)
  expect_setequal(rownames(a$expr), unname(keys))
})
