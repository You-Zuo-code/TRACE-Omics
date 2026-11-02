test_that("exports set fonts on themes and text geoms without changing source plots", {
  d <- data.frame(x = 1:3, y = c(2, 1, 3), label = c("A", "B", "C"))
  g <- ggplot(d, aes(x, y)) + geom_point() + geom_text(aes(label = label)) +
    ggrepel::geom_text_repel(aes(label = label)) + theme_minimal()
  original <- lapply(g$layers, function(x) x$aes_params)
  out <- byod_export_font(g, "Helvetica")
  expect_identical(out$theme$text$family, trace_font("Helvetica"))
  expect_identical(out$layers[[2]]$aes_params$family, trace_font("Helvetica"))
  expect_identical(out$layers[[3]]$aes_params$family, trace_font("Helvetica"))
  expect_identical(lapply(g$layers, function(x) x$aes_params), original)
  expect_identical(out$data, g$data)
  a <- ggplot_build(g)$data
  b <- ggplot_build(out)$data
  for (i in seq_along(a)) {
    keep <- setdiff(names(a[[i]]), "family")
    expect_equal(b[[i]][keep], a[[i]][keep])
  }
})

test_that("an explicit text font is retained and Arial can be selected", {
  g <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
    annotate("text", x = 1, y = 1, label = "custom", family = "Courier") +
    annotate("label", x = 1, y = 2, label = "default")
  out <- byod_export_font(g, "Arial")
  expect_identical(out$theme$text$family, trace_font("Arial"))
  expect_identical(out$layers[[1]]$aes_params$family, "Courier")
  expect_identical(out$layers[[2]]$aes_params$family, trace_font("Arial"))
  mapped <- ggplot(data.frame(x = 1, label = "mapped", font = "Arial"),
                   aes(x, x, label = label, family = font)) + geom_text()
  expect_null(byod_export_font(mapped)$layers[[1]]$aes_params$family)
})

test_that("PDF, SVG and raster exports use the common renderer", {
  g <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point() +
    labs(title = "TRACE-Omics font check", x = "Timepoint", y = "GSVA score")
  dir <- tempfile("export-font-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  for (fmt in c("pdf", "svg", "png")) {
    f <- file.path(dir, paste0("figure.", fmt))
    res <- byod_save_figure(g, f, fmt, 4, 3, 72)
    expect_true(res$ok)
    expect_true(file.exists(f))
    expect_identical(res$format, fmt)
  }
  svg <- paste(readLines(file.path(dir, "figure.svg"), warn = FALSE), collapse = "\n")
  expect_match(svg, BYOD_EXPORT_FONT, fixed = TRUE)
  expect_true(is.function(byod_pdf_device()))
})

test_that("CJK labels select an installed font with glyph coverage", {
  latin <- ggplot(data.frame(x = 1, label = "GFI1B"), aes(x, x)) +
    geom_text(aes(label = label))
  chinese <- latin + labs(title = "\u4f9b\u4f53\u5e74\u9f84")
  expect_false(byod_plot_has_cjk(latin))
  expect_true(byod_plot_has_cjk(chinese))
  categorical <- ggplot(data.frame(x = factor("\u4f9b\u4f53"), y = 1), aes(x, y)) + geom_point()
  expect_true(byod_plot_has_cjk(categorical))
  font <- tryCatch(trace_font("Helvetica", cjk = TRUE), error = function(e) NULL)
  if (!is.null(font)) {
    out <- byod_export_font(chinese, "Helvetica")
    expect_identical(out$theme$text$family, font)
    expect_identical(out$layers[[1]]$aes_params$family, font)
    expect_true(all(systemfonts::glyph_info(c("\u4f9b", "\u4f53", "\u5e74", "\u9f84"), family = font)$index > 0))
  }
  expect_error(trace_font("Helvetica", cjk = TRUE, available = character()), "Install Noto Sans CJK")
  expect_identical(trace_font("Helvetica", available = "Liberation Sans"), "Liberation Sans")
  expect_identical(trace_font("Arial", available = c("Arial", "Helvetica")), "Arial")
})

test_that("Quartz background repair preserves graphics and is idempotent", {
  skip_if_not(isTRUE(capabilities("aqua")))
  first_stream <- function(f) {
    b <- readBin(f, "raw", file.info(f)$size)
    i <- grepRaw("\nstream\n", b, fixed = TRUE)[1] + 8L
    j <- grepRaw("\nendstream", b, offset = i, fixed = TRUE)[1]
    memDecompress(b[seq.int(i, j - 1L)], "gzip", asChar = TRUE)
  }
  dir <- tempfile("pdf-background-"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  for (bg in c("white", "transparent", "#ddeeff")) {
    f <- file.path(dir, paste0(gsub("#", "", bg), ".pdf"))
    grDevices::quartz(type = "pdf", file = f, width = 8, height = 7, bg = bg)
    print(ggplot(data.frame(x = 1:3, y = c(2, 1, 3)), aes(x, y)) +
            geom_point() + labs(title = "Helvetica PDF check \u00b7 \u03c1 \u2265 0") + theme_minimal(base_family = "Helvetica"))
    dev.off()
    before <- first_stream(f)
    expect_true(byod_bound_quartz_background(f, 8, 7))
    after <- first_stream(f)
    expect_true(grepl("0 504.000000 m 576.000000 504.000000 l 576.000000 0 l 0 0 l h f", after, fixed = TRUE, useBytes = TRUE))
    rest <- function(s) charToRaw(sub("(?s)^q\\s+Q\\s+q\\s+.*?h f\\s+Q", "", s, perl = TRUE, useBytes = TRUE))
    expect_identical(rest(after), rest(before))
    hash <- unname(tools::md5sum(f))
    expect_false(byod_bound_quartz_background(f, 8, 7))
    expect_identical(unname(tools::md5sum(f)), hash)
    expect_error(byod_bound_quartz_background(f, 4, 3), "Cannot validate")
    expect_identical(unname(tools::md5sum(f)), hash)
  }
})

test_that("the common PDF renderer bounds its Quartz background before returning", {
  skip_if_not(identical(attr(byod_pdf_device(), "trace_pdf_backend"), "quartz"))
  dir <- tempfile("pdf-export-"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  f <- file.path(dir, "figure.pdf")
  g <- ggplot(data.frame(x = 1:3, y = 1:3), aes(x, y)) + geom_point()
  for (wh in list(c(8, 5), c(8, 7), c(5, 4), c(6, 5), c(8.333, 6.125))) {
    expect_true(byod_save_figure(g, f, "pdf", wh[1], wh[2])$ok)
    expect_false(byod_bound_quartz_background(f, wh[1], wh[2]))
  }
  chinese <- g + labs(title = "\u4f9b\u4f53\u5e74\u9f84 / Donor age \u00b7 \u03c1")
  expect_true(byod_save_figure(chinese, f, "pdf", 8, 5)$ok)
  expect_false(byod_bound_quartz_background(f, 8, 5))
  bad <- file.path(dir, "not-a-pdf.pdf")
  writeBin(charToRaw("not a PDF"), bad)
  expect_error(byod_bound_quartz_background(bad, 8, 5), "Cannot validate")
  expect_identical(readBin(bad, "raw", file.info(bad)$size), charToRaw("not a PDF"))
})
