test_that("input paths with spaces are files, never shell commands", {
  root <- tempfile("input with spaces ")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  path <- file.path(root, "matrix values.csv")
  writeLines(c("feature,sample1", "001,2"), path)
  x <- byod_read(path)
  expect_identical(x[[1]], "001")
  expect_equal(x[[2]], 2)
  expect_error(byod_read(file.path(root, "missing matrix.csv")),
               "Input file does not exist", fixed = TRUE)
})

test_that("gzip inputs preserve feature IDs with or without a file extension", {
  expect_true(requireNamespace("R.utils", quietly = TRUE))
  root <- tempfile("gzip-input-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  expected <- data.frame(feature = c("001", "002"), sample1 = c(2, 3))
  for (name in c("matrix.csv.gz", "uploaded_file")) {
    path <- file.path(root, name)
    con <- gzfile(path, "wt")
    tryCatch(write.csv(expected, con, row.names = FALSE), finally = close(con))
    expect_equal(byod_read(path), expected)
  }
})
