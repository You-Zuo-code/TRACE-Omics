byod_bound_quartz_background <- function(file, width, height) {
  # Quartz can write a page background with both dimensions multiplied by 72.
  # Only fresh, single-page Quartz exports are handled; fonts and other objects stay intact.
  fail <- function() stop("Cannot validate the Quartz PDF background. Export SVG instead.", call. = FALSE)
  b <- readBin(file, "raw", n = file.info(file)$size)
  hits <- grepRaw("startxref", b, fixed = TRUE, all = TRUE)
  if (length(hits) != 1L) fail()
  end <- rawToChar(b[seq.int(hits, length(b))])
  if (!grepl("^startxref\\s+[0-9]+\\s+%%EOF\\s*$", end, perl = TRUE)) fail()
  xref <- as.integer(regmatches(end, regexec("^startxref\\s+([0-9]+)", end, perl = TRUE))[[1]][2])
  if (is.na(xref) || xref < 1L || xref >= hits) fail()
  tail <- rawToChar(b[seq.int(xref + 1L, length(b))])
  lines <- strsplit(tail, "\n", fixed = TRUE)[[1]]
  if (length(lines) < 5L || lines[1] != "xref" || !grepl("^0 [0-9]+$", lines[2])) fail()
  n <- as.integer(sub("^0 ", "", lines[2]))
  if (is.na(n) || n < 2L || length(lines) < n + 3L || lines[n + 3L] != "trailer") fail()
  entries <- lines[seq.int(3L, n + 2L)]
  if (any(!grepl("^[0-9]{10} [0-9]{5} [nf] *$", entries))) fail()
  offsets <- as.integer(substr(entries, 1, 10))
  active <- substr(entries, 18, 18) == "n"
  # Unused Quartz font slots may be marked in-use even though their offset is zero.
  empty <- active & offsets == 0L
  entries[empty] <- "0000000000 00000 f "
  active[empty] <- FALSE
  if (anyNA(offsets) || any(offsets[active] < 1L | offsets[active] >= xref) ||
      anyDuplicated(offsets[active]) || any(substr(entries[active], 12, 16) != "00000")) fail()
  trailer <- paste(lines[seq.int(n + 3L, length(lines))], collapse = "\n")
  if (grepl("/Prev|/Encrypt|/XRefStm", trailer)) fail()
  object <- function(id) {
    if (id < 1L || id >= n || !active[id + 1L]) fail()
    start <- offsets[id + 1L]
    stop <- min(c(offsets[active & offsets > start], xref))
    z <- b[seq.int(start + 1L, stop)]
    if (!identical(z[seq_len(nchar(paste0(id, " 0 obj")))], charToRaw(paste0(id, " 0 obj")))) fail()
    list(raw = z, start = start, end = stop)
  }
  page_ids <- which(vapply(seq_len(n - 1L), function(id) {
    if (!active[id + 1L]) return(FALSE)
    z <- object(id)$raw
    length(grepRaw("/Type /Page ", head(z, 500L), fixed = TRUE)) > 0L
  }, logical(1)))
  if (length(page_ids) != 1L) fail()
  page <- rawToChar(object(page_ids)$raw)
  ref <- regmatches(page, regexec("/Contents ([0-9]+) 0 R", page))[[1]]
  if (length(ref) != 2L) fail()
  id <- as.integer(ref[2]); obj <- object(id)
  stream <- grepRaw("\nstream\n", obj$raw, fixed = TRUE)
  if (length(stream) != 1L) fail()
  header <- rawToChar(obj$raw[seq_len(stream - 1L)])
  len <- regmatches(header, regexec("/Length ([0-9]+)\\s*>>", header, perl = TRUE))[[1]]
  if (length(len) != 2L || !grepl("/Filter /FlateDecode", header, fixed = TRUE)) fail()
  len <- as.integer(len[2]); first <- stream + 8L
  if (is.na(len) || len < 1L || first + len > length(obj$raw)) fail()
  content <- memDecompress(obj$raw[seq.int(first, first + len - 1L)], "gzip", asChar = TRUE)
  number <- "(-?[0-9]+(?:\\.[0-9]+)?)"
  pair <- paste(number, number)
  rect <- paste(pair, "m", pair, "l", pair, "l", pair, "l h f")
  pattern <- paste0("(?s)^(q\\s+Q\\s+q\\s+.{0,200}?)", gsub(" ", "\\s+", rect, fixed = TRUE), "(\\s+Q)")
  m <- regmatches(content, regexec(pattern, content, perl = TRUE, useBytes = TRUE))[[1]]
  if (!length(m)) fail()
  xy <- as.numeric(m[3:10]); w <- width * 72; h <- height * 72
  same <- function(a, expected) all(abs(a - expected) <= pmax(0.01, abs(expected) * 1e-5))
  normal <- c(0, h, w, h, w, 0, 0, 0)
  if (same(xy, normal)) return(invisible(FALSE))
  oversized <- c(0, h, w * 72, h, w * 72, h - h * 72, 0, h - h * 72)
  if (!same(xy, oversized)) fail()
  replacement <- paste0(m[2], sprintf("0 %.6f m %.6f %.6f l %.6f 0 l 0 0 l h f", h, w, h, w), m[11])
  rest <- charToRaw(content)[-seq_len(nchar(m[1], type = "bytes"))]
  compressed <- memCompress(c(charToRaw(replacement), rest), "gzip")
  header <- sub("/Length [0-9]+", paste("/Length", length(compressed)), header)
  rebuilt <- c(charToRaw(paste0(header, "\nstream\n")), compressed,
               obj$raw[seq.int(first + len, length(obj$raw))])
  delta <- length(rebuilt) - length(obj$raw)
  offsets[active & offsets > obj$start] <- offsets[active & offsets > obj$start] + delta
  entries[active] <- paste0(sprintf("%010d", offsets[active]), substring(entries[active], 11L))
  trailer <- sub("startxref\\s+[0-9]+", paste0("startxref\n", xref + delta), trailer, perl = TRUE)
  out <- c(b[seq_len(obj$start)], rebuilt,
           if (obj$end < xref) b[seq.int(obj$end + 1L, xref)] else raw(0),
           charToRaw(paste0("xref\n0 ", n, "\n", paste(entries, collapse = "\n"), "\n", trailer, "\n")))
  writeBin(out, file)
  invisible(TRUE)
}
