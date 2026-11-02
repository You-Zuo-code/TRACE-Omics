BUNDLE_SLUG <- function(x, max_n = 60L) {

  x <- suppressWarnings(as.character(x)[1])
  if (is.null(x) || is.na(x) || !nzchar(x)) return("unnamed")
  original <- enc2utf8(x)
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", gsub("_{2,}", "_", x))
  if (!nzchar(x)) return(paste0("study-", substr(digest::digest(original, algo = "sha256"), 1, 16)))
  if (nchar(x) <= max_n) return(x)
  cut <- substr(x, 1, max_n)
  at <- regexpr("[_.-][^_.-]*$", cut)
  if (at > max_n * 0.6) cut <- substr(cut, 1, at - 1)
  gsub("^_+|_+$", "", cut)
}

bundle_study_id <- function(study) {
  if (is.null(study)) return(NULL)
  study$dataset_id %||% paste0("legacy-", substr(digest::digest(
    list(study$name, study$meta, lapply(study$layers, function(d)
      list(d$expr, d$subject, d$nest, d$tp_levels, d$organism, d$annotation))), algo = "sha256"), 1, 20))
}

bundle_new_study_id <- function(session, slot) {
  n <- (session$userData$study_generation %||% 0L) + 1L
  session$userData$study_generation <- n
  paste0("upload-", BUNDLE_SLUG(slot), "-", sprintf("%06d", n))
}

BUNDLE_PAIR_DATASET <- "_cross-study"

bundle_recorded_layers <- function(inventory)
  unique(inventory$layer[inventory$dataset != BUNDLE_PAIR_DATASET])

BUNDLE_PAIR_SEP <- "--vs--"
bundle_pair_folder <- function(a, b) paste0(BUNDLE_SLUG(a), BUNDLE_PAIR_SEP, BUNDLE_SLUG(b))
bundle_pair_members <- function(folder) strsplit(folder, BUNDLE_PAIR_SEP, fixed = TRUE)[[1]]
stopifnot(identical(BUNDLE_SLUG(paste0("a", BUNDLE_PAIR_SEP, "b")),
                    paste0("a", BUNDLE_PAIR_SEP, "b")))   # the separator must survive slugging

bundle_pair_is_about <- function(folder, dataset)
  BUNDLE_SLUG(dataset) %in% bundle_pair_members(folder)

BUNDLE_FIGURE <- list(width = 8, height = 6, dpi = 300)

bundle_archive_format <- function() {
  if (requireNamespace("zip", quietly = TRUE)) return("zip")
  if (nzchar(Sys.which("zip"))) return("zip")
  "tar.gz"
}

bundle_write_archive <- function(stage, file) {
  if (!length(list.files(stage, recursive = TRUE))) return(FALSE)

  owd <- getwd()
  on.exit(setwd(owd), add = TRUE)
  run_in_stage <- function(expr) tryCatch({ setwd(stage); force(expr); setwd(owd); TRUE },
                                          error = function(e) { setwd(owd); FALSE },
                                          warning = function(w) { setwd(owd); FALSE })
  wrote <- function() file.exists(file) && file.size(file) > 0

  if (requireNamespace("zip", quietly = TRUE)) {
    ok <- tryCatch({ zip::zipr(file, file.path(stage, list.files(stage)), recurse = TRUE); TRUE },
                   error = function(e) FALSE)
    if (ok && wrote()) return(TRUE)
  }
  if (nzchar(Sys.which("zip")) && run_in_stage(utils::zip(file, ".", flags = "-r9Xq")) && wrote())
    return(TRUE)
  if (identical(bundle_archive_format(), "zip")) return(FALSE)
  run_in_stage(utils::tar(file, files = list.files("."), compression = "gzip",
                          tar = "internal")) && wrote()
}

bundle_root <- function(session) {
  ud <- session$userData
  if (is.null(ud$bundle_id))
    ud$bundle_id <- basename(tempfile("session_"))
  if (is.null(ud$bundle_root)) {
    ud$bundle_root <- file.path(tempdir(), paste0("trace_bundle_", ud$bundle_id))
    dir.create(ud$bundle_root, recursive = TRUE, showWarnings = FALSE)
  }
  ud$bundle_root
}

bundle_is_pair <- function(module)
  grepl("^cross-dataset|^Cross-study", module)

bundle_path <- function(session, dataset, module, layer) {
  fold <- analysis_folder(module)
  if (is.na(fold)) return(NA_character_)
  ds <- if (bundle_is_pair(module)) BUNDLE_PAIR_DATASET else BUNDLE_SLUG(dataset)
  p <- file.path(bundle_root(session), ds, fold, BUNDLE_SLUG(layer))
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

bundle_capture <- function(session, panel, e, gg, dataset, layer) {
  if (is.null(e) || is.null(panel)) return(invisible(FALSE))
  i <- match(panel, EXPORT_PANEL_MAP$panel)
  if (is.na(i)) return(invisible(FALSE))
  module <- EXPORT_PANEL_MAP$module[i]; figname <- EXPORT_PANEL_MAP$figure[i]
  dir <- bundle_path(session, dataset, module, layer)
  if (is.na(dir)) return(invisible(FALSE))
  status_file <- file.path(dir, paste0("capture_", BUNDLE_SLUG(panel), ".json"))
  previous <- if (file.exists(status_file)) jsonlite::read_json(status_file) else list()

  old_files <- unlist(previous$files, use.names = FALSE)
  if (length(old_files)) unlink(file.path(dir, basename(old_files)))
  written <- character(); errors <- character()
  attempt <- function(name, expr) {
    tryCatch({ force(expr); written <<- c(written, name) },
      error = function(err) errors <<- c(errors, paste(name, conditionMessage(err), sep = ": ")))
  }
  tabs <- tryCatch(xs_result_tables(e$result), error = function(err) {
    errors <<- c(errors, paste("result tables:", conditionMessage(err))); list()
  })
  for (nm in names(tabs))
    attempt(sprintf("result_%s.csv", BUNDLE_SLUG(nm)),
      utils::write.csv(tabs[[nm]], file.path(dir, sprintf("result_%s.csv", BUNDLE_SLUG(nm))), row.names = FALSE))

  pd <- tryCatch({
    s <- e$summary
    if (!is.null(e$csv) && is.data.frame(e$csv)) e$csv
    else if (is.matrix(s)) data.frame(row = rownames(s), s, check.names = FALSE)
    else if (is.data.frame(s) && nrow(s)) s else NULL
  }, error = function(err) { errors <<- c(errors, paste("plotted data:", conditionMessage(err))); NULL })
  if (!is.null(pd) && nrow(pd))
    attempt(sprintf("%s_plotted_data.csv", BUNDLE_SLUG(figname)),
      utils::write.csv(pd, file.path(dir, sprintf("%s_plotted_data.csv", BUNDLE_SLUG(figname))), row.names = FALSE))

  if (inherits(gg, "ggplot")) {
    stem <- BUNDLE_SLUG(figname)
    for (fmt in c("png", "pdf"))
      attempt(paste0(stem, ".", fmt), byod_save_figure(gg, file.path(dir, paste0(stem, ".", fmt)), fmt,
        BUNDLE_FIGURE$width, BUNDLE_FIGURE$height, BUNDLE_FIGURE$dpi))
  }
  jsonlite::write_json(list(panel = panel, dataset_id = dataset, layer = layer,
    random_seed = trace_seed(session),
    captured_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), files = written,
    status = if (length(errors)) "incomplete export" else "complete", errors = errors),
    status_file, auto_unbox = TRUE, pretty = TRUE)
  if (length(errors)) warning(paste(c("Incomplete session export:", errors), collapse = "\n"), call. = FALSE)
  invisible(!length(errors))
}

bundle_inventory <- function(session) {
  root <- bundle_root(session)
  if (!dir.exists(root)) return(data.frame())
  rows <- list()
  for (ds in list.dirs(root, recursive = FALSE)) {
    for (an in list.dirs(ds, recursive = FALSE)) {
      for (ly in list.dirs(an, recursive = FALSE)) {
        fs <- list.files(ly, full.names = TRUE)
        if (!length(fs)) next
        rows[[length(rows) + 1L]] <- data.frame(
          dataset = basename(ds), analysis = basename(an), layer = basename(ly),
          n_tables = sum(grepl("^result_", basename(fs))),
          n_figures = sum(grepl("\\.png$", basename(fs))),
          bytes = sum(file.size(fs)), stringsAsFactors = FALSE)
      }
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$analysis, out$layer), , drop = FALSE]
}

bundle_filter_log <- function(log_records, dataset) {
  if (!length(log_records)) return(list())
  dataset <- as.character(dataset)[1]
  Filter(function(r) {
    if (isTRUE(bundle_is_pair(r$module %||% ""))) {

      sides <- c(r$params$study_id_A %||% r$params$study_A,
                 r$params$study_id_B %||% r$params$study_B)
      sides <- sides[!vapply(sides, is.null, logical(1))]
      if (!length(sides)) return(FALSE)     # unattributed: cannot claim it belongs
      return(dataset %in% as.character(sides))
    }
    st <- r$params$study_id %||% r$params$study
    !is.null(st) && nzchar(as.character(st)[1]) && identical(as.character(st)[1], dataset)
  }, log_records)
}

bundle_assemble <- function(session, dataset, file, log_records = list(),
                            design_text = NULL, methods_text = NULL,
                            env_text = NULL, study_name = NULL, study = NULL, related_studies = list(),
                            include_comparison_inputs = FALSE) {

  log_records <- bundle_filter_log(log_records, dataset)
  root <- bundle_root(session)
  stage <- tempfile("trace_archive_")
  dir.create(stage, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(stage, recursive = TRUE), add = TRUE)

  srcs <- c(file.path(root, BUNDLE_SLUG(dataset)), file.path(root, BUNDLE_PAIR_DATASET))
  man <- list()
  for (src in srcs) {
    if (!dir.exists(src)) next
    is_pair_src <- identical(basename(src), BUNDLE_PAIR_DATASET)
    for (an in list.dirs(src, recursive = FALSE)) {
      for (ly in list.dirs(an, recursive = FALSE)) {

        if (is_pair_src && !bundle_pair_is_about(basename(ly), dataset)) next
        dest <- file.path(stage, basename(an), basename(ly))
        dir.create(dest, recursive = TRUE, showWarnings = FALSE)
        figs <- 0L; tabs <- character(0)
        for (f in list.files(ly, full.names = TRUE)) {
          b <- basename(f)
          if (!file.copy(f, file.path(dest, b), overwrite = TRUE)) stop("Could not copy archive output: ", b)
          if (grepl("^result_", b)) tabs <- c(tabs, b)
          if (grepl("\\.png$", b)) figs <- figs + 1L
        }

        mod <- ANALYSIS_REGISTRY$module[match(sub("^[0-9]+_", "", basename(an)),
                                              ANALYSIS_REGISTRY$folder)]

        recs <- if (is_pair_src)
          Filter(function(r) identical(bundle_pair_folder(r$params$study_id_A %||% r$params$study_A %||% "",
                                                          r$params$study_id_B %||% r$params$study_B %||% ""),
                                       basename(ly)), log_records)
          else log_records
        runs <- bundle_runs(recs, mod)
        if (nrow(runs)) utils::write.csv(runs, file.path(dest, "runs.csv"), row.names = FALSE)
        man[[length(man) + 1L]] <- data.frame(
          analysis = basename(an), layer = basename(ly),
          module = mod %||% NA_character_,
          result_tables = paste(tabs, collapse = "; "),
          figures = figs, runs_recorded = nrow(runs), stringsAsFactors = FALSE)
      }
    }
  }

  dec <- bundle_declined(log_records)

  if (!length(man) && !nrow(dec) && is.null(study)) {
    writeLines("This session produced no analyses to archive.", file)
    return(invisible(FALSE))
  }
  manifest <- if (length(man)) do.call(rbind, man) else
    data.frame(analysis = character(), layer = character(), module = character(),
               result_tables = character(), figures = integer(), runs_recorded = integer(),
               stringsAsFactors = FALSE)

  manifest$outcome <- if (nrow(manifest)) "completed" else character(0)
  if (nrow(manifest)) for (i in seq_len(nrow(manifest))) {
    cf <- list.files(file.path(stage, manifest$analysis[i], manifest$layer[i]),
      pattern = "^capture_.*\\.json$", full.names = TRUE)
    if (length(cf) && any(vapply(cf, function(f)
      !identical(jsonlite::read_json(f)$status, "complete"), logical(1))))
      manifest$outcome[i] <- "incomplete export; see capture JSON"
  }
  if (nrow(dec)) {
    dec$analysis <- vapply(dec$module, function(m) analysis_folder(m) %||% NA_character_, character(1))
    manifest <- rbind(manifest, data.frame(
      analysis = dec$analysis, layer = "-", module = dec$module,
      result_tables = "", figures = 0L, runs_recorded = 1L,
      outcome = paste("declined:", dec$reason), stringsAsFactors = FALSE))
  }
  manifest <- manifest[order(manifest$analysis, manifest$layer), , drop = FALSE]
  utils::write.csv(manifest, file.path(stage, "MANIFEST.csv"), row.names = FALSE)

  wr <- function(nm, txt) if (!is.null(txt) && length(txt))
    writeLines(as.character(txt), file.path(stage, nm))
  wr("session_log.txt", vapply(log_records, function(r) {
    kv <- if (length(r$params)) paste(sprintf("%s=%s", names(r$params),
            vapply(r$params, function(x) paste(as.character(x), collapse = ","), character(1))),
            collapse = "; ") else ""
    sprintf("[%s] %-34s %s", r$t, r$module, kv) }, character(1)))
  wr("methods.txt", strwrap(methods_text %||% "", 96))
  wr("design.txt", design_text)
  wr("environment.txt", env_text)
  wr("README.txt", bundle_readme(study_name %||% dataset, nrow(manifest)))
  if (!is.null(study)) annotation_write_study(study, file.path(stage, "annotation"))
  pair_dirs <- basename(list.dirs(file.path(root, BUNDLE_PAIR_DATASET), recursive = TRUE))
  for (other in related_studies) {
    id <- bundle_study_id(other)
    if (is.null(id) || identical(id, dataset)) next
    if (any(c(bundle_pair_folder(dataset, id), bundle_pair_folder(id, dataset)) %in% pair_dirs))
      annotation_write_study(other, file.path(stage, "comparison_annotation", BUNDLE_SLUG(id)),
                             include_input = isTRUE(include_comparison_inputs))
  }

  invisible(bundle_write_archive(stage, file))
}

bundle_runs <- function(log_records, module) {
  if (is.na(module) || !length(log_records)) return(data.frame())
  keep <- Filter(function(r) identical(r$module, module), log_records)
  if (!length(keep)) return(data.frame())
  cols <- unique(unlist(lapply(keep, function(r) names(r$params))))
  do.call(rbind, lapply(seq_along(keep), function(i) {
    r <- keep[[i]]
    v <- vapply(cols, function(k) {
      x <- r$params[[k]]
      if (is.null(x)) NA_character_ else paste(as.character(x), collapse = ",") }, character(1))
    cbind(data.frame(run = i, time = r$t, stringsAsFactors = FALSE),
          as.data.frame(as.list(v), stringsAsFactors = FALSE))
  }))
}

bundle_declined <- function(log_records) {
  if (!length(log_records)) return(data.frame())
  d <- Filter(function(r) identical(r$params$outcome, "declined"), log_records)
  if (!length(d)) return(data.frame())
  out <- do.call(rbind, lapply(d, function(r) data.frame(
    module = r$module, reason = as.character(r$params$reason %||% "no reason recorded")[1],
    stringsAsFactors = FALSE)))
  out[!duplicated(out$module), , drop = FALSE]
}

bundle_readme <- function(study, n) c(
  "TRACE-Omics session archive",
  "===========================",
  "",
  sprintf("Study: %s", study),
  sprintf("Written: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("Analyses recorded: %d", n),
  "",
  "One folder per analysis. The two-digit number belongs to the ANALYSIS and not",
  "to the order the analyses were run in, so an archive from another study has",
  "the same folder names and the two can be compared directly - `diff -r`, or two",
  "windows side by side.",
  "",
  "Inside each analysis folder, one folder per omics layer, containing:",
  "",
  "  result_*.csv              what the analysis produced - the full tables and",
  "                            matrices, not the values the figure draws",
  "  *_plotted_data.csv        the values the figure draws",
  sprintf("  *.png, *.pdf              the figures, as displayed, at %g x %g in and %d DPI.",
          BUNDLE_FIGURE$width, BUNDLE_FIGURE$height, BUNDLE_FIGURE$dpi),
  "                            For a figure at a different size, use the export",
  "                            toolbar on that panel in the application.",
  "  runs.csv                  every parameter set this analysis was run with.",
  "                            The folder holds the most recent run.",
  "",
  "At the top level:",
  "",
  "  MANIFEST.csv              every analysis, its layer, its tables and figures,",
  "                            and its outcome. An analysis the engine REFUSED",
  "                            appears here with the reason and has no folder.",
  "  design.txt                the design that was declared for this study",
  "  session_log.txt           every analysis run, with its parameters",
  "  methods.txt               a Methods paragraph assembled from those parameters",
  "  environment.txt           R and package versions for this instance",
  "",
  "Numbers are append-only across versions of the platform: a later release may",
  "add an analysis at the end, and will not renumber an existing one.")
