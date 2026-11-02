ANNOTATION_SPECIES <- c(Human = 9606L, Pig = 9823L, Mouse = 10090L, Rat = 10116L)
ANNOTATION_ORTHOLOGS <- readRDS("data/species_orthologs.rds")
ANNOTATION_RESOURCES <- list(
  ortholog_snapshot = jsonlite::read_json("data/species_orthologs_provenance.json", simplifyVector = TRUE),
  human_identifier_database = paste("org.Hs.eg.db", utils::packageVersion("org.Hs.eg.db")),
  human_gene_sets_sha256 = digest::digest(GENESETS, algo = "sha256"))
if (!identical(digest::digest(file = "data/species_orthologs.rds", algo = "sha256"),
               ANNOTATION_RESOURCES$ortholog_snapshot$sha256)) stop("Bundled ortholog snapshot checksum mismatch.")

annotation_organism <- function(x) {
  if (is.null(x) || !length(x)) return(NULL)
  if (length(x) != 1L || !is.atomic(x)) stop("Declare one source organism per study.")
  if (is.na(x) || !nzchar(trimws(as.character(x)))) return(NULL)
  trimws(as.character(x))
}

annotation_resource_hint <- function(organism, id_type = "symbol", id_species = "native",
                                     mode = "human", mapping_choice = "builtin") {
  if (identical(mode, "none")) return("Native feature IDs: no human mapping. Add a custom GMT for gene-set analyses.")
  if (identical(mapping_choice, "custom"))
    return("Use your mapping table and record its source/version. The built-in mapping is not used for this layer.")
  sp <- if (identical(id_species, "human")) "Human" else organism
  resource <- if (identical(sp, "Human") && id_type == "symbol") {
    "Human gene symbols are used as declared; no conversion table is needed."
  } else if (identical(sp, "Human") && id_type %in% c("ensembl", "entrez", "uniprot")) {
    paste("Human identifier lookup:", ANNOTATION_RESOURCES$human_identifier_database, "(offline).")
  } else if (!is.null(sp) && sp %in% c("Pig", "Mouse", "Rat") && id_type %in% c("symbol", "ensembl", "entrez")) {
    paste0(sp, " to human orthologs: bundled HGNC HCOP / babelgene 22.9 snapshot (offline).")
  } else {
    return("No built-in mapping for this organism/ID combination. Upload a custom mapping, or choose Native IDs only.")
  }
  paste(resource, "Built-in gene sets:", paste(COLLECTIONS, collapse = ", "), "— no GMT upload needed.")
}

annotation_read_gmt <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- lines[nzchar(trimws(lines))]
  rows <- strsplit(lines, "\t", fixed = TRUE)
  if (!length(rows) || any(lengths(rows) < 3L))
    stop("GMT requires a set name, description and at least one feature, separated by tabs.")
  nm <- vapply(rows, `[`, character(1), 1L)
  if (any(!nzchar(nm)) || anyDuplicated(toupper(gsub("_", " ", nm))))
    stop("GMT set names must be nonempty and unique (ignoring case and replacing underscores with spaces).")
  sets <- setNames(lapply(rows, function(r) unique(r[-c(1, 2)][nzchar(r[-c(1, 2)])])), nm)
  if (any(lengths(sets) == 0L)) stop("GMT contains an empty feature set.")
  list(sets = sets, lines = lines,
       hash = digest::digest(lapply(sets[order(names(sets))], sort), algo = "sha256"))
}

annotation_custom_map <- function(path, source) {
  if (is.null(source) || !nzchar(trimws(source))) stop("Provide the custom mapping source and version.")
  m <- byod_read(path)
  if (!all(c("source_id", "target_id") %in% names(m)))
    stop("Mapping CSV requires source_id and target_id columns; target_id must be a human gene symbol.")
  m$source_id <- trimws(as.character(m$source_id)); m$target_id <- trimws(as.character(m$target_id))
  if (anyNA(m[, c("source_id", "target_id")]) || any(!nzchar(m$source_id)) || any(!nzchar(m$target_id)))
    stop("Mapping source_id and target_id cannot be blank.")
  m <- unique(m)
  pairs <- unique(m[, c("source_id", "target_id")])
  f <- table(pairs$source_id); r <- table(pairs$target_id)
  m$eligible <- f[m$source_id] == 1L & r[m$target_id] == 1L
  if ("relationship" %in% names(m)) m$eligible <- m$eligible & m$relationship == "one_to_one"
  m$eligible[is.na(m$eligible)] <- FALSE
  m$relationship <- ifelse(m$eligible, "unique_in_user_mapping", "ambiguous_in_user_mapping")
  m$source <- source; m$version <- digest::digest(file = path, algo = "sha256")
  m
}

annotation_config <- function(organism, config = NULL) {
  base <- list(mode = if (!is.null(organism) && organism %in% names(ANNOTATION_SPECIES)) "human" else "none",
               id_type = "symbol", id_species = "native", source = NULL, mapping = NULL, gmt = NULL)
  if (!is.null(config)) base[names(config)] <- config
  if (!base$mode %in% c("human", "none")) stop("Unknown annotation target.")
  if (!base$id_type %in% c("symbol", "ensembl", "entrez", "uniprot", "other")) stop("Unknown identifier type.")
  if (!base$id_species %in% c("native", "human")) stop("Unknown identifier namespace.")
  base
}

annotation_build <- function(ids, organism = NULL, config = NULL) {
  organism <- annotation_organism(organism)
  cfg <- annotation_config(organism, config)
  sp <- if (cfg$id_species == "human") "Human" else organism
  q <- as.character(ids)
  if (cfg$id_type == "ensembl" && is.null(cfg$mapping)) q <- sub("\\.[0-9]+$", "", q)
  out <- data.frame(source_id = ids, query_id = q, target_id = NA_character_,
                    status = "unmapped", relationship = NA_character_, source = NA_character_,
                    version = NA_character_, stringsAsFactors = FALSE)
  candidates <- NULL
  if (cfg$mode == "human") {
    if (!is.null(cfg$mapping)) {
      candidates <- cfg$mapping
    } else if (identical(sp, "Human") && cfg$id_type == "symbol") {
      candidates <- data.frame(source_id = q, target_id = q, eligible = TRUE,
        relationship = "declared_human_symbol", source = "user-declared human symbols", version = "not converted")
    } else if (identical(sp, "Human") && cfg$id_type %in% c("ensembl", "entrez", "uniprot")) {
      keytype <- c(ensembl = "ENSEMBL", entrez = "ENTREZID", uniprot = "UNIPROT")[[cfg$id_type]]
      db <- org.Hs.eg.db::org.Hs.eg.db
      keys <- intersect(q, AnnotationDbi::keys(db, keytype = keytype))
      if (length(keys)) {
        z <- suppressMessages(AnnotationDbi::select(db, keys = keys, keytype = keytype, columns = "SYMBOL"))
        candidates <- data.frame(source_id = z[[keytype]], target_id = z$SYMBOL, eligible = TRUE,
          relationship = "human_identifier", source = "org.Hs.eg.db",
          version = as.character(utils::packageVersion("org.Hs.eg.db")))
      }
    } else if (!is.null(sp) && sp %in% c("Pig", "Mouse", "Rat") && cfg$id_type %in% c("symbol", "ensembl", "entrez")) {
      z <- ANNOTATION_ORTHOLOGS[ANNOTATION_ORTHOLOGS$organism == sp, ]
      candidates <- data.frame(source_id = as.character(z[[cfg$id_type]]), target_id = z$human_symbol,
        eligible = z$relationship == "one_to_one_in_HCOP_snapshot" & z$support_n >= 3L,
        relationship = z$relationship, source = paste0("HCOP via babelgene; ", z$support_n, " resources"), version = "22.9")
    }
  }
  if (!is.null(candidates) && nrow(candidates)) {
    candidates <- candidates[, c("source_id", "target_id", "eligible", "relationship", "source", "version")]
    candidates <- candidates[!is.na(candidates$source_id) & !is.na(candidates$target_id), ]
    pairs <- unique(candidates[, c("source_id", "target_id")])
    ambiguous <- names(which(table(pairs$source_id) > 1L))
    candidates$eligible <- candidates$eligible & !candidates$source_id %in% ambiguous
    candidates <- candidates[candidates$source_id %in% q, ]
    if (nrow(candidates)) {
      candidates <- candidates[order(!candidates$eligible), ]
      candidates <- candidates[!duplicated(candidates[, c("source_id", "target_id")]), ]
      # Keep every candidate in the report, including excluded relationships.
      out <- merge(out[, c("source_id", "query_id")], candidates,
        by.x = "query_id", by.y = "source_id", all.x = TRUE, sort = FALSE)
      out$eligible[is.na(out$eligible)] <- FALSE
      out$status <- ifelse(is.na(out$target_id), "unmapped", ifelse(out$eligible, "mapped", "ambiguous"))
      out$status[out$status == "ambiguous" & !is.na(out$relationship) &
        out$relationship == "one_to_one_in_HCOP_snapshot"] <- "insufficient_support"
      out$eligible <- NULL
      links <- unique(out[out$status == "mapped", c("source_id", "target_id")])
      duplicate_targets <- names(which(table(links$target_id) > 1L))
      out$status[out$status == "mapped" & out$target_id %in% duplicate_targets] <- "duplicate_target"
      out <- unique(out)
    }
  }
  out <- out[order(match(out$source_id, ids), out$target_id), ]
  rownames(out) <- NULL
  if (cfg$mode == "none") out$status <- "annotation_disabled"
  list(config = cfg, organism = organism, target = if (cfg$mode == "human") "Human gene symbol" else "native IDs",
       mapping = out, gmt = cfg$gmt)
}

annotation_get <- function(d) {
  if (!is.null(d$annotation)) return(d$annotation)
  # Legacy internal DS layers already carry the atlas's human-symbol annotation.
  if (length(d$species) == 1L && d$species %in% c("Pig", "Human"))
    return(annotation_build(rownames(d$expr), d$species, list(id_species = "human", mode = "human")))
  annotation_build(rownames(d$expr), d$organism)
}

annotation_keys <- function(d) {
  m <- annotation_get(d)$mapping
  m <- unique(m[m$status == "mapped" & m$source_id %in% rownames(d$expr), c("source_id", "target_id")])
  setNames(m$source_id, m$target_id)
}

annotation_collections <- function(d) {
  a <- annotation_get(d)
  c(if (a$config$mode == "human") COLLECTIONS, if (!is.null(a$gmt)) "Custom GMT")
}

annotation_sets <- function(d = NULL, collection) {
  if (is.null(d)) return(GENESETS[[collection]])
  a <- annotation_get(d)
  if (identical(collection, "Custom GMT")) return(a$gmt$sets %||% list())
  if (a$config$mode != "human" || !collection %in% COLLECTIONS) return(list())
  key <- annotation_keys(d)
  lapply(GENESETS[[collection]], function(s) unname(key[intersect(s, names(key))]))
}

annotation_universe <- function(d, collection) {
  if (is.null(d)) return(NULL)
  if (identical(collection, "Custom GMT")) return(rownames(d$expr))
  unname(annotation_keys(d))
}

annotation_compatible <- function(dA, dB, collection) {
  if (!collection %in% intersect(annotation_collections(dA), annotation_collections(dB))) return(FALSE)
  if (identical(collection, "Custom GMT")) {
    a <- annotation_get(dA); b <- annotation_get(dB)
    return(!is.null(dA$organism) && identical(dA$organism, dB$organism) && identical(a$gmt$hash, b$gmt$hash))
  }
  TRUE
}

annotation_summary <- function(d) {
  a <- annotation_get(d); m <- a$mapping
  sprintf("Source organism: %s; IDs: %s (%s). Annotation: %s; %d/%d features mapped uniquely; %d ambiguous or duplicate-target features. %s",
    a$organism %||% "not declared", a$config$id_type, a$config$id_species, a$target,
    length(unique(m$source_id[m$status == "mapped"])), nrow(d$expr),
    length(unique(m$source_id[m$status %in% c("ambiguous", "duplicate_target")])),
    if (!is.null(a$gmt)) "Custom GMT uses the original feature IDs." else "Expression feature IDs are unchanged.")
}

annotation_write_study <- function(study, directory, include_input = TRUE) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(ANNOTATION_RESOURCES, file.path(directory, "resources.json"), auto_unbox = TRUE, pretty = TRUE)
  for (i in seq_along(study$layers)) {
    d <- study$layers[[i]]; a <- annotation_get(d)
    path <- file.path(directory, sprintf("layer_%02d", i)); dir.create(path, showWarnings = FALSE)
    utils::write.csv(a$mapping, file.path(path, "mapping.csv"), row.names = FALSE, na = "")
    cfg <- a$config; cfg$mapping <- NULL; cfg$gmt <- NULL
    jsonlite::write_json(list(layer = names(study$layers)[i], organism = d$organism, settings = cfg,
      summary = annotation_summary(d), gmt_sha256 = a$gmt$hash),
      file.path(path, "annotation.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
    if (!is.null(a$gmt)) writeLines(a$gmt$lines, file.path(path, "custom.gmt"), useBytes = TRUE)
    if (!is.null(a$config$mapping)) utils::write.csv(a$config$mapping, file.path(path, "supplied_mapping.csv"), row.names = FALSE)
    if (include_input && !is.null(d$input_expression)) {
      con <- gzfile(file.path(path, "input_expression.csv.gz"), "wt")
      tryCatch(utils::write.csv(d$input_expression, con, row.names = FALSE), finally = close(con))
    }
  }
  if (include_input && !is.null(study$input_metadata))
    utils::write.csv(study$input_metadata, file.path(directory, "input_metadata.csv"), row.names = FALSE)
  invisible(directory)
}
