# Labels and exports require a UTF-8 locale.
if (!isTRUE(l10n_info()[["UTF-8"]])) {
  got <- ""
  for (lc in c("en_US.UTF-8", "C.UTF-8", "en_GB.UTF-8", "UTF-8"))
    if (nzchar(got <- suppressWarnings(Sys.setlocale("LC_CTYPE", lc)))) break
  if (!nzchar(got))
    warning("No UTF-8 locale available; dashes and symbols in labels may render as escapes.",
            call. = FALSE)
}

suppressPackageStartupMessages({
  library(stats); library(utils); library(grDevices)
  library(shiny); library(bslib); library(plotly); library(DT)
  library(limma); library(matrixStats); library(scales); library(fgsea); library(ggplot2)
  library(colourpicker); library(viridisLite); library(ggrepel)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

TRACE_DEFAULT_SEED <- 1234L
trace_validate_seed <- function(seed) {
  x <- suppressWarnings(as.numeric(seed))
  if (length(x) != 1L || !is.finite(x) || x < 0 || x > 2147483646 || x != floor(x))
    stop("Random seed must be an integer from 0 to 2147483646.")
  as.integer(x)
}
trace_seed <- function(session = shiny::getDefaultReactiveDomain()) {
  getter <- if (!is.null(session)) session$userData$analysis_seed else NULL
  trace_validate_seed(if (is.function(getter)) getter() else TRACE_DEFAULT_SEED)
}

# Restore the caller's RNG state so one session cannot change another's draws.
trace_rng_scope <- function(seed = trace_seed()) {
  seed <- trace_validate_seed(seed)
  kind <- RNGkind()
  had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  previous <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
  set.seed(seed)
  function() {
    do.call(RNGkind, as.list(kind))
    if (had) assign(".Random.seed", previous, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }
}

trace_compute_umap <- function(pca, seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  nnb <- max(2L, min(15L, nrow(pca) - 1L))
  um <- if (requireNamespace("uwot", quietly = TRUE))
    tryCatch(uwot::umap(pca, n_neighbors = nnb, min_dist = 0.3,
      n_threads = 1L, n_sgd_threads = 1L, verbose = FALSE),
      error = function(e) pca[, 1:2, drop = FALSE]) else pca[, 1:2, drop = FALSE]
  data.frame(UMAP1 = um[, 1], UMAP2 = um[, 2], row.names = rownames(pca))
}
trace_umap <- function(d, seed = trace_seed()) {
  seed <- trace_validate_seed(seed)
  if (identical(d$umap_seed, seed)) return(d$umap)
  session <- shiny::getDefaultReactiveDomain()
  key <- paste(d$species, d$omic, sep = " / ")
  cache <- if (!is.null(session)) session$userData$umap_cache else NULL
  entry <- if (!is.null(cache)) cache[[key]] else NULL
  if (!is.null(entry) && identical(entry$seed, seed) && identical(entry$pca, d$pca))
    return(entry$coordinates)
  coordinates <- trace_compute_umap(d$pca, seed)
  if (!is.null(session)) {
    cache[[key]] <- list(seed = seed, pca = d$pca, coordinates = coordinates)
    session$userData$umap_cache <- cache
  }
  coordinates
}

source("R/fonts.R", local = TRUE)
BYOD_EXPORT_FONT <- trace_font(trimws(Sys.getenv("TRACE_PLOT_FONT", "")))
BYOD_UI_FONT <- tryCatch(trace_font(cjk = TRUE), error = function(e) BYOD_EXPORT_FONT)
byod_theme_bw <- function(base_size = 11, base_family = BYOD_UI_FONT, ...)
  ggplot2::theme_bw(base_size = base_size, base_family = base_family, ...)

atlas <- readRDS("data/atlas.rds")
DS <- atlas$datasets

# The deployment can override the per-file upload limit in megabytes.
BYOD_MAX_UPLOAD_MB <- as.numeric(Sys.getenv("NMP_MAX_UPLOAD_MB", "2048"))
options(shiny.maxRequestSize = BYOD_MAX_UPLOAD_MB * 1024^2)

# Convert bundled ln(TPM+1) to log2(TPM+1); contrasts remain differences on this scale.
for (k in c("pig_mrna", "hum_mrna")) {
  DS[[k]]$expr <- DS[[k]]$expr / log(2)
  DS[[k]]$unit <- "log2(TPM+1)"; atlas$datasets[[k]]$unit <- "log2(TPM+1)"
}

NAVY <- "#1E3A5F"; TEAL <- "#2C8C99"; PIG <- "#B0413E"; HUMAN <- "#2E6CA4"
INK <- "#22303C"; MUTED <- "#6E7E8D"; GOLD <- "#E0A03A"; LIGHT <- "#F5F7FA"
TP_COL <- c(in_vivo = "#9AA7B2", preNMP = "#5B7FA6", "1h" = "#1E6E8C",
            "2h" = "#1F8A8A", "3h" = "#3AA76D", "6h" = "#E0A03A")
WI_COL <- c(WI = "#B0413E", NWI = "#2E6CA4")

TPWI_PAIRED <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
                 "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
TPWI_COL <- NULL
for (.k in names(DS)) {
  .d <- DS[[.k]]
  if (!is.null(.d$meta$WI_status)) {
    .lev <- as.vector(t(outer(.d$tp_levels, c("NWI", "WI"), paste, sep = "_")))
    .tw  <- ifelse(is.na(.d$meta$WI_status), NA,
                   paste0(as.character(.d$meta$Timepoint), "_", .d$meta$WI_status))
    DS[[.k]]$meta$Timepoint_WI <- factor(.tw, levels = .lev)
    DS[[.k]]$tpwi_levels <- .lev
    if (is.null(TPWI_COL)) TPWI_COL <- setNames(TPWI_PAIRED[seq_along(.lev)], .lev)
  }
}
rm(.k, .d, .lev, .tw)

# The bundled map records one-to-one and ambiguous human-pig orthologs.
ORTHOLOG_MAP <- tryCatch(
  { m <- utils::read.csv("data/ortholog_map.csv", stringsAsFactors = FALSE)
    stats::setNames(split(m, seq_len(nrow(m))), m$human_symbol) },
  error = function(e) NULL)
ortholog_note <- function(gene) {          # short provenance string for a queried gene
  if (is.null(ORTHOLOG_MAP) || is.null(ORTHOLOG_MAP[[gene]])) return("")
  r <- ORTHOLOG_MAP[[gene]]
  if (identical(r$class, "one_to_one")) sprintf("1:1 ortholog (%d DBs)", r$support_n)
  else "⚠ ambiguous ortholog (see Table S5)"
}

ds_key <- function(species, omic)
  paste0(tolower(substr(species, 1, if (species == "Pig") 3 else 3)),
         "_", tolower(omic))

KEYMAP <- c("Pig|mRNA" = "pig_mrna", "Pig|Protein" = "pig_protein",
            "Human|mRNA" = "hum_mrna", "Human|Protein" = "hum_protein")
get_ds <- function(species, omic) DS[[ KEYMAP[[paste(species, omic, sep = "|")]] ]]

SPECIES <- c("Pig", "Human"); OMICS <- c("mRNA", "Protein")

features_of <- function(species, omic) get_ds(species, omic)$features
has_feature <- function(d, g) !is.null(g) && g %in% rownames(d$expr)
expr_vec <- function(d, g) if (has_feature(d, g)) as.numeric(d$expr[g, ]) else NULL

grp_factor <- function(d, var) {
  if (var == "Timepoint") return(factor(d$meta$Timepoint, levels = d$tp_levels))
  if (var == "Timepoint_WI" && !is.null(d$meta$Timepoint_WI))
    return(factor(d$meta$Timepoint_WI, levels = d$tpwi_levels))
  x <- d$meta[[var]]
  if (is.null(x)) return(factor(rep(NA_character_, nrow(d$meta))))   # var absent for this species (e.g. WI_status on human)
  if (is.numeric(x)) {
    u <- sort(unique(x[!is.na(x)]))
    if (length(u) <= 2) return(factor(as.character(x), levels = as.character(u)))
    med <- stats::median(x, na.rm = TRUE)
    f <- factor(ifelse(is.na(x), NA, ifelse(x > med, "high", "low")), levels = c("low", "high"))

    if (nlevels(droplevels(f)) < 2)
      f <- factor(ifelse(is.na(x), NA, ifelse(x >= med, "high", "low")), levels = c("low", "high"))
    return(f)
  }

  if (is.factor(x)) return(droplevels(x))
  factor(as.character(x))
}
is_num_var <- function(d, var) var != "Timepoint" && is.numeric(d$meta[[var]])

group_vars <- function(d) {
  v <- "Timepoint"
  if (!is.null(d$meta$WI_status)) v <- c(v, "WI_status", "Timepoint_WI")
  if (!is.null(d$clinical_cat)) v <- c(v, d$clinical_cat)
  if (!is.null(d$clinical_num)) v <- c(v, d$clinical_num)
  v
}

col_groups <- function(d) {
  v <- "Timepoint"
  if (!is.null(d$meta$WI_status)) v <- c(v, "WI_status", "Timepoint_WI")
  if (!is.null(d$clinical_cat)) v <- c(v, d$clinical_cat)
  v
}
group_levels <- function(d, var) levels(droplevels(grp_factor(d, var)))

byod_nest_ok <- function(meta, subject_col, nest_col) {
  if (is.null(subject_col) || is.null(nest_col)) return(FALSE)
  if (!all(c(subject_col, nest_col) %in% names(meta))) return(FALSE)
  s <- as.character(meta[[subject_col]]); b <- as.character(meta[[nest_col]])
  if (anyNA(s) || anyNA(b)) return(FALSE)
  ns <- length(unique(s)); nb <- length(unique(b))
  if (nb < 2 || nb >= ns) return(FALSE)                       # must be strictly coarser
  all(tapply(b, s, function(z) length(unique(z))) == 1L)
}

byod_nest_candidates <- function(meta, subject_col, exclude = character()) {
  if (is.null(subject_col)) return(character())
  cand <- setdiff(names(meta), c(subject_col, "Timepoint", "Timepoint_WI", exclude))
  cand <- cand[!vapply(meta[cand], is.numeric, logical(1))]
  cand[vapply(cand, function(cc) byod_nest_ok(meta, subject_col, cc), logical(1))]
}

byod_nested_paired_fit <- function(E, subj, grp, nest, min_units = MIN_UNITS$group_contrast) {
  subj <- droplevels(factor(subj)); grp <- droplevels(factor(grp)); nest <- factor(nest)
  if (nlevels(grp) != 2L) return(NULL)
  lv <- levels(grp)
  keep <- vapply(levels(subj), function(s) {
    g <- as.character(grp[subj == s]); all(lv %in% g)
  }, logical(1))
  use <- levels(subj)[keep]
  if (length(use) < 2L) return(NULL)

  D <- vapply(use, function(s) {
    i1 <- which(subj == s & grp == lv[1]); i2 <- which(subj == s & grp == lv[2])
    rowMeans(E[, i2, drop = FALSE], na.rm = TRUE) - rowMeans(E[, i1, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(E)))
  dimnames(D) <- list(rownames(E), use)

  nest_of <- vapply(use, function(s) as.character(nest[subj == s][1]), character(1))
  units <- unique(nest_of)
  if (length(units) < min_units)
    return(list(declined = sprintf(
      "This contrast rests on %d independent unit%s (%s), below the floor of %d. It has %d paired subject%s, but they are not independent of one another.",
      length(units), if (length(units) == 1L) "" else "s", paste(head(units, 4), collapse = ", "),
      min_units, length(use), if (length(use) == 1L) "" else "s")))
  U <- vapply(units, function(u) rowMeans(D[, nest_of == u, drop = FALSE], na.rm = TRUE), numeric(nrow(D)))
  dimnames(U) <- list(rownames(D), units)

  des <- matrix(1, ncol(U), 1, dimnames = list(colnames(U), "(Intercept)"))
  fit <- tryCatch(eBayes(lmFit(U, des)), error = function(e) e)
  if (inherits(fit, "error"))
    return(list(declined = sprintf("The nested one-sample model could not be fitted across %d %s: %s",
                                   length(units), "units", conditionMessage(fit))))
  tt <- topTable(fit, coef = 1, number = Inf, sort.by = "P")

  list(tab = tt, n_units = length(units), n_subjects = length(use),
       n_profiles = sum(subj %in% use), excluded_subjects = setdiff(levels(subj), use))
}

run_limma <- function(d, groupvar, g1, g2, paired = TRUE, at_tp = "All", mat = NULL) {
  E0 <- if (is.null(mat)) d$expr else mat
  gf <- grp_factor(d, groupvar)
  sel <- which(as.character(gf) %in% c(g1, g2) & !is.na(gf))
  if (!is.null(at_tp) && at_tp != "All") sel <- sel[d$meta$Timepoint[sel] == at_tp]
  sel <- sel[colnames(d$expr)[sel] %in% colnames(E0)]
  mt <- d$meta[sel, , drop = FALSE]; E <- E0[, colnames(d$expr)[sel], drop = FALSE]
  grp <- factor(ifelse(as.character(gf[sel]) == g2, "g2", "g1"), levels = c("g1", "g2"))
  if (nlevels(droplevels(grp)) < 2 || min(table(grp)) < 2)
    return(data.frame(Feature = character(), logFC = numeric(), AveExpr = numeric(),
                      t = numeric(), P.Value = numeric(), adj.P.Val = numeric()))

  subj <- if (groupvar == "Timepoint") mt$Kidney_ID
          else if (groupvar == "WI_status") mt$Pig_Number else NULL
  use_pair <- paired && !is.null(subj) && !any(is.na(subj))

  nest_col <- d$nest %||% NULL
  if (use_pair && groupvar == "Timepoint" && !is.null(nest_col) &&
      byod_nest_ok(mt, "Kidney_ID", nest_col)) {
    nf <- byod_nested_paired_fit(E, as.character(subj), grp, as.character(mt[[nest_col]]))
    if (!is.null(nf) && is.null(nf$declined)) {
      tt <- nf$tab
      return(data.frame(Feature = rownames(tt),
                        logFC = round(tt$logFC, 3), AveExpr = round(tt$AveExpr, 3),
                        t = round(tt$t, 2), P.Value = signif(tt$P.Value, 3),
                        adj.P.Val = signif(tt$adj.P.Val, 3),
                        row.names = NULL, stringsAsFactors = FALSE))
    }
  }
  if (use_pair) {
    subj <- factor(subj)
    tab <- table(subj, grp); keepS <- rownames(tab)[tab[, "g1"] > 0 & tab[, "g2"] > 0]
    keep <- subj %in% keepS
    E <- E[, keep, drop = FALSE]; grp <- droplevels(grp[keep]); subj <- droplevels(subj[keep])
    use_pair <- nlevels(subj) >= 2 && nlevels(grp) == 2
  }
  design <- if (use_pair) model.matrix(~ subj + grp) else model.matrix(~ grp)
  fit <- eBayes(lmFit(E, design))
  coefn <- grep("^grp", colnames(design), value = TRUE)[1]
  tt <- topTable(fit, coef = coefn, number = Inf, sort.by = "P")
  data.frame(Feature = rownames(tt),
             logFC = round(tt$logFC, 3), AveExpr = round(tt$AveExpr, 3),
             t = round(tt$t, 2), P.Value = signif(tt$P.Value, 3),
             adj.P.Val = signif(tt$adj.P.Val, 3),
             row.names = NULL, stringsAsFactors = FALSE)
}

trend_p <- function(d, g) {
  v <- expr_vec(d, g); if (is.null(v)) return(NA_real_)
  tp <- d$meta$Timepoint
  if (length(unique(tp)) < 2) return(NA_real_)
  summary(aov(v ~ tp))[[1]][["Pr(>F)"]][1]
}

module_score <- function(d, genes) {
  genes <- intersect(genes, rownames(d$expr))
  if (length(genes) < 1) return(NULL)
  z <- t(scale(t(d$expr[genes, , drop = FALSE])))
  z[!is.finite(z)] <- NA
  list(score = colMeans(z, na.rm = TRUE), n = length(genes), genes = genes)
}

ply <- function(p, ...) layout(p, font = list(family = "Helvetica, Arial, sans-serif",
        color = INK), paper_bgcolor = "white", plot_bgcolor = "white",
        margin = list(l = 60, r = 20, t = 40, b = 50), ...) %>%
  config(displaylogo = FALSE,
         modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"))

summ_stats <- function(v, grp) {
  s <- tapply(v, grp, function(x) c(m = mean(x, na.rm = TRUE),
        se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))), n = sum(!is.na(x))))
  data.frame(group = names(s), mean = sapply(s, `[`, 1),
             se = sapply(s, `[`, 2), n = sapply(s, `[`, 3), row.names = NULL)
}

plot_timecourse <- function(species, gene, palette = PAL_NAMES[1], cols = NULL, splitby = "WI_status") {

  is_split <- !(splitby %in% c("Timepoint", "Timepoint_WI"))
  mk <- function(omic) {
    d <- get_ds(species, omic); v <- expr_vec(d, gene)
    if (is.null(v)) return(NULL)
    mt <- d$meta; use_split <- is_split && !is.null(mt[[splitby]])
    if (use_split) {
      g <- interaction(factor(mt$Timepoint, levels = d$tp_levels), grp_factor(d, splitby), sep = "@", drop = TRUE)
      st <- summ_stats(v, g)
      parts <- do.call(rbind, strsplit(as.character(st$group), "@", fixed = TRUE))
      st$Timepoint <- factor(parts[, 1], levels = d$tp_levels); st$Group <- parts[, 2]
    } else {
      st <- summ_stats(v, factor(mt$Timepoint, levels = d$tp_levels))
      st$Timepoint <- factor(st$group, levels = d$tp_levels); st$Group <- "all"
    }
    st$omic <- omic; st[order(st$Timepoint), c("omic", "Timepoint", "Group", "mean", "se", "n")]
  }
  rna <- mk("mRNA"); pro <- mk("Protein")
  if (is.null(rna) && is.null(pro)) return(NULL)
  summary <- do.call(rbind, list(rna, pro))
  summary$omic <- factor(summary$omic, levels = c("mRNA", "Protein")); rownames(summary) <- NULL
  splits <- unique(as.character(summary$Group))
  if (all(c("NWI", "WI") %in% splits)) splits <- c("NWI", "WI")
  splitby_eff <- if (identical(splits, "all")) NULL else splitby
  if (is.null(cols)) cols <- resolve_appearance(splits, palette, splitby = splitby_eff)
  else { base <- resolve_appearance(splits, palette, splitby = splitby_eff)
    keep <- names(cols)[names(cols) %in% splits]; base[keep] <- cols[keep]; cols <- base[splits] }

  one <- function(st, ttl, unit, showleg) {
    if (is.null(st)) return(plot_ly(x = numeric(), y = numeric(), type = "scatter", mode = "markers") %>% ply() %>%
      layout(annotations = list(text = paste0(ttl, ": not detected"), showarrow = FALSE,
             x = 0.5, y = 0.5, xref = "paper", yref = "paper", font = list(color = MUTED))))
    tpl <- levels(st$Timepoint); p <- plot_ly()
    grps <- intersect(splits, unique(st$Group))
    for (w in grps) {
      sub <- st[st$Group == w, , drop = FALSE]; col <- cols[[w]]; nm <- if (w == "all") ttl else w
      p <- add_trace(p, data = sub, x = ~as.numeric(Timepoint), y = ~mean,
        type = "scatter", mode = "lines+markers", name = nm, legendgroup = nm, showlegend = showleg,
        line = list(color = col, width = 3), marker = list(color = col, size = 8),
        error_y = list(array = sub$se, color = col, thickness = 1.3, width = 4),
        hovertemplate = "%{text}<br>mean=%{y:.2f}<extra></extra>",
        text = ~paste0(Timepoint, if (w != "all") paste0(" (", w, ")") else ""))
    }
    p %>% ply(title = list(text = ttl, x = 0.02, font = list(size = 14)),
      xaxis = list(tickvals = seq_along(tpl), ticktext = tpl, title = ""),
      yaxis = list(title = unit), showlegend = (showleg && length(grps) > 1),
      legend = list(orientation = "h", y = 1.12, x = 0))
  }
  pr <- one(rna, "mRNA", "log2 expression", showleg = !is.null(rna))
  pp <- one(pro, "Protein", "log2 intensity", showleg = is.null(rna))
  plt <- subplot(pr, pp, nrows = 1, margin = 0.06, titleY = TRUE) %>%
    layout(title = list(text = paste0("<b>", gene, "</b> — ", species,
           " kidney across perfusion"), x = 0.02, font = list(size = 16)))
  list(plot = plt, summary = summary, splits = splits, ylab = "expression",
       ttl = paste0(gene, " — ", species, " kidney across perfusion"), cols = cols)
}

plot_violin <- function(species, omic, gene, groupvar, type = "violin", cols = NULL) {
  d <- get_ds(species, omic); v <- expr_vec(d, gene)
  ylab <- paste0(omic, " (", d$unit, ")"); ttl <- paste0(gene, " — ", species, " ", omic)
  fail <- function(msg) list(plot = emptyplot(msg),
    summary = data.frame(Group = character(), mean = numeric(), sd = numeric(),
                         se = numeric(), median = numeric(), n = integer()),
    splits = character(0), ylab = ylab, ttl = msg, cols = setNames(character(0), character(0)))
  if (is.null(v)) return(fail(paste0(gene, " not detected in ", species, " ", omic)))
  g <- droplevels(grp_factor(d, groupvar)); keep <- !is.na(g) & !is.na(v)
  df <- data.frame(v = v[keep], g = g[keep]); if (!nrow(df)) return(fail("no data"))
  lv <- levels(df$g); agg <- function(f) tapply(df$v, df$g, f)
  summary <- data.frame(Group = factor(lv, levels = lv),
    mean = as.numeric(agg(function(x) mean(x, na.rm = TRUE))[lv]),
    sd   = as.numeric(agg(function(x) sd(x, na.rm = TRUE))[lv]),
    median = as.numeric(agg(function(x) median(x, na.rm = TRUE))[lv]),
    n = as.integer(agg(function(x) sum(!is.na(x)))[lv]), row.names = NULL, stringsAsFactors = FALSE)
  summary$se <- summary$sd / sqrt(pmax(summary$n, 1))
  attr(summary, "raw") <- data.frame(Group = factor(as.character(df$g), levels = lv), value = df$v)
  base <- if (groupvar == "Timepoint") TP_COL[lv] else if (groupvar == "WI_status") WI_COL[lv] else
          if (groupvar == "Timepoint_WI") TPWI_COL[lv] else resolve_appearance(lv, PAL_NAMES[1], splitby = NULL)
  base <- setNames(as.character(base), lv); if (!is.null(cols)) { ok <- intersect(names(cols), lv); base[ok] <- cols[ok] }
  cols <- base
  p <- plot_ly(df, x = ~g, y = ~v, color = ~g, colors = unname(cols[lv]),
    type = type, box = list(visible = TRUE), meanline = list(visible = TRUE),
    points = "all", jitter = 0.3, pointpos = 0, marker = list(size = 4, opacity = .5)) %>%
    ply(xaxis = list(title = groupvar), yaxis = list(title = ylab), showlegend = FALSE,
        title = list(text = paste0("<b>", gene, "</b> — ", species, " ", omic), x = 0.02, font = list(size = 15)))
  list(plot = p, summary = summary, splits = lv, ylab = ylab, ttl = ttl, cols = cols)
}

dr_kind_of <- function(d, colorby) {
  if (identical(colorby, "Gene expression")) return("gradient")
  if (!is.null(colorby) && is_num_var(d, colorby)) return("gradient")
  "discrete"
}
plot_dimred <- function(species, omic, type = "UMAP", colorby = "Timepoint", gene = NULL,
                        cols = NULL, scale_stops = NULL) {
  d  <- get_ds(species, omic)
  co <- if (type == "UMAP") trace_umap(d) else as.data.frame(d$pca[, 1:2])
  ax <- if (type == "UMAP") c("UMAP1", "UMAP2") else
        c(paste0("PC1 (", d$varexp[1], "%)"), paste0("PC2 (", d$varexp[2], "%)"))
  df <- data.frame(sample = rownames(co), x = co[, 1], y = co[, 2],
                   d$meta[rownames(co), , drop = FALSE], stringsAsFactors = FALSE)
  ttl <- paste0(species, " ", omic, " — ", type); kind <- dr_kind_of(d, colorby)
  if (kind == "gradient") {
    if (identical(colorby, "Gene expression")) {
      if (!has_feature(d, gene)) return(list(plot = emptyplot(paste(gene, "not detected")),
        summary = data.frame(), kind = "gradient", levels = character(0), ttl = ttl, ax = ax, scale_title = gene %||% "value"))
      cv <- as.numeric(d$expr[gene, rownames(co)]); sttl <- gene
    } else { cv <- df[[colorby]]; sttl <- colorby }
    df$cval <- cv; stops <- scale_stops %||% heat_stops("YlOrRd")
    p <- plot_ly(df, x = ~x, y = ~y, type = "scatter", mode = "markers",
        marker = list(color = ~cval, colorscale = heat_colorscale(stops), showscale = TRUE,
          colorbar = list(title = sttl, len = .6), size = 8, line = list(width = .5, color = "white")),
        text = ~paste0(sample, "<br>", Timepoint, "<br>", sttl, " = ", signif(cval, 3)), hoverinfo = "text") %>%
      ply(xaxis = list(title = ax[1]), yaxis = list(title = ax[2]),
          title = list(text = ttl, x = .02, font = list(size = 15)))
    return(list(plot = p, summary = data.frame(sample = df$sample, x = df$x, y = df$y, value = df$cval,
        stringsAsFactors = FALSE), kind = "gradient", levels = character(0), ttl = ttl, ax = ax, scale_title = sttl))
  }
  cb <- colorby
  g  <- droplevels(grp_factor(d, cb)[match(rownames(co), rownames(d$meta))])
  df$grp <- g; keep <- !is.na(df$grp); df <- df[keep, , drop = FALSE]; g <- droplevels(df$grp)
  levels <- levels(g)
  if (is.null(cols)) {
    cols <- if (cb == "Timepoint") TP_COL[levels] else if (cb == "WI_status") WI_COL[levels] else
            if (cb == "Timepoint_WI") TPWI_COL[levels] else resolve_appearance(levels, PAL_NAMES[1], splitby = cb)
  } else cols <- cols[levels]
  cols <- setNames(as.character(cols), levels)
  p <- plot_ly(df, x = ~x, y = ~y, color = g, colors = unname(cols[levels]),
      type = "scatter", mode = "markers", marker = list(size = 8, line = list(width = .5, color = "white")),
      text = ~paste0(sample, "<br>", Timepoint), hoverinfo = "text") %>%
    ply(xaxis = list(title = ax[1]), yaxis = list(title = ax[2]),
        title = list(text = ttl, x = .02, font = list(size = 15)))
  list(plot = p, summary = data.frame(sample = df$sample, x = df$x, y = df$y,
       group = as.character(g), stringsAsFactors = FALSE),
       kind = "discrete", levels = levels, ttl = ttl, ax = ax, scale_title = cb)
}

plot_volcano <- function(tt, lab, fdr = 0.05, lfc = 1, cols = NULL, n_label = 0) {
  tt$sig <- factor(ifelse(tt$adj.P.Val < fdr & abs(tt$logFC) >= lfc,
            ifelse(tt$logFC > 0, "Up", "Down"), "n.s."), levels = VOLCANO_LEVELS)
  tt$y <- -log10(pmax(tt$P.Value, 1e-300))
  cc <- VOLCANO_DEFAULT; if (!is.null(cols)) { ok <- intersect(names(cols), names(cc)); cc[ok] <- cols[ok] }
  lab_df <- if (n_label > 0) { ss <- tt[tt$sig %in% c("Up","Down"), , drop = FALSE]
              head(ss[order(ss$adj.P.Val, -abs(ss$logFC)), ], n_label) } else tt[0, , drop = FALSE]
  ymax <- max(tt$y, na.rm = TRUE)
  p <- plot_ly(tt, x = ~logFC, y = ~y, color = ~sig, colors = cc[VOLCANO_LEVELS],
    type = "scatter", mode = "markers", marker = list(size = 5, opacity = .7),
    text = ~paste0(Feature, "<br>logFC=", round(logFC, 3), "<br>FDR=", signif(adj.P.Val, 3)), hoverinfo = "text") %>%
    ply(xaxis = list(title = "Mean difference (analysis scale)"), yaxis = list(title = "-log10 P"),
        title = list(text = lab, x = .02, font = list(size = 14)), legend = list(orientation = "h", y = 1.1),
        shapes = list(list(type = "line", x0 = lfc, x1 = lfc, y0 = 0, y1 = ymax, line = list(dash = "dot", color = MUTED, width = 1)),
                      list(type = "line", x0 = -lfc, x1 = -lfc, y0 = 0, y1 = ymax, line = list(dash = "dot", color = MUTED, width = 1))))
  if (nrow(lab_df)) p <- add_annotations(p, data = lab_df, x = ~logFC, y = ~y, text = ~Feature,
      showarrow = TRUE, arrowhead = 0, ax = 12, ay = -12, font = list(size = 10, color = INK),
      bgcolor = "rgba(255,255,255,.6)")
  summary <- tt[, c("Feature","logFC","AveExpr","t","P.Value","adj.P.Val","sig","y")]
  attr(summary, "thr") <- c(fdr = fdr, lfc = lfc); attr(summary, "lab") <- as.character(lab_df$Feature)
  list(plot = p, summary = summary, splits = VOLCANO_LEVELS, ylab = "-log10 P", ttl = lab, cols = cc)
}

plot_mp <- function(species, groupvar, g1, g2, paired = TRUE, lfc = 1,
                    cols = NULL, palette = PAL_NAMES[1], topn = 200) {
  QUAD <- c("both up", "both down", "discordant"); ylab <- "Protein mean difference"
  rt <- run_limma(get_ds(species, "mRNA"), groupvar, g1, g2, paired)
  pt <- run_limma(get_ds(species, "Protein"), groupvar, g1, g2, paired)
  fail <- function(msg, m = data.frame()) list(plot = emptyplot(msg), summary = m, splits = character(0),
    r = NA_real_, p = NA_real_, n = nrow(m), ylab = ylab, ttl = msg, cols = setNames(character(0), character(0)), tab = m)
  if (!nrow(rt) || !nrow(pt)) return(fail("too few samples in one omic"))
  m <- merge(rt[, c("Feature","logFC","adj.P.Val")], pt[, c("Feature","logFC","adj.P.Val")],
             by = "Feature", suffixes = c("_mRNA", "_protein"))
  if (nrow(m) < 3) return(fail("fewer than 3 shared features", m))
  ct <- suppressWarnings(cor.test(m$logFC_mRNA, m$logFC_protein))
  m$quad <- ifelse(m$logFC_mRNA > 0 & m$logFC_protein > 0, "both up",
            ifelse(m$logFC_mRNA < 0 & m$logFC_protein < 0, "both down", "discordant"))
  m <- head(m[order(-abs(m$logFC_mRNA + m$logFC_protein)), ], topn); m$quad <- factor(m$quad, levels = QUAD)
  splits <- intersect(QUAD, as.character(unique(m$quad))); cols <- fb_cols(splits, cols, palette)
  ttl <- paste0("mRNA vs protein change (", g2, " vs ", g1, ")  ·  r = ", round(ct$estimate, 2),
                ", p = ", signif(ct$p.value, 2), ", n = ", nrow(m))
  p <- plot_ly(m, x = ~logFC_mRNA, y = ~logFC_protein, color = ~factor(quad, splits), colors = unname(cols[splits]),
      type = "scatter", mode = "markers", marker = list(size = 6, opacity = .75),
      text = ~paste0(Feature, "<br>mRNA=", logFC_mRNA, " protein=", logFC_protein), hoverinfo = "text") %>%
    ply(xaxis = list(title = "mRNA mean difference", zeroline = TRUE), yaxis = list(title = ylab, zeroline = TRUE),
        legend = list(orientation = "h", y = 1.1), title = list(text = ttl, x = .02, font = list(size = 13)))
  list(plot = p, summary = m, splits = splits, r = unname(ct$estimate), p = ct$p.value, n = nrow(m),
       ylab = ylab, ttl = ttl, cols = cols, tab = m)
}

plot_heatmap <- function(species, omic, genes, by = "Timepoint", colorscale = NULL) {
  d <- get_ds(species, omic); genes <- intersect(genes, rownames(d$expr))
  ttl <- paste0(species, " ", omic, " — ", length(genes), " genes (row z-score)")
  if (length(genes) < 2) return(list(plot = emptyplot("Enter ≥2 valid genes"), found = genes, mat = NULL, ttl = ttl, by = by))
  g <- d$meta[[by]]; if (by == "Timepoint") g <- factor(g, levels = d$tp_levels)
  M <- sapply(levels(factor(g)), function(lv) rowMeans(d$expr[genes, g == lv, drop = FALSE], na.rm = TRUE))
  Z <- t(scale(t(M))); Z[!is.finite(Z)] <- 0; Z <- Z[hclust(dist(Z))$order, , drop = FALSE]
  cs <- colorscale %||% list(c(0, HUMAN), c(.5, "#FFFFFF"), c(1, PIG))
  p <- plot_ly(x = colnames(Z), y = rownames(Z), z = Z, type = "heatmap", colorscale = cs, zmid = 0,
      hovertemplate = "%{y} @ %{x}<br>z=%{z:.2f}<extra></extra>") %>%
    ply(xaxis = list(title = by), yaxis = list(title = ""), title = list(text = ttl, x = .02, font = list(size = 14)))
  list(plot = p, found = genes, mat = Z, ttl = ttl, by = by)
}

plot_modulescore <- function(species, omic, genes, splitby = "None", palette = PAL_NAMES[1]) {
  d <- get_ds(species, omic); ms <- module_score(d, genes)
  if (is.null(ms)) return(list(plot = emptyplot("No matching genes"), n = 0, p = NA,
    summary = data.frame(Timepoint = factor(levels = d$tp_levels), Group = character(),
                         mean = numeric(), se = numeric(), n = numeric()),
    splits = character(0), ylab = "module score (mean z)", ttl = "No matching genes",
    cols = setNames(character(0), character(0))))
  val <- setNames(ms$score, colnames(d$expr))
  td <- tc_data(d, val, splitby); cols <- resolve_cols(td$splits, palette, splitby)
  pval <- tryCatch(summary(aov(ms$score ~ d$meta$Timepoint))[[1]][["Pr(>F)"]][1], error = function(e) NA)
  ttl <- paste0(ms$n, " genes · ", species, " ", omic, " · ANOVA p = ", signif(pval, 2))
  ylab <- "module score (mean z)"
  list(plot = tc_plotly(td, ylab, ttl, cols), n = ms$n, p = pval, genes = ms$genes,
       summary = td$summary, cols = cols, ylab = ylab, ttl = ttl)
}

plot_crossspecies <- function(omic, gene, cols = NULL, palette = PAL_NAMES[1],
                              display = c("raw", "relative"), baseline = NULL,
                              measure = "difference", align = c("label", "position")) {
  display <- match.arg(display); align <- match.arg(align)
  allt <- atlas$tp_pig
  on <- ortholog_note(gene)          # versioned 1:1 orthology provenance for this gene
  ylab <- paste0(omic, " expression")
  ttl <- paste0("<b>", gene, "</b> — pig vs human (", omic, ")", if (nzchar(on)) paste0(" · ", on) else "")
  fail <- function(msg, yl = ylab) list(plot = emptyplot(msg),
    summary = data.frame(Timepoint = factor(levels = allt), Group = character(),
                         mean = numeric(), se = numeric(), n = numeric()),
    splits = character(0), ylab = yl, ttl = msg, cols = setNames(character(0), character(0)),
    display = display, ok = FALSE)

  dP <- get_ds("Pig", omic); dH <- get_ds("Human", omic)
  if (!has_feature(dP, gene) && !has_feature(dH, gene)) return(fail(paste(gene, "not detected")))

  if (identical(display, "raw")) {
    mk <- function(sp) { d <- get_ds(sp, omic); v <- expr_vec(d, gene); if (is.null(v)) return(NULL)
      st <- summ_stats(v, factor(d$meta$Timepoint, levels = d$tp_levels)); st <- st[!is.na(st$group), , drop = FALSE]
      data.frame(Timepoint = factor(as.character(st$group), levels = allt), Group = sp,
                 mean = st$mean, se = st$se, n = st$n, stringsAsFactors = FALSE) }
    summary <- rbind(mk("Pig"), mk("Human"))
    note <- NULL
  } else {

    al <- cs_align_timepoints(dP, dH, align)
    if (al$n < 2 && identical(align, "label")) { al <- cs_align_timepoints(dP, dH, "position"); align <- "position" }
    if (al$n < 2) return(fail(al$note %||% "The two cohorts have too few comparable timepoints."))
    i <- match(baseline, al$label); if (is.na(i)) i <- match(baseline, al$A); if (is.na(i)) i <- 1L
    if (i >= al$n) return(fail("Choose a baseline with at least one later comparable timepoint."))
    lab <- if (measure %in% names(CS_MEASURE_LAB)) CS_MEASURE_LAB[[measure]] else measure
    ylab <- paste0(lab, "  (vs ", al$label[i], ")")

    one <- function(d, tps, sp) {
      rows <- lapply(seq_along(tps), function(j) {
        if (j == i) {                       # the baseline is zero against itself, by definition
          n0 <- length(unique(as.character(d$meta[[d$subject %||% "Kidney_ID"]])))
          return(data.frame(Timepoint = tps[j], Group = sp, mean = 0, se = 0,
                            n = if (is.finite(n0) && n0 > 0) n0 else sum(d$meta$Timepoint == tps[j]),
                            stringsAsFactors = FALSE))
        }
        dd <- cs_subject_delta(d, tps[i], tps[j], measure, d$subject)
        if (is.null(dd) || !(gene %in% rownames(dd$delta))) return(NULL)
        x <- as.numeric(dd$delta[gene, ]); x <- x[is.finite(x)]
        if (!length(x)) return(NULL)
        data.frame(Timepoint = tps[j], Group = sp, mean = mean(x),
                   se = if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else 0,
                   n = length(x), stringsAsFactors = FALSE)
      })
      out <- do.call(rbind, rows)
      if (is.null(out)) return(NULL)
      out$pairing <- cs_subject_delta(d, tps[i], tps[min(i + 1L, length(tps))], measure,
                                      d$subject)$pairing %||% "paired"
      out
    }
    pg <- if (has_feature(dP, gene)) one(dP, al$A, "Pig") else NULL
    hu <- if (has_feature(dH, gene)) one(dH, al$B, "Human") else NULL
    if (is.null(pg) && is.null(hu)) return(fail(paste(gene, "has no usable baseline in either cohort")))

    relab <- function(z, own) { if (is.null(z)) return(z); z$Timepoint <- al$label[match(z$Timepoint, own)]; z }
    summary <- rbind(relab(pg, al$A), relab(hu, al$B))
    allt <- al$label
    unp <- unique(c(if (!is.null(pg)) pg$pairing[1], if (!is.null(hu)) hu$pairing[1]))
    note <- if (any(unp != "paired"))
      "At least one cohort has no subject spanning the baseline, so its values are referenced to the baseline group mean and carry between-subject variation."
    else NULL
    ttl <- paste0(ttl, " · relative to ", al$label[i],
                  if (identical(align, "position")) " (timepoints matched by position)" else "")
    summary$pairing <- NULL
  }

  if (is.null(summary) || !nrow(summary)) return(fail(paste(gene, "not detected"), ylab))
  summary$Timepoint <- factor(as.character(summary$Timepoint), levels = allt)
  splits <- intersect(c("Pig", "Human"), unique(summary$Group))
  summary$Group <- factor(summary$Group, levels = splits)
  summary <- summary[order(summary$Group, summary$Timepoint), ]
  if (is.null(cols)) cols <- resolve_appearance(splits, palette, splitby = NULL) else cols <- cols[splits]
  p <- tc_plotly(list(summary = summary, splits = splits), ylab, ttl, cols)

  if (identical(display, "relative"))
    p <- plotly::layout(p, shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
      yref = "y", y0 = 0, y1 = 0, line = list(color = MUTED, width = 1, dash = "dot"))))
  list(plot = p, summary = summary, splits = splits, ylab = ylab, ttl = ttl, cols = cols,
       display = display, note = note, ok = TRUE)
}

GSVA_MATS  <- atlas$gsva
PATHWAYS   <- sort(unique(unlist(lapply(GSVA_MATS, rownames))))
GENESETS   <- atlas$genesets
COLLECTIONS<- names(GENESETS)
source("R/annotation.R", local = TRUE)
gsva_of <- function(species, omic) GSVA_MATS[[ KEYMAP[[paste(species, omic, sep = "|")]] ]]

emptyplot <- function(msg) {
  blank <- list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
                showline = FALSE, ticks = "", title = "")
  plot_ly(x = numeric(), y = numeric(), type = "scatter", mode = "markers") %>% ply() %>%
    layout(xaxis = blank, yaxis = blank,
           annotations = list(text = paste(strwrap(msg, width = 60), collapse = "<br>"),
                              showarrow = FALSE, xref = "paper", yref = "paper",
                              x = .5, y = .5, align = "center",
                              font = list(color = MUTED)))
}

emptyplot_gg <- function(msg) {
  ggplot() +
    annotate("text", x = 0, y = 0, label = paste(strwrap(msg, width = 54), collapse = "\n"),
             size = 3.6, colour = "grey35", lineheight = 1.2) +
    theme_void(base_size = 11) +
    theme(plot.margin = margin(24, 24, 24, 24),
          panel.border = element_rect(colour = "grey85", fill = NA, linewidth = .4))
}

plot_gsva_heatmap <- function(species, omic, by = "Timepoint", colorscale = NULL) {
  d <- get_ds(species, omic); G <- gsva_of(species, omic)
  ttl <- paste0(species, " ", omic, " — Hallmark activity (GSVA)")
  if (is.null(G)) return(list(plot = emptyplot("GSVA not available"), mat = NULL, ttl = ttl, by = by))
  gf <- grp_factor(d, by); names(gf) <- rownames(d$meta); g <- droplevels(gf[colnames(G)])
  keep <- !is.na(g); G <- G[, keep, drop = FALSE]; g <- droplevels(g[keep])
  M <- sapply(levels(g), function(lv) rowMeans(G[, g == lv, drop = FALSE], na.rm = TRUE))
  M <- M[hclust(dist(M))$order, , drop = FALSE]
  cs <- colorscale %||% list(c(0, HUMAN), c(.5, "#FFFFFF"), c(1, PIG))
  p <- plot_ly(x = colnames(M), y = rownames(M), z = M, type = "heatmap", colorscale = cs, zmid = 0,
    hovertemplate = "%{y}<br>%{x}: %{z:.2f}<extra></extra>") %>%
    ply(xaxis = list(title = by), yaxis = list(title = "", tickfont = list(size = 9)), margin = list(l = 230),
        title = list(text = ttl, x = .02, font = list(size = 14)))
  list(plot = p, mat = M, ttl = ttl, by = by)
}

line_palette <- function(n) rep(c(NAVY, PIG, TEAL, GOLD, HUMAN, "#6a3d9a", "#33a02c",
  "#e31a1c", "#ff7f00", "#b15928"), length.out = n)
PAL_NAMES <- c("Brand (WI = blue/red)", "Okabe-Ito (colourblind-safe)", "Viridis",
               "Set 2", "Dark 3", "Warm", "Cold", "Spectral")
resolve_cols <- function(levels, palette, splitby = NULL) {
  n <- length(levels)
  if (grepl("^Brand", palette) && !is.null(splitby) && splitby == "WI_status")
    return(WI_COL[levels])
  if (grepl("^Brand", palette) && length(levels) && all(levels %in% c("Pig", "Human")))
    return(setNames(c(Pig = PIG, Human = HUMAN)[levels], levels))
  cols <- if (grepl("^Brand", palette)) line_palette(n)
    else if (grepl("^Okabe", palette)) unname(grDevices::palette.colors(n + 1, "Okabe-Ito", recycle = TRUE)[-1])
    else grDevices::hcl.colors(n, palette)
  setNames(cols, levels)
}

resolve_appearance <- function(levels, palette, splitby = NULL, manual = NULL, order = NULL) {
  levels <- levels[!is.na(levels) & nzchar(levels)]
  if (!length(levels)) return(setNames(character(0), character(0)))
  if (!is.null(order)) {
    order  <- order[order %in% levels]
    levels <- c(order, setdiff(levels, order))
  }
  cols <- resolve_cols(levels, palette, splitby)
  if (!is.null(manual) && length(manual)) {
    ok <- names(manual) %in% levels & !is.na(manual) &
          grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", manual)
    if (any(ok)) cols[names(manual)[ok]] <- manual[ok]
  }
  cols
}

EXPORT_FORMATS <- c("PNG (raster)" = "png", "TIFF (raster, journals)" = "tiff",
                    "JPEG (raster)" = "jpeg", "BMP (raster)" = "bmp",
                    "PDF (vector)" = "pdf", "SVG (vector, editable)" = "svg",
                    "EPS (vector)" = "eps", "PS (vector)" = "ps")
EXPORT_VECTOR  <- c("pdf", "svg", "eps", "ps")

byod_host_memory_bytes <- function() {
  num <- function(x) { x <- suppressWarnings(as.numeric(x)); if (length(x) && is.finite(x)) x else NA_real_ }

  for (f in c("/sys/fs/cgroup/memory.max",
              "/sys/fs/cgroup/memory/memory.limit_in_bytes")) {
    if (file.exists(f)) {
      v <- tryCatch(readLines(f, n = 1, warn = FALSE), error = function(e) NA)
      if (length(v) && !identical(trimws(v), "max")) {
        b <- num(v)

        if (!is.na(b) && b > 0 && b < 2^60) return(b)
      }
    }
  }
  if (file.exists("/proc/meminfo")) {
    l <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character(0))
    a <- grep("^MemAvailable:", l, value = TRUE)
    if (length(a)) { b <- num(gsub("[^0-9]", "", a[1])) * 1024; if (!is.na(b)) return(b) }
  }
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    b <- num(tryCatch(system("sysctl -n hw.memsize", intern = TRUE, ignore.stderr = TRUE),
                      error = function(e) NA))
    if (!is.na(b)) return(b)
  }
  NA_real_
}

byod_export_budget <- function() {
  cached <- getOption("byod.export.max_px")
  if (!is.null(cached) && is.finite(cached)) return(cached)
  b <- byod_host_memory_bytes()
  px <- if (is.na(b)) 60e6 else max(20e6, min(600e6, (b * 0.25) / (4 * 3)))
  options(byod.export.max_px = px)
  px
}

EXPORT_MAX_PX  <- byod_export_budget()
EXPORT_LIM <- list(w = c(1, 100), h = c(1, 100), dpi = c(36, 1200))

byod_plot_has_cjk <- function(g) {
  has <- function(x) {
    if (is.factor(x)) x <- levels(x)
    if (is.character(x)) return(any(grepl("[\\p{Han}\\p{Hiragana}\\p{Katakana}\\p{Hangul}]", x, perl = TRUE), na.rm = TRUE))
    if (is.list(x)) return(any(vapply(x, has, logical(1))))
    FALSE
  }
  has(list(g$labels, g$data, lapply(g$layers, function(x) x$data),
           lapply(g$scales$scales, function(x) list(x$name, x$labels, x$limits))))
}

byod_export_font <- function(g, family = BYOD_EXPORT_FONT) {
  if (!inherits(g, "ggplot")) return(g)
  family <- trace_font(family, cjk = byod_plot_has_cjk(g))
  g <- g + theme(text = element_text(family = family))
  # Text geoms do not all inherit the theme's font; clone before changing defaults.
  g$layers <- lapply(g$layers, function(layer) {
    if ("family" %in% names(layer$geom$default_aes) &&
        is.null(layer$mapping$family) && is.null(g$mapping$family) &&
        (is.null(layer$aes_params$family) || identical(layer$aes_params$family, ""))) {
      updated <- ggplot2::ggproto(NULL, layer)
      updated$aes_params <- utils::modifyList(layer$aes_params, list(family = family))
      return(updated)
    }
    layer
  })
  g
}

source("R/pdf_export.R", local = TRUE)

byod_pdf_device <- function() {
  cv <- suppressWarnings(tryCatch(grDevices::grSoftVersion()[["cairo"]],
                                 error = function(e) NA_character_))
  if (length(cv) == 1L && !is.na(cv) && nzchar(cv)) return(grDevices::cairo_pdf)
  if (isTRUE(capabilities("aqua"))) return(structure(function(filename, width, height, ...)
    grDevices::quartz(type = "pdf", file = filename, width = width, height = height, ...),
    trace_pdf_backend = "quartz"))
  grDevices::pdf
}

byod_save_figure <- function(g, file, fmt = "png", w = 8, h = 5, dpi = 300) {
  num <- function(x, d) { x <- suppressWarnings(as.numeric(x))
                          if (length(x) != 1L || !is.finite(x)) d else x }
  clamp <- function(x, lo, hi) max(lo, min(x, hi))
  fmt <- if (fmt %in% names(EXPORT_FORMATS) || fmt %in% EXPORT_FORMATS)
           (if (fmt %in% EXPORT_FORMATS) fmt else unname(EXPORT_FORMATS[fmt])) else "png"
  w <- clamp(num(w, 8), EXPORT_LIM$w[1], EXPORT_LIM$w[2])
  h <- clamp(num(h, 5), EXPORT_LIM$h[1], EXPORT_LIM$h[2])
  d <- clamp(num(dpi, 300), EXPORT_LIM$dpi[1], EXPORT_LIM$dpi[2])

  capped <- FALSE; requested_dpi <- d
  budget <- byod_export_budget()
  if (!fmt %in% EXPORT_VECTOR && w * h * d^2 > budget) {
    d <- max(EXPORT_LIM$dpi[1], floor(sqrt(budget / (w * h))))
    capped <- TRUE
  }

  attempts <- 0L; err <- NULL
  g <- byod_export_font(g)
  device <- if (identical(fmt, "pdf")) byod_pdf_device() else fmt
  repeat {
    attempts <- attempts + 1L
    ok <- tryCatch({ withCallingHandlers(ggsave(file, g, device = device, width = w, height = h,
                            dpi = d, units = "in", limitsize = FALSE), warning = function(wrn) {
      if (grepl("conversion failure|mbcsToSbcs|failed to load|unable to start device", conditionMessage(wrn)))
        stop(conditionMessage(wrn))
    })
    if (identical(fmt, "pdf") && identical(attr(device, "trace_pdf_backend"), "quartz"))
      byod_bound_quartz_background(file, w, h)
    TRUE },
                   error = function(e) { err <<- conditionMessage(e); FALSE })
    if (ok) break
    if (fmt %in% EXPORT_VECTOR || d <= EXPORT_LIM$dpi[1] || attempts >= 6L) break
    d <- max(EXPORT_LIM$dpi[1], floor(d / 2)); capped <- TRUE

    # Retry at a lower resolution for this export only.
  }
  if (!ok) {
    if (file.exists(file)) unlink(file)
    stop("Figure export failed in the requested ", fmt, " format: ", err, call. = FALSE)
  }
  invisible(list(format = fmt, width = w, height = h, dpi = d,
                 font = if (inherits(g, "ggplot")) g$theme$text$family else NA_character_,
                 requested_dpi = requested_dpi, capped = capped,
                 attempts = attempts, error = err, ok = ok,
                 budget_mp = budget / 1e6,
                 megapixels = if (fmt %in% EXPORT_VECTOR) NA_real_ else w * h * d^2 / 1e6))
}

apprUI <- function(id) {
  ns <- NS(id)
  tagList(

    conditionalPanel(sprintf("output['%s'] === true", ns("has_levels")),
                     selectInput(ns("palette"), "Colour palette", PAL_NAMES)),
    tags$details(style = "margin:.25rem 0 .5rem;",
      tags$summary(style = "cursor:pointer;font-weight:600;color:#1E3A5F;",
                   "Advanced: custom colours & order"),
      div(style = "padding-top:.4rem;",
        uiOutput(ns("order_ui")), uiOutput(ns("pickers")),
        helpText("Pickers start at the palette colour; change one to override it. ",
                 "Switching palette resets them."))),
    fluidRow(column(6, numericInput(ns("w"), "Width (in)", 8, EXPORT_LIM$w[1], EXPORT_LIM$w[2], .5)),
             column(6, numericInput(ns("h"), "Height (in)", 5, EXPORT_LIM$h[1], EXPORT_LIM$h[2], .5))),
    fluidRow(
      column(6, numericInput(ns("dpi"), "Resolution (DPI)", 300,
                             EXPORT_LIM$dpi[1], EXPORT_LIM$dpi[2], 25)),
      column(6, selectInput(ns("fmt"), "Format", EXPORT_FORMATS, selected = "png"))),
    uiOutput(ns("size_note")),
    div(class = "d-flex gap-1 flex-wrap",
        downloadButton(ns("dl_csv"), "Plotted data", class = "btn-sm"),
        conditionalPanel(sprintf("output['%s'] === true", ns("has_result")),
                         downloadButton(ns("dl_result"), "Analysis result", class = "btn-sm")),
        downloadButton(ns("dl_fig"), "Figure", class = "btn-sm btn-outline-primary")))
}

xs_result_tables <- function(x) {
  as_df <- function(z, nm) {
    if (is.null(z)) return(NULL)
    if (is.matrix(z) || is.array(z)) {
      d <- as.data.frame(z, stringsAsFactors = FALSE, check.names = FALSE)
      rn <- rownames(z)
      if (!is.null(rn)) d <- cbind(stats::setNames(data.frame(rn, stringsAsFactors = FALSE),
                                                   nm %||% "row"), d)
      return(d)
    }
    if (is.data.frame(z)) return(z)
    if (is.atomic(z) && length(z))
      return(data.frame(name = names(z) %||% seq_along(z), value = as.vector(z),
                        stringsAsFactors = FALSE))
    NULL
  }
  if (is.null(x)) return(list())
  if (is.data.frame(x) || is.matrix(x) || (is.atomic(x) && !is.list(x)))
    return(Filter(Negate(is.null), list(result = as_df(x, "row"))))
  out <- lapply(names(x), function(nm) as_df(x[[nm]], "row"))
  names(out) <- names(x)
  Filter(function(z) !is.null(z) && nrow(z) > 0, out)
}

xs_write_result <- function(tabs, file, prefix = "result") {
  if (!length(tabs)) { utils::write.csv(data.frame(), file, row.names = FALSE); return(invisible()) }
  if (length(tabs) == 1L) { utils::write.csv(tabs[[1]], file, row.names = FALSE); return(invisible()) }

  if (bundle_archive_format() == "zip") {
    dir <- tempfile("res_")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    fs <- vapply(names(tabs), function(nm) {
      f <- file.path(dir, sprintf("%s_%s.csv", prefix, gsub("[^A-Za-z0-9]+", "_", nm)))
      utils::write.csv(tabs[[nm]], f, row.names = FALSE); f }, character(1))
    ok <- if (requireNamespace("zip", quietly = TRUE))
            tryCatch({ zip::zipr(file, fs); TRUE }, error = function(e) FALSE)
          else tryCatch({ utils::zip(file, fs, flags = "-j9Xq"); TRUE },
                        error = function(e) FALSE, warning = function(w) FALSE)
    unlink(dir, recursive = TRUE)
    if (ok && file.exists(file) && file.size(file) > 0) return(invisible())
  }

  cols <- unique(unlist(lapply(tabs, names)))
  stacked <- do.call(rbind, lapply(names(tabs), function(nm) {
    d <- tabs[[nm]]
    for (m in setdiff(cols, names(d))) d[[m]] <- NA
    cbind(table = nm, d[, cols, drop = FALSE])
  }))
  utils::write.csv(stacked, file, row.names = FALSE)
  invisible()
}

apprServer <- function(id, levels_r, export_r) {
  moduleServer(id, function(input, output, session) {

    output$size_note <- renderUI({
      num <- function(x, d) { x <- suppressWarnings(as.numeric(x))
                              if (length(x) != 1L || !is.finite(x)) d else x }
      w <- num(input$w, 8); h <- num(input$h, 5); d <- num(input$dpi, 300)
      fmt <- input$fmt %||% "png"
      if (fmt %in% EXPORT_VECTOR)
        return(helpText(sprintf("Vector format: resolution-independent, so DPI is ignored. Size %g\u2013%g in.",
                                EXPORT_LIM$w[1], EXPORT_LIM$w[2])))
      mp <- w * h * d^2 / 1e6
      budget <- byod_export_budget()
      over <- mp > budget / 1e6
      helpText(class = if (over) "text-warning" else NULL,
        sprintf("%.0f megapixels at %g \u00d7 %g in, %g DPI. %s Ranges: %g\u2013%g in, %g\u2013%g DPI.",
                mp, w, h, d,
                if (over) sprintf("More than this machine's memory allows (about %.0f megapixels here), so it will render at about %d DPI. A larger machine would render it in full; a vector format has no such limit.",
                                  budget / 1e6, as.integer(floor(sqrt(budget / (w * h)))))
                else "Within what this machine can render.",
                EXPORT_LIM$w[1], EXPORT_LIM$w[2], EXPORT_LIM$dpi[1], EXPORT_LIM$dpi[2]))
    })

    ns <- session$ns; MAXPICK <- 12L
    base_cols <- reactiveVal(NULL)
    recompute_base <- function() {
      lv <- levels_r()
      if (!length(lv)) { base_cols(setNames(character(0), character(0))); return(invisible()) }
      ord <- isolate(input$order); ord <- if (!is.null(ord) && setequal(ord, lv)) ord else lv
      sb  <- isolate(export_r()$splitby) %||% NULL
      base_cols(resolve_appearance(lv, isolate(input$palette) %||% PAL_NAMES[1], splitby = sb, order = ord))
    }
    observeEvent(levels_r(), recompute_base(), ignoreNULL = FALSE)
    observeEvent(input$palette, recompute_base(), ignoreInit = TRUE)
    observeEvent(input$order, recompute_base(), ignoreInit = TRUE)   # keep swatches in sync after a drag

    output$order_ui <- renderUI({
      lv <- levels_r(); req(length(lv) >= 1); if (length(lv) < 2) return(NULL)
      cur <- isolate(input$order); sel <- if (!is.null(cur) && setequal(cur, lv)) cur else lv
      selectizeInput(ns("order"), "Legend / colour order (drag)", choices = lv,
                     selected = sel, multiple = TRUE, options = list(plugins = list("drag_drop")))
    })
    output$pickers <- renderUI({
      bc <- base_cols(); req(bc); if (length(bc) < 2) return(NULL)
      if (length(bc) > MAXPICK)
        return(helpText(sprintf("%d categories - manual colours hidden; use the palette.", length(bc))))
      ord <- isolate(input$order); lv <- if (!is.null(ord) && setequal(ord, names(bc))) ord else names(bc)

      tagList(lapply(lv, function(L)
        colourInput(ns(paste0("col_", match(L, names(bc)))), L, value = unname(bc[[L]]), showColour = "both")))
    })
    manual_r <- reactive({
      bc <- base_cols(); if (is.null(bc) || length(bc) < 2 || length(bc) > MAXPICK) return(NULL)
      nm <- names(bc)
      vals <- vapply(seq_along(nm), function(i) {
        cur <- input[[paste0("col_", i)]]
        if (is.null(cur) || !nzchar(cur)) return(NA_character_)
        if (tolower(trimws(cur)) == tolower(trimws(bc[[i]]))) return(NA_character_)
        cur }, character(1))
      names(vals) <- nm; vals <- vals[!is.na(vals)]; if (!length(vals)) NULL else vals
    })
    output$has_levels <- reactive(length(levels_r()) > 0)
    outputOptions(output, "has_levels", suspendWhenHidden = FALSE)

    cols_r <- reactive({
      lv <- levels_r()
      if (!length(lv)) return(setNames(character(0), character(0)))
      sb <- export_r()$splitby %||% NULL
      resolve_appearance(lv, input$palette %||% PAL_NAMES[1], splitby = sb, manual = manual_r(), order = input$order)
    })
    pref <- reactive(export_r()$prefix %||% "plot")
    result_r <- reactive(tryCatch(xs_result_tables(export_r()$result), error = function(e) list()))
    output$has_result <- reactive(length(result_r()) > 0)
    outputOptions(output, "has_result", suspendWhenHidden = FALSE)
    output$dl_result <- downloadHandler(
      filename = function() {
        n <- length(isolate(result_r()))
        sprintf("%s_result.%s", pref(),
                if (n > 1 && bundle_archive_format() == "zip") "zip" else "csv")
      },
      content = function(f) xs_write_result(result_r(), f, pref()))
    .num <- function(x, d) { x <- suppressWarnings(as.numeric(x)); if (length(x) != 1 || !is.finite(x)) d else x }
    output$dl_csv <- downloadHandler(function() paste0(pref(), "_data.csv"), function(f) {
      e <- export_r()

      if (!is.null(e$csv) && is.data.frame(e$csv)) { write.csv(e$csv, f, row.names = FALSE); return(invisible()) }
      s <- e$summary
      out <- data.frame(Timepoint = as.character(s$Timepoint), Group = as.character(s$Group),
                        mean = s$mean, se = s$se, n = s$n)
      if (!is.null(s$omic)) out <- cbind(omic = as.character(s$omic), out)
      write.csv(out, f, row.names = FALSE) })
    save_one <- function(f, fmt) {
      e <- export_r(); cols <- cols_r()
      g <- if (is.function(e$gg)) e$gg(e$summary, e$ylab, e$ttl, cols) else tc_gg(e$summary, e$ylab, e$ttl, cols)
      byod_save_figure(g, f, fmt, input$w, input$h, input$dpi)
    }
    output$dl_fig <- downloadHandler(
      filename = function() paste0(pref(), ".", input$fmt %||% "png"),
      content  = function(f) save_one(f, input$fmt %||% "png"))

    .cap <- shiny::debounce(reactive({
      e <- tryCatch(export_r(), error = function(err) NULL)
      if (is.null(e)) return(NULL)
      list(e = e, cols = tryCatch(cols_r(), error = function(err) NULL),
           context = tryCatch(session$userData$bundle_ctx(id), error = function(err) NULL))
    }), 1500)
    observeEvent(.cap(), {
      z <- .cap(); ctxf <- session$userData$bundle_ctx
      if (is.null(z) || is.null(ctxf)) return(invisible())
      cc <- z$context
      if (is.null(cc)) return(invisible())
      gg <- tryCatch(if (is.function(z$e$gg)) z$e$gg(z$e$summary, z$e$ylab, z$e$ttl, z$cols)
                     else tc_gg(z$e$summary, z$e$ylab, z$e$ttl, z$cols),
                     error = function(err) NULL)
      bundle_capture(session, id, z$e, gg, cc$dataset, cc$layer)
    }, ignoreNULL = TRUE)
    cols_r
  })
}

resolved_unit_col <- function(d, meta = NULL, keep = NULL) {
  mt <- meta %||% d$meta
  keep <- unique(c("Timepoint", keep[!is.na(keep) & nzchar(keep)]))
  keep <- intersect(keep, names(mt))
  subj <- d$subject %||% NULL
  nest <- d$nest %||% NULL
  split_varies_within_nest <- FALSE
  if (!is.null(nest) && length(keep) > 1) {
    for (k in setdiff(keep, "Timepoint")) {
      if (k %in% names(mt) &&
          any(tapply(as.character(mt[[k]]), as.character(mt[[nest]]),
                     function(x) length(unique(x)) > 1))) split_varies_within_nest <- TRUE
    }
  }
  unit_col <- if (!is.null(nest) && !split_varies_within_nest) nest else subj
  if (is.null(unit_col) || !unit_col %in% names(mt)) NULL else unit_col
}

collapse_to_resolved_unit <- function(d, val_named, keep = NULL) {
  samp <- intersect(names(val_named), rownames(d$meta))
  if (!length(samp)) return(list(value = val_named[0], meta = d$meta[0, , drop = FALSE],
                                 unit = "profile", n_units = 0L))
  v  <- val_named[samp]
  mt <- d$meta[samp, , drop = FALSE]

  keep <- unique(c("Timepoint", keep[!is.na(keep) & nzchar(keep)]))
  keep <- intersect(keep, names(mt))
  unit_col <- resolved_unit_col(d, mt, keep)
  if (is.null(unit_col))
    return(list(value = v, meta = mt, unit = "profile", n_units = length(v)))

  key <- interaction(lapply(c(unit_col, keep), function(k) as.character(mt[[k]])),
                     drop = TRUE, sep = "\r")
  agg <- tapply(v, key, mean, na.rm = TRUE)
  parts <- do.call(rbind, strsplit(names(agg), "\r", fixed = TRUE))
  nm <- make.unique(names(agg))
  meta2 <- as.data.frame(parts, stringsAsFactors = FALSE)
  names(meta2) <- c(unit_col, keep)
  rownames(meta2) <- nm
  meta2$Timepoint <- factor(meta2$Timepoint, levels = d$tp_levels)
  list(value = setNames(as.numeric(agg), nm), meta = meta2,
       unit = unit_col, n_units = length(unique(meta2[[unit_col]])))
}

byod_time_test <- function(cu, d, split_col = NULL) {
  no <- function(why) list(ok = FALSE, p = NA_real_, text = paste("no timepoint test:", why))
  if (is.null(cu$unit) || identical(cu$unit, "profile"))
    return(no("the design declares no subject or nesting level"))
  mt <- cu$meta
  if (!cu$unit %in% names(mt)) return(no("the resolved unit is not in the metadata"))
  df <- data.frame(y = as.numeric(cu$value),
                   tp = factor(as.character(mt$Timepoint), levels = d$tp_levels),
                   unit = factor(as.character(mt[[cu$unit]])), stringsAsFactors = FALSE)

  has_split <- !is.null(split_col) && !identical(split_col, "None") &&
               split_col %in% names(mt) && length(unique(as.character(mt[[split_col]]))) >= 2L
  if (has_split) df$grp <- droplevels(factor(as.character(mt[[split_col]])))
  df <- df[is.finite(df$y) & !is.na(df$tp) & !is.na(df$unit), , drop = FALSE]
  df$tp <- droplevels(df$tp); df$unit <- droplevels(df$unit)
  if (nlevels(df$tp) < 2L) return(no("only one timepoint is present"))
  if (nlevels(df$unit) < 3L)
    return(no(sprintf("only %d independent %s", nlevels(df$unit), cu$unit)))
  term <- if (has_split) "tp:grp" else "tp"
  fml  <- if (has_split) y ~ tp * grp + Error(unit) else y ~ tp + Error(unit)
  fit <- tryCatch(stats::aov(fml, data = df), error = function(e) NULL)
  if (is.null(fit)) return(no("the within-unit model is not estimable on this design"))
  pv <- NA_real_
  for (st in summary(fit)) {
    tb <- tryCatch(st[[1]], error = function(e) NULL)
    if (is.null(tb) || is.null(rownames(tb))) next
    hit <- trimws(rownames(tb)) == term
    if (any(hit) && "Pr(>F)" %in% colnames(tb)) pv <- tb[hit, "Pr(>F)"][1]
  }
  if (!is.finite(pv))
    return(no(sprintf("the within-unit %s is not estimable on this design",
                      if (has_split) "interaction" else "time effect")))
  list(ok = TRUE, p = pv, unit = cu$unit, n_units = nlevels(df$unit),
       term = term,
       text = sprintf("repeated-measures ANOVA, %s p = %s across %d %s",
                      if (has_split) sprintf("%s x time", split_col) else "time",
                      signif(pv, 2), nlevels(df$unit), cu$unit))
}

tc_data <- function(d, val_named, splitby) {
  use_split <- !is.null(splitby) && splitby != "None" && !is.null(d$meta[[splitby]])
  cu <- collapse_to_resolved_unit(d, val_named, keep = if (use_split) splitby else NULL)
  val_named <- cu$value; mt <- cu$meta
  tp <- factor(mt$Timepoint, levels = d$tp_levels)
  sf <- if (use_split) droplevels(factor(mt[[splitby]])) else NULL
  grp <- if (use_split) interaction(tp, sf, sep = "@@", drop = TRUE) else droplevels(tp)
  st <- summ_stats(val_named, grp)

  if (!nrow(st)) return(list(
    summary = data.frame(Timepoint = factor(character(0), levels = d$tp_levels),
                         Group = character(0), mean = numeric(0), se = numeric(0), n = integer(0)),
    splits = character(0)))
  if (use_split) { pr <- do.call(rbind, strsplit(as.character(st$group), "@@", fixed = TRUE))
    st$Timepoint <- factor(pr[, 1], levels = d$tp_levels); st$Group <- pr[, 2]
  } else { st$Timepoint <- factor(st$group, levels = d$tp_levels); st$Group <- "all" }
  st <- st[order(st$Timepoint), c("Timepoint", "Group", "mean", "se", "n")]
  list(summary = st, splits = unique(as.character(st$Group)),
       unit = cu$unit, n_units = cu$n_units)
}
tc_plotly <- function(td, ylab, ttl, cols) {
  st <- td$summary
  if (is.null(st) || !nrow(st)) return(emptyplot(gsub("<[^>]+>", "", if (length(ttl)) ttl else "no data")))
  tpl <- levels(st$Timepoint); p <- plot_ly()
  for (s in td$splits) { sub <- st[st$Group == s, ]; col <- if (s %in% names(cols)) cols[[s]] else MUTED
    p <- add_trace(p, data = sub, x = ~as.numeric(Timepoint), y = ~mean, type = "scatter",
      mode = "lines+markers", name = if (s == "all") "all" else s,
      line = list(color = col, width = 3), marker = list(color = col, size = 8),
      error_y = list(array = sub$se, color = col, width = 4)) }

  p %>% ply(xaxis = list(tickvals = seq_along(tpl), ticktext = tpl, title = ""),
    yaxis = list(title = ylab), showlegend = (length(td$splits) > 1),
    margin = list(l = 60, r = 20, t = 58, b = 50),
    legend = list(orientation = "h", y = 1.10, yanchor = "bottom",
                  x = 1, xanchor = "right", font = list(size = 11)),
    title = list(text = ttl, x = 0, xanchor = "left", y = 1, yanchor = "top",
                 yref = "container", pad = list(t = 8, l = 8), font = list(size = 13)))
}
tc_gg <- function(summary, ylab, ttl, cols) {
  splits <- unique(as.character(summary$Group))
  ggplot(summary, aes(Timepoint, mean, colour = Group, group = Group)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .15, linewidth = .5) +
    geom_line(linewidth = 1) + geom_point(size = 2.6) +
    scale_colour_manual(values = cols) +
    labs(x = NULL, y = ylab, title = gsub("·", "-", gsub("<[^>]+>", "", ttl)), colour = NULL) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = if (length(splits) > 1) "top" else "none",
          plot.title = element_text(size = 12), panel.grid.minor = element_blank())
}
save_gg <- function(file, summary, ylab, ttl, cols, fmt, w, h, dpi) {
  byod_save_figure(tc_gg(summary, ylab, ttl, cols), file, fmt, w, h, dpi)
}

tc_gg_facet <- function(summary, ylab, ttl, cols) {
  summary <- summary[!is.na(summary$mean), , drop = FALSE]
  if (!is.factor(summary$omic)) summary$omic <- factor(summary$omic, levels = c("mRNA", "Protein"))
  summary$omic <- droplevels(summary$omic)
  splits <- unique(as.character(summary$Group))
  ggplot(summary, aes(Timepoint, mean, colour = Group, group = Group)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = .15, linewidth = .5) +
    geom_line(linewidth = 1) + geom_point(size = 2.6) +
    scale_colour_manual(values = cols, drop = FALSE) +
    facet_wrap(~ omic, scales = "free_y") +
    labs(x = NULL, y = "mRNA = log2 · Protein = log2",
         title = gsub("·", "-", gsub("<[^>]+>", "", ttl)), colour = NULL) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = if (length(splits) > 1) "top" else "none",
          plot.title = element_text(size = 12),
          strip.background = element_rect(fill = "#EAF1F4", colour = NA),
          strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
}

HEAT_SCALE_NAMES <- c("Blue-White-Red (brand)", "RdBu", "RdYlBu", "Spectral",
                      "YlOrRd", "Viridis", "Magma", "Cividis")
HEAT_DIVERGING   <- c("Blue-White-Red (brand)", "RdBu", "RdYlBu", "Spectral")
heat_stops <- function(scale_name, reverse = FALSE, n = 11L) {
  st <- switch(scale_name,
    "Blue-White-Red (brand)" = grDevices::colorRampPalette(c(HUMAN, "#FFFFFF", PIG))(n),
    "RdBu"     = grDevices::hcl.colors(n, "RdBu"),
    "RdYlBu"   = grDevices::hcl.colors(n, "RdYlBu"),
    "Spectral" = grDevices::hcl.colors(n, "Spectral"),
    "YlOrRd"   = grDevices::hcl.colors(n, "YlOrRd"),
    "Viridis"  = viridisLite::viridis(n),
    "Magma"    = viridisLite::magma(n),
    "Cividis"  = viridisLite::cividis(n),
    grDevices::colorRampPalette(c(HUMAN, "#FFFFFF", PIG))(n))
  st <- st[!is.na(st)]                       # hcl.colors(11, RdBu/RdYlBu/Spectral) has NA at centre
  if (!length(st)) st <- grDevices::colorRampPalette(c(HUMAN, "#FFFFFF", PIG))(n)
  if (isTRUE(reverse)) st <- rev(st)
  st
}
is_diverging_scale <- function(scale_name) scale_name %in% HEAT_DIVERGING
heat_colorscale <- function(stops) {
  pos <- seq(0, 1, length.out = length(stops))
  lapply(seq_along(stops), function(i) list(pos[i], stops[i]))
}

apprUIheat <- function(id, selected = HEAT_SCALE_NAMES[1]) {
  ns <- NS(id)
  tagList(
    selectInput(ns("scale"), "Colour scale", HEAT_SCALE_NAMES, selected = selected),
    checkboxInput(ns("reverse"), "Reverse scale", FALSE),
    fluidRow(column(6, numericInput(ns("w"), "Width (in)", 8, EXPORT_LIM$w[1], EXPORT_LIM$w[2], .5)),
             column(6, numericInput(ns("h"), "Height (in)", 6, EXPORT_LIM$h[1], EXPORT_LIM$h[2], .5))),
    fluidRow(
      column(6, numericInput(ns("dpi"), "Resolution (DPI)", 300,
                             EXPORT_LIM$dpi[1], EXPORT_LIM$dpi[2], 25)),
      column(6, selectInput(ns("fmt"), "Format", EXPORT_FORMATS, selected = "png"))),
    uiOutput(ns("size_note")),

    div(class = "d-flex gap-1 flex-wrap",
        downloadButton(ns("dl_csv"), "Plotted data", class = "btn-sm"),
        conditionalPanel(sprintf("output['%s'] === true", ns("has_result")),
                         downloadButton(ns("dl_result"), "Analysis result", class = "btn-sm")),
        downloadButton(ns("dl_fig"), "Figure", class = "btn-sm btn-outline-primary")))
}
apprServerHeat <- function(id, export_r) {
  moduleServer(id, function(input, output, session) {

    output$size_note <- renderUI({
      num <- function(x, d) { x <- suppressWarnings(as.numeric(x))
                              if (length(x) != 1L || !is.finite(x)) d else x }
      w <- num(input$w, 8); h <- num(input$h, 5); d <- num(input$dpi, 300)
      fmt <- input$fmt %||% "png"
      if (fmt %in% EXPORT_VECTOR)
        return(helpText(sprintf("Vector format: resolution-independent, so DPI is ignored. Size %g\u2013%g in.",
                                EXPORT_LIM$w[1], EXPORT_LIM$w[2])))
      mp <- w * h * d^2 / 1e6
      budget <- byod_export_budget()
      over <- mp > budget / 1e6
      helpText(class = if (over) "text-warning" else NULL,
        sprintf("%.0f megapixels at %g \u00d7 %g in, %g DPI. %s Ranges: %g\u2013%g in, %g\u2013%g DPI.",
                mp, w, h, d,
                if (over) sprintf("More than this machine's memory allows (about %.0f megapixels here), so it will render at about %d DPI. A larger machine would render it in full; a vector format has no such limit.",
                                  budget / 1e6, as.integer(floor(sqrt(budget / (w * h)))))
                else "Within what this machine can render.",
                EXPORT_LIM$w[1], EXPORT_LIM$w[2], EXPORT_LIM$dpi[1], EXPORT_LIM$dpi[2]))
    })

    .num <- function(x, d) { x <- suppressWarnings(as.numeric(x)); if (length(x) != 1 || !is.finite(x)) d else x }
    pref <- reactive(export_r()$prefix %||% "heatmap")
    spec_r <- reactive({
      nm <- input$scale %||% HEAT_SCALE_NAMES[1]; rv <- isTRUE(input$reverse)
      stops <- heat_stops(nm, rv, n = 11L)
      div <- export_r()$diverging %||% is_diverging_scale(nm)
      list(stops = stops, name = nm, reverse = rv, diverging = isTRUE(div),
           colorscale = heat_colorscale(stops))
    })
    output$dl_csv <- downloadHandler(function() paste0(pref(), "_data.csv"), function(f) {
      e <- export_r(); s <- e$summary
      out <- if (is.function(e$csv)) e$csv(s)
             else if (is.matrix(s)) data.frame(row = rownames(s), s, check.names = FALSE)
             else as.data.frame(s)
      write.csv(out, f, row.names = FALSE)
    })
    result_r <- reactive(tryCatch(xs_result_tables(export_r()$result), error = function(e) list()))
    output$has_result <- reactive(length(result_r()) > 0)
    outputOptions(output, "has_result", suspendWhenHidden = FALSE)
    output$dl_result <- downloadHandler(
      filename = function() {
        n <- length(isolate(result_r()))
        sprintf("%s_result.%s", pref(),
                if (n > 1 && bundle_archive_format() == "zip") "zip" else "csv")
      },
      content = function(f) xs_write_result(result_r(), f, pref()))
    .cap <- shiny::debounce(reactive({
      e <- tryCatch(export_r(), error = function(err) NULL)
      if (is.null(e)) return(NULL)
      list(e = e, sp = tryCatch(spec_r(), error = function(err) NULL),
           context = tryCatch(session$userData$bundle_ctx(id), error = function(err) NULL))
    }), 1500)
    observeEvent(.cap(), {
      z <- .cap(); ctxf <- session$userData$bundle_ctx
      if (is.null(z) || is.null(ctxf)) return(invisible())
      cc <- z$context
      if (is.null(cc)) return(invisible())
      gg <- tryCatch(z$e$gg(z$e$summary, z$sp$stops, z$e$ttl, z$sp$diverging),
                     error = function(err) NULL)
      bundle_capture(session, id, z$e, gg, cc$dataset, cc$layer)
    }, ignoreNULL = TRUE)

    save_one <- function(f, fmt) {
      e <- export_r(); sp <- spec_r()
      g <- e$gg(e$summary, sp$stops, e$ttl, sp$diverging)
      byod_save_figure(g, f, fmt, input$w, input$h, input$dpi)
    }
    output$dl_fig <- downloadHandler(
      filename = function() paste0(pref(), ".", input$fmt %||% "png"),
      content  = function(f) save_one(f, input$fmt %||% "png"))
    spec_r
  })
}

VOLCANO_LEVELS  <- c("Up", "Down", "n.s.")
VOLCANO_DEFAULT <- c(Up = PIG, Down = HUMAN, "n.s." = "#C7D0D8")
DEOT_LV  <- c("up in comparison", "down in comparison")
GDIFF_LV <- c("increased", "decreased")
GSEA_LV  <- c("NES > 0", "NES < 0")
ORA_LV   <- "enriched"

fb_cols <- function(levels, cols, palette = PAL_NAMES[1], splitby = NULL) {
  base <- resolve_appearance(levels, palette, splitby = splitby)
  if (is.null(cols)) return(base)
  keep <- names(cols)[names(cols) %in% levels]; base[keep] <- cols[keep]; base[levels]
}

gg_heatmap <- function(mat, stops, ttl, diverging = TRUE) {
  if (is.null(mat) || !length(mat) || !nrow(mat))
    return(ggplot() + theme_void() + labs(title = "no data"))
  df <- expand.grid(row = rownames(mat), col = colnames(mat),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df$z   <- as.vector(mat)
  df$row <- factor(df$row, levels = rev(rownames(mat)))
  df$col <- factor(df$col, levels = colnames(mat))
  sc <- if (isTRUE(diverging)) {
    lim <- max(abs(mat), na.rm = TRUE); if (!is.finite(lim) || lim == 0) lim <- 1
    scale_fill_gradientn(colours = stops, limits = c(-lim, lim), name = "z")
  } else scale_fill_gradientn(colours = stops, name = "value")
  ggplot(df, aes(col, row, fill = z)) + geom_tile() + sc +
    labs(x = NULL, y = NULL, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 12), axis.text.y = element_text(size = 8),
          panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
}

gg_violin <- function(summary, ylab, ttl, cols) {
  lv  <- levels(summary$Group) %||% as.character(summary$Group)
  raw <- attr(summary, "raw")
  base <- ggplot() + scale_fill_manual(values = cols, drop = FALSE) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = ylab, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "none", plot.title = element_text(size = 12), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = if (length(lv) > 4) 30 else 0, hjust = if (length(lv) > 4) 1 else .5))
  if (!is.null(raw) && nrow(raw)) {
    raw$Group <- factor(as.character(raw$Group), levels = lv)
    base + geom_violin(data = raw, aes(Group, value, fill = Group), colour = NA, alpha = .35, scale = "width", trim = FALSE) +
      geom_boxplot(data = raw, aes(Group, value, colour = Group), width = .15, outlier.shape = NA, fill = "white", linewidth = .5) +
      geom_point(data = raw, aes(Group, value, colour = Group),
                 position = position_jitter(width = .12, seed = trace_seed()), size = .9, alpha = .5)
  } else base + geom_crossbar(data = summary,
      aes(Group, mean, ymin = mean - se, ymax = mean + se, colour = Group), width = .5, fatten = 2)
}

gg_clinbox <- function(summary, ylab, ttl, cols) {
  lv  <- levels(summary$Group) %||% as.character(summary$Group)
  raw <- attr(summary, "raw")
  base <- ggplot() + scale_colour_manual(values = cols, drop = FALSE) + scale_fill_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = ylab, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "none", plot.title = element_text(size = 12), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = if (length(lv) > 4) 30 else 0, hjust = if (length(lv) > 4) 1 else .5))
  if (!is.null(raw) && nrow(raw)) {
    raw$Group <- factor(as.character(raw$Group), levels = lv)
    base + geom_boxplot(data = raw, aes(Group, value, colour = Group), width = .55, outlier.shape = NA, fill = "white", linewidth = .5) +
      geom_point(data = raw, aes(Group, value, colour = Group),
                 position = position_jitter(width = .14, seed = trace_seed()), size = 1, alpha = .55)
  } else base + geom_crossbar(data = summary,
      aes(Group, mean, ymin = mean - se, ymax = mean + se, colour = Group), width = .5, fatten = 2)
}

gg_clin_scatter <- function(summary, ylab, ttl, cols, xlab = "value") {
  summary$group <- factor(summary$group, levels = names(cols))
  ggplot(summary, aes(x, y, colour = group)) + geom_point(size = 2.4, alpha = .85) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(x = xlab, y = ylab, colour = NULL, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "right", plot.title = element_text(size = 11), panel.grid.minor = element_blank())
}

gg_volcano <- function(summary, ylab, ttl, cols) {
  thr <- attr(summary, "thr"); lfc <- if (is.null(thr)) 1 else unname(thr["lfc"])
  labs_keep <- attr(summary, "lab") %||% character(0)
  summary$sig <- factor(as.character(summary$sig), levels = VOLCANO_LEVELS)
  cc <- VOLCANO_DEFAULT; ok <- intersect(names(cols), names(cc)); cc[ok] <- cols[ok]
  g <- ggplot(summary, aes(logFC, y, colour = sig)) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dotted", colour = MUTED, linewidth = .4) +
    geom_point(size = 1.1, alpha = .7) +
    scale_colour_manual(values = cc, drop = FALSE, name = NULL) +
    labs(x = "Mean difference (analysis scale)", y = ylab, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(size = 11), panel.grid.minor = element_blank())
  if (length(labs_keep)) {
    lab_df <- summary[summary$Feature %in% labs_keep, , drop = FALSE]
    g <- g + ggrepel::geom_text_repel(data = lab_df, aes(label = Feature), seed = trace_seed(), size = 3, colour = INK,
      max.overlaps = Inf, min.segment.length = 0, box.padding = .3, segment.size = .3, show.legend = FALSE)
  }
  g
}

gg_mp <- function(summary, ylab, ttl, cols, n_label = 15L) {
  if (is.null(summary) || !nrow(summary)) return(ggplot() + theme_void() + labs(title = "no data"))
  if (!length(cols)) cols <- c("both up" = "#B0413E", "both down" = "#2C6E9B", discordant = "#C9D2DA")
  lev <- names(cols); summary$quad <- factor(as.character(summary$quad), levels = lev)
  lab <- head(summary[order(-abs(summary$logFC_mRNA + summary$logFC_protein)), ], n_label)
  ggplot(summary, aes(logFC_mRNA, logFC_protein, colour = quad)) +
    geom_hline(yintercept = 0, colour = MUTED, linewidth = .3) + geom_vline(xintercept = 0, colour = MUTED, linewidth = .3) +
    geom_point(size = 2, alpha = .8) +
    ggrepel::geom_text_repel(data = lab, aes(label = Feature), seed = trace_seed(), size = 3, max.overlaps = 20, show.legend = FALSE) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(x = "mRNA mean difference", y = ylab, colour = NULL, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = if (length(lev) > 1) "top" else "none", plot.title = element_text(size = 11), panel.grid.minor = element_blank())
}

gg_dimred_discrete <- function(summary, ylab, ttl, cols, ax = c("Dim1", "Dim2")) {
  if (is.null(summary) || !nrow(summary) || !all(c("x","y","group") %in% names(summary)))
    return(ggplot() + theme_void() + labs(title = gsub("·", "-", gsub("<[^>]+>", "", ttl))))
  if (!length(cols)) cols <- resolve_cols(unique(as.character(summary$group)), PAL_NAMES[1])
  summary$group <- factor(summary$group, levels = names(cols))
  ggplot(summary, aes(x, y, colour = group)) + geom_point(size = 2.2, alpha = .9) +
    scale_colour_manual(values = cols, drop = FALSE) +
    labs(x = ax[1], y = ax[2], colour = NULL, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "right", plot.title = element_text(size = 12), panel.grid.minor = element_blank())
}
gg_dimred_gradient <- function(summary, stops, ttl, diverging = FALSE, ax = c("Dim1", "Dim2"), scale_title = "value") {
  if (is.null(summary) || !nrow(summary) || !all(c("x","y","value") %in% names(summary)))
    return(ggplot() + theme_void() + labs(title = gsub("·", "-", gsub("<[^>]+>", "", ttl))))
  ggplot(summary, aes(x, y, colour = value)) + geom_point(size = 2.2, alpha = .9) +
    scale_colour_gradientn(colours = stops, name = scale_title) +
    labs(x = ax[1], y = ax[2], title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "right", plot.title = element_text(size = 12), panel.grid.minor = element_blank())
}

gg_gsea <- function(summary, ylab, ttl, cols) {
  if (is.null(summary) || !nrow(summary)) return(ggplot() + theme_void() + labs(title = "no enrichment"))
  if (!length(cols)) cols <- c("NES > 0" = "#B0413E", "NES < 0" = "#2C6E9B")
  s <- summary[order(summary$NES), ]; s$pathway <- factor(s$pathway, levels = s$pathway)
  s$sign <- factor(ifelse(s$NES > 0, "NES > 0", "NES < 0"), levels = names(cols))
  ggplot(s, aes(NES, pathway, colour = sign, size = -log10(padj + 1e-10))) +
    geom_vline(xintercept = 0, colour = MUTED, linewidth = .3) + geom_point() +
    scale_colour_manual(values = cols, drop = FALSE) + scale_size_continuous(range = c(2, 7), guide = "none") +
    labs(x = ylab, y = NULL, colour = NULL, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 12) +
    theme(legend.position = if (nlevels(s$sign) > 1) "top" else "none", plot.title = element_text(size = 11), panel.grid.minor = element_blank())
}

ora_ratio <- function(s) {
  if (!is.null(s$n_input) && all(is.finite(s$n_input)) && all(s$n_input > 0))
    list(x = s$overlap / s$n_input,
         lab = sprintf("gene ratio (of the %s submitted)",
                       format(s$n_input[1], big.mark = ",")))
  else list(x = s$overlap / pmax(1L, s$set_size), lab = "gene ratio (of the set)")
}
gg_ora <- function(summary, ylab, ttl, cols, sig_cut = 0.05, faded_alpha = 0.35,
                   fade = FALSE) {
  if (is.null(summary) || !nrow(summary)) return(ggplot() + theme_void() + labs(title = "no enrichment"))
  bar <- if (length(cols)) unname(cols[[1]]) else TEAL
  pale <- grDevices::colorRampPalette(c("#FFFFFF", bar))(10)[3]
  s <- summary[order(summary$padj), ]
  gr <- ora_ratio(s)
  s$ratio <- gr$x
  s$pathway <- factor(s$pathway, levels = rev(s$pathway))
  s$padj_plot <- pmax(s$padj, .Machine$double.xmin)
  g <- ggplot(s, aes(ratio, pathway, fill = padj_plot))
  if (isTRUE(fade)) {
    s$passes <- s$padj < sig_cut
    g <- ggplot(s, aes(ratio, pathway, fill = padj_plot, alpha = passes)) +
      geom_col(width = .7) +
      scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = faded_alpha),
                         breaks = c(TRUE, FALSE),
                         labels = c(sprintf("FDR < %g", sig_cut),
                                    sprintf("FDR >= %g", sig_cut)),
                         name = NULL, drop = FALSE) +
      guides(alpha = guide_legend(order = 2, override.aes = list(fill = bar)))
  } else g <- g + geom_col(width = .7)
  g +

    scale_fill_gradient(low = bar, high = pale, trans = "log10", name = "FDR",
                        labels = function(z) format(signif(z, 2), scientific = TRUE),
                        guide = guide_colourbar(order = 1, barwidth = unit(8, "pt"),
                                                barheight = unit(46, "pt"))) +
    labs(x = gr$lab, y = NULL, title = gsub("\u00b7", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 11), panel.grid.minor = element_blank(),
          legend.position = "right")
}

gg_de_over_time <- function(summary, ylab, ttl, cols) {
  if (is.null(summary) || !nrow(summary)) return(ggplot() + theme_void() + labs(title = "no DE results"))
  lv <- names(cols)
  long <- rbind(data.frame(tp = summary$tp, dir = "up in comparison",   n = summary$up),
                data.frame(tp = summary$tp, dir = "down in comparison", n = summary$down))
  long$dir <- factor(long$dir, levels = lv)
  if (!is.factor(long$tp)) long$tp <- factor(long$tp, levels = as.character(summary$tp))
  ggplot(long, aes(tp, n, fill = dir)) + geom_col(position = position_dodge(width = .8), width = .7) +
    scale_fill_manual(values = cols, drop = FALSE) +
    labs(x = NULL, y = ylab, fill = NULL, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(size = 11), panel.grid.minor = element_blank())
}

gg_gsva_diff <- function(summary, ylab, ttl, cols) {
  if (is.null(summary) || !nrow(summary)) return(ggplot() + theme_void() + labs(title = "no data"))
  s <- summary[order(summary$logFC), ]; s$Feature <- factor(s$Feature, levels = s$Feature)
  s$dir <- factor(ifelse(s$logFC > 0, "increased", "decreased"), levels = names(cols))
  ggplot(s, aes(logFC, Feature, fill = dir)) + geom_col(width = .7) + scale_fill_manual(values = cols, drop = FALSE) +
    labs(x = ylab, y = NULL, fill = NULL, title = gsub("·", "-", gsub("<[^>]+>", "", ttl))) +
    byod_theme_bw(base_size = 11) +
    theme(legend.position = if (nlevels(s$dir) > 1) "top" else "none", plot.title = element_text(size = 11),
          panel.grid.minor = element_blank(), axis.text.y = element_text(size = 8))
}

plot_pathway_tc <- function(species, omic, pathway, splitby = "None", palette = PAL_NAMES[1]) {
  d <- get_ds(species, omic); G <- gsva_of(species, omic)
  if (is.null(G) || !(pathway %in% rownames(G)))
    return(list(plot = emptyplot("pathway not available"),
      summary = data.frame(Timepoint = factor(levels = d$tp_levels), Group = character(),
                           mean = numeric(), se = numeric(), n = numeric()),
      splits = character(0), ylab = "", ttl = "pathway not available",
      cols = setNames(character(0), character(0))))
  val <- setNames(as.numeric(G[pathway, ]), colnames(G))
  td <- tc_data(d, val, splitby); cols <- resolve_cols(td$splits, palette, splitby)
  pval <- tryCatch(summary(aov(val ~ d$meta[colnames(G), "Timepoint"]))[[1]][["Pr(>F)"]][1], error = function(e) NA)
  ttl <- paste0(pathway, " — ", species, " ", omic, " · ANOVA p = ", signif(pval, 2))
  ylab <- "Hallmark activity (GSVA)"
  list(plot = tc_plotly(td, ylab, ttl, cols), summary = td$summary, cols = cols, ylab = ylab, ttl = ttl)
}

plot_gsva_diff <- function(species, omic, groupvar, g1, g2, paired = TRUE, at_tp = "All",
                           cols = NULL, palette = PAL_NAMES[1]) {
  d <- get_ds(species, omic); G <- gsva_of(species, omic)
  ttl <- paste0(species, " ", omic, " — Hallmark differential activity")
  xlab <- paste0("Δ activity (", g2, " − ", g1, ")")
  fail <- function(msg, tt = data.frame()) list(plot = emptyplot(msg), summary = tt, tab = tt,
    splits = character(0), ylab = xlab, ttl = msg, cols = setNames(character(0), character(0)))
  if (is.null(G)) return(fail("GSVA not available"))
  tt <- run_limma(d, groupvar, g1, g2, paired, at_tp, mat = G); if (!nrow(tt)) return(fail("too few samples"))
  tt <- tt[order(tt$logFC), ]; tt$dir <- factor(ifelse(tt$logFC > 0, "increased", "decreased"), levels = GDIFF_LV)
  splits <- intersect(GDIFF_LV, as.character(unique(tt$dir))); cols <- fb_cols(splits, cols, palette)
  tt$col <- unname(cols[as.character(tt$dir)]); yl <- factor(tt$Feature, levels = tt$Feature)
  p <- plot_ly(tt, x = ~logFC, y = yl, type = "bar", orientation = "h", marker = list(color = ~col),
      text = ~paste0("FDR=", adj.P.Val), hoverinfo = "text") %>%
    ply(yaxis = list(title = "", tickfont = list(size = 9)), xaxis = list(title = xlab), margin = list(l = 230),
        showlegend = FALSE, title = list(text = ttl, x = .02, font = list(size = 13)))
  tab <- tt[order(tt$adj.P.Val), c("Feature","logFC","t","P.Value","adj.P.Val")]
  list(plot = p, summary = tab, tab = tab, splits = splits, ylab = xlab, ttl = ttl, cols = cols)
}

run_gsea <- function(stats_named, collection, minSize = 10, maxSize = 500, seed = trace_seed(), d = NULL) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  sets <- annotation_sets(d, collection)
  if (!is.null(d)) stats_named <- stats_named[names(stats_named) %in% annotation_universe(d, collection)]
  stats_named <- stats_named[is.finite(stats_named)]
  stats_named <- sort(stats_named, decreasing = TRUE)

  scoreable <- sum(vapply(sets, function(g) sum(g %in% names(stats_named)) >= minSize, logical(1)))
  if (!scoreable) {
    out <- data.frame()
    attr(out, "declined") <- sprintf(
      paste("No %s gene set has %d or more eligible features. Ranked IDs include %s.",
            "Check the layer's organism, ID type, mapping coverage and selected GMT."),
      collection, minSize,
      paste(utils::head(names(stats_named), 3), collapse = ", "))
    return(out)
  }
  res <- tryCatch(suppressWarnings({
    fgsea::fgseaMultilevel(sets, stats_named, minSize = minSize, maxSize = maxSize, nproc = 1L) }), error = identity)
  if (inherits(res, "error")) {
    out <- data.frame()
    attr(out, "declined") <- paste("GSEA could not be completed:", conditionMessage(res))
    return(out)
  }
  if (is.null(res) || !nrow(res)) {
    out <- data.frame()
    attr(out, "declined") <- sprintf(
      "%d %s set%s met the overlap threshold, but GSEA returned no test results.",
      scoreable, collection, if (scoreable == 1L) "" else "s")
    return(out)
  }
  res <- as.data.frame(res)[, c("pathway", "NES", "pval", "padj", "size")]
  res$NES <- round(res$NES, 2); res$pval <- signif(res$pval, 3); res$padj <- signif(res$padj, 3)
  res$pathway <- gsub("_", " ", res$pathway); res[order(res$padj), ]
}
plot_gsea <- function(res, n = 20, lab = "", cols = NULL, palette = PAL_NAMES[1]) {
  if (!nrow(res)) return(list(
    plot = emptyplot(attr(res, "declined") %||% "no enrichment"),
    summary = res, splits = character(0),
    ylab = "NES", ttl = lab, cols = setNames(character(0), character(0))))
  res <- head(res[order(res$padj), ], n); res <- res[order(res$NES), ]
  res$sign <- factor(ifelse(res$NES > 0, "NES > 0", "NES < 0"), levels = GSEA_LV)
  splits <- intersect(GSEA_LV, as.character(unique(res$sign))); cols <- fb_cols(splits, cols, palette)
  res$col <- unname(cols[as.character(res$sign)])
  p <- plot_ly(res, x = ~NES, y = ~factor(pathway, levels = pathway), type = "scatter", mode = "markers",
      marker = list(size = ~pmin(22, -log10(padj + 1e-10) * 2.5 + 6), color = ~col),
      text = ~paste0(pathway, "<br>NES=", round(NES, 2), "<br>padj=", signif(padj, 2)), hoverinfo = "text") %>%
    ply(xaxis = list(title = "NES", zeroline = TRUE), yaxis = list(title = "", tickfont = list(size = 9)),
        margin = list(l = 280), title = list(text = lab, x = .02, font = list(size = 13)))
  list(plot = p, summary = res[, c("pathway","NES","pval","padj","size")], splits = splits, ylab = "NES", ttl = lab, cols = cols)
}

run_ora <- function(genes, collection, universe, d = NULL) {
  sets <- annotation_sets(d, collection)
  universe <- unique(universe)
  if (!is.null(d)) universe <- intersect(universe, annotation_universe(d, collection))
  genes <- intersect(genes, universe)
  N <- length(universe); k <- length(genes); if (k < 3) return(data.frame())

  rows <- lapply(names(sets), function(p) {
    S <- intersect(sets[[p]], universe); m <- length(S); if (m < 3) return(NULL)
    x <- length(intersect(genes, S))

    data.frame(pathway = p, overlap = x, set_size = m, n_input = k,
               p = phyper(x - 1, m, N - m, k, lower.tail = FALSE)) })
  res <- do.call(rbind, rows); if (is.null(res) || !nrow(res)) return(data.frame())
  res$padj <- p.adjust(res$p, "BH")
  res <- res[res$overlap >= 2, , drop = FALSE]; if (!nrow(res)) return(data.frame())
  res <- res[order(res$padj), ]
  res$p <- signif(res$p, 3); res$padj <- signif(res$padj, 3)
  res$pathway <- gsub("_", " ", res$pathway); res
}
plot_ora <- function(res, n = 20, lab = "", cols = NULL, palette = PAL_NAMES[1],
                     sig_cut = 0.05, faded_alpha = 0.35, fade = FALSE) {
  if (!nrow(res)) return(list(plot = emptyplot("no enrichment (need >=3 matched genes)"), summary = res,
    splits = character(0), ylab = "gene ratio", ttl = lab, cols = setNames(character(0), character(0))))
  res <- head(res, n); res <- res[order(-res$padj), ]; cols <- fb_cols(ORA_LV, cols, palette)
  bar <- unname(cols[[ORA_LV]]); pale <- grDevices::colorRampPalette(c("#FFFFFF", bar))(10)[3]
  gr <- ora_ratio(res); res$ratio <- gr$x

  nlp  <- -log10(pmax(res$padj, .Machine$double.xmin))
  ramp <- grDevices::colorRampPalette(c(pale, bar))(64)
  idx  <- if (diff(range(nlp, finite = TRUE)) > 0)
            1 + round(63 * (nlp - min(nlp)) / diff(range(nlp))) else rep(64L, nrow(res))
  op <- if (isTRUE(fade)) ifelse(res$padj < sig_cut, 1, faded_alpha) else rep(1, nrow(res))
  note <- if (isTRUE(fade)) ifelse(res$padj < sig_cut, "", " - below the cut-off") else rep("", nrow(res))
  res$hover <- paste0(res$overlap, "/", res$n_input, " submitted, set size ", res$set_size,
                      "<br>FDR ", signif(res$padj, 3), note)
  p <- plot_ly(res, x = ~ratio, y = ~factor(pathway, levels = pathway), type = "bar",
      orientation = "h", marker = list(color = ramp[idx], opacity = op),
      text = ~hover, hoverinfo = "text") %>%
    ply(xaxis = list(title = gr$lab), yaxis = list(title = "", tickfont = list(size = 9)),
        margin = list(l = 280), title = list(text = lab, x = .02, font = list(size = 13)))
  list(plot = p, summary = res[, c("pathway", "overlap", "set_size", "n_input", "p", "padj")],
       splits = ORA_LV, ylab = gr$lab, ttl = lab, cols = cols)
}

plot_clinical_assoc <- function(species, omic, feature, ftype, clin, at_tp = "All", cols = NULL) {
  d <- get_ds(species, omic)
  fail <- function(msg) list(plot = emptyplot(msg), summary = data.frame(),
    kind = "discrete", splits = character(0), ttl = msg, ylab = "", mode = "none", xlab = "")
  if (ftype == "Pathway") { G <- gsva_of(species, omic)
    if (is.null(G) || !(feature %in% rownames(G))) return(fail("pathway not available"))
    samp <- colnames(G); val <- as.numeric(G[feature, ]); ylab <- paste0(feature, " (GSVA)")
  } else { v0 <- expr_vec(d, feature); if (is.null(v0)) return(fail(paste(feature, "not detected")))
    samp <- colnames(d$expr); val <- v0; ylab <- paste(feature, omic) }
  mt <- d$meta[samp, , drop = FALSE]
  df <- data.frame(sample = samp, val = val, clin = mt[[clin]],
                   tp = factor(mt$Timepoint, levels = d$tp_levels), stringsAsFactors = FALSE)
  if (at_tp != "All") df <- df[df$tp == at_tp, ]
  df <- df[!is.na(df$clin) & !is.na(df$val), ]; if (nrow(df) < 3) return(fail("too few samples"))
  suffix <- if (at_tp != "All") paste0(" · @", at_tp) else " · all timepoints"
  if (is.numeric(df$clin)) {                                   # scatter, colour = Timepoint (discrete)
    df$tp <- droplevels(df$tp); lv <- levels(df$tp)
    base <- setNames(as.character(TP_COL[lv]), lv); if (!is.null(cols)) { ok <- intersect(names(cols), lv); base[ok] <- cols[ok] }; cols <- base
    ct <- suppressWarnings(cor.test(df$clin, df$val))
    ttl <- paste0(feature, " vs ", clin, " · r = ", round(ct$estimate, 2), ", p = ", signif(ct$p.value, 2), suffix)
    p <- plot_ly(df, x = ~clin, y = ~val, color = ~tp, colors = unname(cols[lv]),
        type = "scatter", mode = "markers", marker = list(size = 7), text = ~paste0(sample, "<br>", tp), hoverinfo = "text") %>%
      ply(xaxis = list(title = clin), yaxis = list(title = ylab), title = list(text = ttl, x = .02, font = list(size = 13)))
    return(list(plot = p, summary = data.frame(sample = df$sample, x = df$clin, y = df$val,
        group = as.character(df$tp), stringsAsFactors = FALSE), kind = "discrete", splits = lv,
        ttl = ttl, ylab = ylab, mode = "scatter", xlab = clin, cols = cols))
  }
  df$clin <- droplevels(factor(df$clin)); lv <- levels(df$clin)                  # box, colour = clinical (discrete)
  base <- setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = NULL)), lv)
  if (!is.null(cols)) { ok <- intersect(names(cols), lv); base[ok] <- cols[ok] }; cols <- base
  pv <- tryCatch(kruskal.test(val ~ clin, df)$p.value, error = function(e) NA)
  ttl <- paste0(feature, " by ", clin, " · Kruskal p = ", signif(pv, 2), suffix)
  agg <- function(f) tapply(df$val, df$clin, f)
  summary <- data.frame(Group = factor(lv, levels = lv),
    mean = as.numeric(agg(function(x) mean(x, na.rm = TRUE))[lv]),
    sd = as.numeric(agg(function(x) sd(x, na.rm = TRUE))[lv]),
    median = as.numeric(agg(function(x) median(x, na.rm = TRUE))[lv]),
    n = as.integer(agg(function(x) sum(!is.na(x)))[lv]), row.names = NULL, stringsAsFactors = FALSE)
  summary$se <- summary$sd / sqrt(pmax(summary$n, 1))
  attr(summary, "raw") <- data.frame(Group = factor(as.character(df$clin), levels = lv), value = df$val)
  p <- plot_ly(df, x = ~clin, y = ~val, color = ~clin, colors = unname(cols[lv]),
      type = "box", boxpoints = "all", jitter = .3, pointpos = 0, marker = list(size = 4)) %>%
    ply(xaxis = list(title = clin), yaxis = list(title = ylab), showlegend = FALSE,
        title = list(text = ttl, x = .02, font = list(size = 13)))
  list(plot = p, summary = summary, kind = "discrete", splits = lv, ttl = ttl, ylab = ylab, mode = "box", xlab = clin, cols = cols)
}

contrast_vars <- function(d) setdiff(group_vars(d), c("Timepoint", "Timepoint_WI"))

plot_split_trajectory <- function(species, omic, feature, ftype, var, g1, g2, cols = NULL, palette = PAL_NAMES[1]) {
  d <- get_ds(species, omic)
  fail <- function(msg) list(plot = emptyplot(msg),
    summary = data.frame(Timepoint = factor(levels = d$tp_levels), Group = character(),
                         mean = numeric(), se = numeric(), n = numeric()),
    splits = character(0), ylab = "", ttl = msg, cols = setNames(character(0), character(0)))
  if (ftype == "Pathway") { G <- gsva_of(species, omic)
    if (is.null(G) || !(feature %in% rownames(G))) return(fail("pathway not available"))
    samp <- colnames(G); val <- setNames(as.numeric(G[feature, ]), samp); ylab <- paste0(feature, " (GSVA)")
  } else { v0 <- expr_vec(d, feature); if (is.null(v0)) return(fail(paste(feature, "not detected")))
    samp <- colnames(d$expr); val <- setNames(v0, samp); ylab <- paste(feature, omic) }
  gf <- grp_factor(d, var); names(gf) <- rownames(d$meta)
  tp <- factor(d$meta$Timepoint, levels = d$tp_levels); names(tp) <- rownames(d$meta)
  df <- data.frame(val = val[samp], grp = as.character(gf[samp]), tp = tp[samp], stringsAsFactors = FALSE)
  df <- df[df$grp %in% c(g1, g2) & !is.na(df$grp) & !is.na(df$tp) & !is.na(df$val), ]
  if (!nrow(df)) return(fail("no data for these groups"))
  df$gt <- interaction(df$grp, df$tp, sep = "@@", drop = TRUE)
  s <- summ_stats(df$val, df$gt)
  pr <- do.call(rbind, strsplit(as.character(s$group), "@@", fixed = TRUE))
  splits <- intersect(c(g1, g2), unique(pr[, 1]))
  summary <- data.frame(Timepoint = factor(pr[, 2], levels = d$tp_levels),
                        Group = factor(pr[, 1], levels = splits),
                        mean = s$mean, se = s$se, n = s$n, stringsAsFactors = FALSE)
  summary <- summary[order(summary$Group, summary$Timepoint), ]
  if (is.null(cols)) cols <- resolve_appearance(splits, palette, splitby = var) else cols <- cols[splits]
  ttl <- paste0("<b>", feature, "</b> — ", g2, " vs ", g1, " across perfusion (", species, " ", omic, ")")
  list(plot = tc_plotly(list(summary = summary, splits = splits), ylab, ttl, cols),
       summary = summary, splits = splits, ylab = ylab, ttl = ttl, cols = cols)
}

de_over_time <- function(d, var, g1, g2, fdr = 0.05, lfc = 1) {
  paired <- var == "WI_status"
  rows <- lapply(d$tp_levels, function(tp) {
    tt <- tryCatch(run_limma(d, var, g1, g2, paired, at_tp = tp), error = function(e) NULL)
    if (is.null(tt) || !nrow(tt)) return(data.frame(tp = tp, up = 0, down = 0))
    data.frame(tp = tp, up = sum(tt$adj.P.Val < fdr & tt$logFC >= lfc, na.rm = TRUE),
               down = sum(tt$adj.P.Val < fdr & tt$logFC <= -lfc, na.rm = TRUE)) })
  res <- do.call(rbind, rows); res$tp <- factor(res$tp, levels = d$tp_levels); res
}
plot_de_over_time <- function(res, lab, cols = NULL, palette = PAL_NAMES[1]) {
  ylab <- "# DE genes (FDR<0.05, |log2FC|>=1)"
  if (is.null(res) || !nrow(res)) return(list(plot = emptyplot("no DE results"), summary = data.frame(),
    splits = character(0), ylab = ylab, ttl = lab, cols = setNames(character(0), character(0))))
  splits <- DEOT_LV; cols <- fb_cols(splits, cols, palette)
  p <- plot_ly(res, x = ~tp) %>%
    add_bars(y = ~up, name = "up in comparison", marker = list(color = unname(cols[["up in comparison"]]))) %>%
    add_bars(y = ~down, name = "down in comparison", marker = list(color = unname(cols[["down in comparison"]]))) %>%
    ply(barmode = "group", xaxis = list(title = ""), yaxis = list(title = ylab),
        legend = list(orientation = "h", y = 1.1, x = .42), title = list(text = lab, x = .02, font = list(size = 13)))
  list(plot = p, summary = res, splits = splits, ylab = ylab, ttl = lab, cols = cols)
}

gsea_over_time <- function(d, var, g1, g2, collection = "Hallmark") {
  paired <- var == "WI_status"
  mats <- lapply(d$tp_levels, function(tp) {
    tt <- tryCatch(run_limma(d, var, g1, g2, paired, at_tp = tp), error = function(e) NULL)
    if (is.null(tt) || !nrow(tt)) return(NULL)
    r <- run_gsea(setNames(tt$t, tt$Feature), collection); if (!nrow(r)) return(NULL)
    setNames(r$NES, r$pathway) })
  names(mats) <- d$tp_levels; mats <- mats[!vapply(mats, is.null, logical(1))]
  if (!length(mats)) return(NULL)
  allp <- unique(unlist(lapply(mats, names)))
  M <- sapply(mats, function(m) m[allp]); rownames(M) <- allp; M[is.na(M)] <- 0; M
}
plot_gsea_over_time <- function(M, lab, colorscale = NULL) {
  if (is.null(M)) return(list(plot = emptyplot("no enrichment / too few samples"), mat = NULL, ttl = lab))
  M <- M[order(-rowMeans(abs(M))), , drop = FALSE]; M <- head(M, 30)
  if (nrow(M) > 2) M <- M[hclust(dist(M))$order, , drop = FALSE]
  cs <- colorscale %||% list(c(0, HUMAN), c(.5, "#FFFFFF"), c(1, PIG))
  p <- plot_ly(x = colnames(M), y = rownames(M), z = M, type = "heatmap", colorscale = cs, zmid = 0,
    hovertemplate = "%{y}<br>%{x}: NES %{z:.2f}<extra></extra>") %>%
    ply(xaxis = list(title = ""), yaxis = list(title = "", tickfont = list(size = 9)), margin = list(l = 270),
        title = list(text = lab, x = .02, font = list(size = 12)))
  list(plot = p, mat = M, ttl = lab)
}

byod_rowvars <- function(M) {
  n  <- rowSums(!is.na(M)); mu <- rowMeans(M, na.rm = TRUE)
  v  <- ifelse(n > 1, rowSums((M - mu)^2, na.rm = TRUE) / (n - 1), 0)
  v[!is.finite(v)] <- 0; v
}

byod_read <- function(path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path))
    stop("Input file does not exist: ", paste(path, collapse = ", "), call. = FALSE)
  if (requireNamespace("data.table", quietly = TRUE))
    return(as.data.frame(data.table::fread(file = path, header = TRUE, check.names = FALSE, colClasses = list(character = 1L)),
                         stringsAsFactors = FALSE, check.names = FALSE))
  con <- gzfile(path, "rt"); header <- tryCatch(readLines(con, n = 1L, warn = FALSE), finally = close(con))
  sep <- if (length(header) && grepl("\t", header, fixed = TRUE)) "\t" else ","
  n <- ncol(read.delim(path, sep = sep, header = TRUE, nrows = 0, check.names = FALSE))
  classes <- rep(NA_character_, n); classes[1] <- "character"
  read.delim(path, sep = sep, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE, colClasses = classes)
}

BYOD_NORM_METHODS <- c("Auto-detect", "None", "log2(x+1)", "CPM to log2", "Quantile (log2)")
byod_diagnose <- function(M) {
  v <- suppressWarnings(as.numeric(M)); v <- v[is.finite(v)]
  if (!length(v)) return(list(min = NA, max = NA, pct_zero = NA, integer = FALSE,
                              negative = FALSE, scale = "empty", recommend = "None"))
  mx <- max(v); mn <- min(v); fz <- mean(v == 0)
  intish <- mean(abs(v - round(v)) < 1e-8) > 0.98
  neg <- mn < -1e-6
  scale <- if (neg) "centered / log-ratio / z-scored"
           else if (intish && mx > 30) "raw counts"
           else if (mx > 100) "linear intensities (not log-transformed)"
           else if (mx > 30) "large positive values (possibly un-logged)"
           else "log-scale (looks normalized)"
  recommend <- if (neg) "None" else if (intish && mx > 30) "CPM to log2" else if (mx > 30) "log2(x+1)" else "None"
  list(min = round(mn, 2), max = round(mx, 2), pct_zero = round(100 * fz, 1),
       integer = intish, negative = neg, scale = scale, recommend = recommend)
}

BYOD_OMICS_TYPES <- list(
  transcriptomics = list(
    label = "Transcriptomics (RNA-seq / array)", name = "mRNA",
    unit = "raw counts", norm = "Auto-detect", id_map = "none",
    feature = "gene", genesets = TRUE,
    hint = "Rows are genes. Declare their identifier type and source organism under Species annotation & gene sets."),
  proteomics = list(
    label = "Proteomics (LC-MS/MS)", name = "Protein",
    unit = "log2 intensity", norm = "None", id_map = "none",
    feature = "gene", genesets = TRUE,
    hint = "Rows are proteins. Human UniProt accessions can be annotated automatically; other accessions and platform IDs need a custom mapping or native-ID GMT."),
  mirna = list(
    label = "miRNA / small RNA", name = "miRNA",
    unit = "raw counts", norm = "Auto-detect", id_map = "none",
    feature = "mirna", genesets = FALSE,
    hint = "Rows are mature miRNAs. Gene-set collections are keyed on protein-coding gene symbols, so pathway and enrichment modules do not apply; everything else does."),
  metabolomics = list(
    label = "Metabolomics", name = "Metabolite",
    unit = "linear intensity", norm = "log2(x+1)", id_map = "none",
    feature = "metabolite", genesets = FALSE,
    hint = "Rows are metabolites. Peak intensities are usually linear and heavily skewed, so a log transform is the default. Gene-set modules do not apply."),
  lipidomics = list(
    label = "Lipidomics", name = "Lipid",
    unit = "linear intensity", norm = "log2(x+1)", id_map = "none",
    feature = "lipid", genesets = FALSE,
    hint = "Rows are lipid species. Same handling as metabolomics."),
  methylation = list(
    label = "DNA methylation", name = "Methylation",
    unit = "beta value (0-1)", norm = "None", id_map = "none",
    feature = "cpg", genesets = FALSE,
    hint = "Beta values are bounded 0-1 and must NOT be log-transformed. M-values are already on a log scale; declare whichever you have."),
  atac = list(
    label = "Chromatin accessibility (ATAC)", name = "ATAC",
    unit = "raw counts", norm = "Auto-detect", id_map = "none",
    feature = "region", genesets = FALSE,
    hint = "Rows are peaks or regions. Gene-set modules apply only if you have already summarised peaks to gene symbols."),
  other = list(
    label = "Other / custom", name = "Layer",
    unit = "expression", norm = "Auto-detect", id_map = "none",
    feature = "other", genesets = FALSE,
    hint = "Any numeric feature x sample matrix. Declare the unit accurately; it decides which relative measures are meaningful."))

BYOD_OMICS_CHOICES <- setNames(names(BYOD_OMICS_TYPES),
                               vapply(BYOD_OMICS_TYPES, `[[`, character(1), "label"))

byod_genesets_apply <- function(d, collection = "Hallmark", min_hit = 50L) {
  any(vapply(annotation_sets(d, collection), function(s)
    length(intersect(s, rownames(d$expr))) >= 3L, logical(1)))
}

byod_detect_ids <- function(rn, omics = NULL) {
  rn <- as.character(rn); rn <- rn[!is.na(rn)]

  if (!length(rn))
    return(list(kind = "none", versioned = FALSE, label = "no features",
                mappable = NA_character_, species = NA_character_))
  pct <- function(re) mean(grepl(re, rn))
  ver <- pct("^ENS[A-Z]*G[0-9]{6,}\\.[0-9]+$")
  if (pct("^ENS[A-Z]*G[0-9]{6,}") > 0.5) {
    hu <- pct("^ENSG") > 0.5
    return(list(kind = "ensembl", versioned = ver > 0.25,
                label = sprintf("%s Ensembl gene IDs%s", if (hu) "human" else "non-human",
                                if (ver > 0.25) " with version suffixes" else ""),

                mappable = if (hu) "ensembl_hs" else NA_character_,
                species = if (hu) "human" else "non-human"))
  }
  if (pct("^AT[1-5CM]G[0-9]{5}") > 0.5)
    return(list(kind = "agi", versioned = FALSE, label = "Arabidopsis AGI locus codes",
                mappable = NA_character_, species = "Arabidopsis"))

  if (pct("^[0-9]+_at$") > 0.5)
    return(list(kind = "entrez_at", versioned = FALSE,
                label = "numeric probe-set IDs (may be Entrez gene IDs from a custom CDF)",
                mappable = "entrez_hs", species = "human"))

  if (pct("^[0-9]+$") > 0.8)
    return(list(kind = "numeric_probe", versioned = FALSE,
                label = "numeric probe IDs (no gene names)",
                mappable = NA_character_, species = NA_character_))
  if (pct("_at$") > 0.5)
    return(list(kind = "probe", versioned = FALSE, label = "microarray probe IDs",
                mappable = NA_character_, species = NA_character_))
  if (pct("^[a-z]{3}-(let|miR|mir)") > 0.5)
    return(list(kind = "mirna", versioned = FALSE, label = "miRBase mature miRNA names",
                mappable = NA_character_, species = NA_character_))

  if (pct("^(cg|ch|rs)[0-9]{6,}") > 0.5)
    return(list(kind = "cpg", versioned = FALSE,
                label = "Illumina methylation probe IDs (CpG sites, not genes)",
                mappable = NA_character_, species = NA_character_))

  if (pct("^(chr)?[0-9XYMT]+[.:_-][0-9]+[.:_-][0-9]+$") > 0.5)
    return(list(kind = "region", versioned = FALSE,
                label = "genomic intervals (peaks or regions, not genes)",
                mappable = NA_character_, species = NA_character_))
  if (pct("^A_[0-9]{2}_P[0-9]+$") > 0.5)
    return(list(kind = "probe", versioned = FALSE,
                label = "Agilent probe IDs (not gene names)",
                mappable = NA_character_, species = NA_character_))

  if (pct("^seq[.][0-9]+[.][0-9]+$") > 0.5)
    return(list(kind = "somamer", versioned = FALSE,
                label = "SomaScan aptamer IDs (not gene symbols)",
                mappable = NA_character_, species = NA_character_))
  if (pct("^[A-Za-z][A-Za-z0-9-]*[.][0-9]{3,6}[.][0-9]{1,3}$") > 0.5)
    return(list(kind = "somamer_symbol", versioned = TRUE,
                label = "gene symbols with SomaScan aptamer suffixes (SYMBOL.seqid.version) - enable \"strip trailing numeric suffixes\" to recover the symbol",
                mappable = NA_character_, species = NA_character_))
  if (pct("^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$") > 0.5)
    return(list(kind = "uniprot", versioned = FALSE, label = "UniProt accessions",
                mappable = NA_character_, species = NA_character_))

  if (pct("^(chr)?[0-9XYMT]{1,2}[.:_-][0-9]{4,}$") > 0.5)
    return(list(kind = "locus", versioned = FALSE,
                label = "genomic positions (single sites, not genes)",
                mappable = NA_character_, species = NA_character_))

  if (pct("^[0-9]+([.][0-9]+)?@[0-9]+([.][0-9]+)?$|^[0-9]+[.][0-9]+_[0-9]+([.][0-9]+)?$") > 0.5)
    return(list(kind = "mz_rt", versioned = FALSE,
                label = "untargeted MS features (m/z @ retention time, not compound names)",
                mappable = NA_character_, species = NA_character_))

  symbolish <- pct("^[A-Za-z0-9][A-Za-z0-9._-]{0,17}$")

  if (!is.null(omics) && omics %in% c("metabolomics", "lipidomics")) {
    what <- if (identical(omics, "lipidomics")) "lipid species names" else "metabolite names"
    return(list(kind = "compound", versioned = FALSE,
                label = sprintf("%s (not gene symbols)", what),
                mappable = NA_character_, species = NA_character_))
  }
  if (symbolish < 0.5)
    return(list(kind = "unknown", versioned = FALSE,
                label = "unrecognised identifiers (not gene symbols)",
                mappable = NA_character_, species = NA_character_))
  list(kind = "symbol", versioned = FALSE, label = "gene symbols",
       mappable = NA_character_, species = NA_character_)
}

byod_id_coverage <- function(rn, collection = "Hallmark") {
  sets <- tryCatch(GENESETS[[collection]], error = function(e) NULL)
  if (is.null(sets)) return(list(n = length(rn), hit = NA_integer_, pct = NA_real_))
  u <- toupper(unique(unlist(sets)))
  hit <- sum(toupper(rn) %in% u)
  list(n = length(rn), hit = hit, pct = round(100 * hit / max(length(rn), 1), 1))
}

byod_map_ids <- function(M, how = "none") {
  if (!is.null(how) && how != "none") stop("Use layer annotation with an explicit organism; expression IDs are no longer rewritten.")
  list(M = M, applied = FALSE, msg = "identifiers retained", n_in = nrow(M), n_out = nrow(M))
}

BYOD_DUP_RULES <- c("Keep the highest-mean row" = "mean",
                    "Keep the most variable row" = "var",
                    "Sum the rows" = "sum",
                    "Average the rows" = "avg")

byod_clean_defaults <- function()
  list(trim_names = TRUE, strip_version = FALSE, upper_names = FALSE,
       dup_rule = "mean", min_value = 0, min_samples = 0, max_missing = 100,
       drop_constant = FALSE, drop_failed = FALSE, failed_frac = 0.01)

byod_clean_matrix <- function(M, opt = NULL, is_counts = FALSE) {
  o <- utils::modifyList(byod_clean_defaults(), opt %||% list())
  steps <- list()
  note <- function(what, n_before, n_after, detail = "")
    steps[[length(steps) + 1L]] <<- list(step = what, removed = n_before - n_after,
                                         kept = n_after, detail = detail)

  rn <- rownames(M)
  if (isTRUE(o$trim_names)) rn <- trimws(rn)

  if (isTRUE(o$strip_version)) {
    repeat {
      nxt <- sub("[.][0-9]+$", "", rn)
      if (identical(nxt, rn)) break
      rn <- nxt
    }
  }
  if (isTRUE(o$upper_names)) rn <- toupper(rn)
  n_renamed <- sum(rn != rownames(M), na.rm = TRUE)
  rownames(M) <- rn
  if (n_renamed) note("tidy feature names", nrow(M), nrow(M),
                      sprintf("%s name(s) changed", format(n_renamed, big.mark = ",")))

  bad <- is.na(rn) | !nzchar(rn)
  if (any(bad)) { n0 <- nrow(M); M <- M[!bad, , drop = FALSE]
                  note("drop unnamed features", n0, nrow(M)) }

  if (anyDuplicated(rownames(M))) {
    n0 <- nrow(M)
    M <- switch(o$dup_rule,
      var = { v <- byod_rowvars(M); ord <- order(v, decreasing = TRUE)
              Z <- M[ord, , drop = FALSE]; Z[!duplicated(rownames(Z)), , drop = FALSE] },
      sum = { g <- rownames(M); Z <- rowsum(M, g, na.rm = TRUE); Z },
      avg = { g <- rownames(M); Z <- rowsum(M, g, na.rm = TRUE) / as.vector(table(g)[rownames(rowsum(M, g))]); Z },
      { ord <- order(rowMeans(M, na.rm = TRUE), decreasing = TRUE)   # "mean", the default
        Z <- M[ord, , drop = FALSE]; Z[!duplicated(rownames(Z)), , drop = FALSE] })
    note("collapse duplicate features", n0, nrow(M),
         sprintf("rule: %s", names(BYOD_DUP_RULES)[match(o$dup_rule, BYOD_DUP_RULES)] %||% o$dup_rule))
  }

  if (is.finite(o$max_missing) && o$max_missing < 100) {
    pct_na <- 100 * rowMeans(is.na(M))
    n0 <- nrow(M); M <- M[pct_na <= o$max_missing, , drop = FALSE]
    note("drop features above the missingness limit", n0, nrow(M),
         sprintf("limit: %.0f%% missing", o$max_missing))
  }

  if (o$min_samples > 0) {
    ok <- rowSums(M >= o$min_value, na.rm = TRUE) >= o$min_samples
    n0 <- nrow(M); M <- M[ok, , drop = FALSE]
    note("drop low-signal features", n0, nrow(M),
         sprintf("kept features with value >= %g in >= %d sample(s)", o$min_value, o$min_samples))
  }

  if (isTRUE(o$drop_constant)) {
    n0 <- nrow(M); M <- M[byod_rowvars(M) > 0, , drop = FALSE]
    note("drop features with no variance", n0, nrow(M))
  }

  dropped_samples <- character(0)
  if (isTRUE(o$drop_failed) && is_counts && ncol(M) > 3) {
    ls <- colSums(M, na.rm = TRUE); med <- stats::median(ls)
    if (is.finite(med) && med > 0) {
      keep <- ls / med >= o$failed_frac
      dropped_samples <- colnames(M)[!keep]
      if (any(!keep)) { n0 <- ncol(M); M <- M[, keep, drop = FALSE]
        steps[[length(steps) + 1L]] <- list(step = "drop failed libraries",
          removed = n0 - ncol(M), kept = ncol(M),
          detail = sprintf("below %g x the median library size: %s", o$failed_frac,
                           paste(dropped_samples, collapse = ", "))) }
    }
  }

  list(M = M, steps = steps, dropped_samples = dropped_samples, opt = o,
       any_applied = length(steps) > 0)
}

byod_clean_summary <- function(cl) {
  if (is.null(cl) || !length(cl$steps)) return("no cleaning steps were applied")
  vapply(cl$steps, function(s)
    sprintf("%s: %s removed, %s kept%s", s$step,
            format(s$removed, big.mark = ","), format(s$kept, big.mark = ","),
            if (nzchar(s$detail)) paste0(" (", s$detail, ")") else ""),
    character(1))
}

byod_normalize <- function(M, method = "Auto-detect", diag = NULL) {
  if (identical(method, "Auto-detect")) { if (is.null(diag)) diag <- byod_diagnose(M); method <- diag$recommend }
  applied <- method
  out <- switch(method,
    "log2(x+1)"       = log2(pmax(M, 0) + 1),
    "CPM to log2"     = { cs <- colSums(M, na.rm = TRUE); cs[cs == 0 | is.na(cs)] <- 1
                          log2(sweep(pmax(M, 0), 2, cs, "/") * 1e6 + 1) },
    "Quantile (log2)" = { L <- log2(pmax(M, 0) + 1)
                          if (requireNamespace("limma", quietly = TRUE)) limma::normalizeQuantiles(L) else L },
    "None"            = M,
    { applied <- "None"; M })
  list(M = out, applied = applied)
}

byod_combat <- function(M, batch, protect = NULL) {
  if (!requireNamespace("sva", quietly = TRUE)) return(list(M = M, applied = FALSE, msg = "sva not installed"))
  b <- factor(as.character(batch)); tab <- table(b)

  singles <- names(tab)[tab < 2]
  if (length(singles))
    return(list(M = M, applied = FALSE, msg = sprintf("not applied — %d batch(es) have a single sample (%s). Exclude those samples or merge batches; ComBat was skipped to avoid mixing corrected and uncorrected scales.",
                length(singles), paste(head(singles, 6), collapse = ", "))))
  if (nlevels(droplevels(b)) < 2 || ncol(M) < 4)
    return(list(M = M, applied = FALSE, msg = "need >=2 batches with >=2 samples each"))
  mod <- NULL
  if (!is.null(protect)) { pf <- factor(as.character(protect))
    if (nlevels(pf) >= 2) mod <- tryCatch(model.matrix(~ pf), error = function(e) NULL) }
  cor <- tryCatch(suppressMessages(sva::ComBat(dat = as.matrix(M), batch = b, mod = mod,
                  par.prior = TRUE, prior.plots = FALSE)), error = function(e) NULL)
  if (is.null(cor)) return(list(M = M, applied = FALSE, msg = "ComBat failed (batch may be confounded with the protected timepoint variable)"))
  list(M = cor, applied = TRUE, msg = sprintf("%d batches corrected (all %d samples)", nlevels(b), ncol(M)))
}

byod_order_timepoints <- function(labels) {
  if (length(labels) < 2L) return(labels)
  z <- tolower(trimws(labels))
  number <- "[+-]?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:e[+-]?[0-9]+)?"
  units <- c(s = 1, sec = 1, secs = 1, second = 1, seconds = 1,
             min = 60, mins = 60, minute = 60, minutes = 60,
             h = 3600, hr = 3600, hrs = 3600, hour = 3600, hours = 3600,
             d = 86400, day = 86400, days = 86400,
             w = 604800, wk = 604800, wks = 604800, week = 604800, weeks = 604800)
  elapsed <- vapply(z, function(label) {
    m <- regmatches(label, regexec(paste0("^(", number, ")\\s*([a-z]+)$"), label, perl = TRUE))[[1]]
    if (length(m) == 3L && m[3] %in% names(units))
      return(as.numeric(m[2]) * unname(units[m[3]]))
    m <- regmatches(label, regexec(paste0("^([a-z]+)\\s*(", number, ")$"), label, perl = TRUE))[[1]]
    if (length(m) == 3L && m[2] %in% names(units))
      return(as.numeric(m[3]) * unname(units[m[2]]))
    NA_real_
  }, numeric(1))
  if (all(is.finite(elapsed))) return(labels[order(elapsed)])

  parts <- regmatches(z, regexec(paste0("^([^0-9]*?)(", number, ")([^0-9]*)$"), z, perl = TRUE))
  if (all(lengths(parts) == 4L) &&
      length(unique(vapply(parts, `[`, character(1), 2))) == 1L &&
      length(unique(vapply(parts, `[`, character(1), 4))) == 1L) {
    values <- as.numeric(vapply(parts, `[`, character(1), 3))
    if (all(is.finite(values))) return(labels[order(values)])
  }
  labels
}

byod_prep_meta <- function(meta_df, tp_col = NULL, tp_order = NULL, subject_col = NULL) {
  if (ncol(meta_df) < 2) stop("Metadata file needs a sample-ID column plus at least one variable column.")
  meta <- meta_df[, -1, drop = FALSE]; rownames(meta) <- as.character(meta_df[[1]])
  if (is.null(tp_col) || !(tp_col %in% names(meta))) tp_col <- names(meta)[1]

  tpv <- trimws(as.character(meta[[tp_col]]))
  tpv[!nzchar(tpv) | tpv %in% c("NA", "na", "N/A", "null", "NULL", "-")] <- NA_character_
  meta[[tp_col]] <- tpv
  tp_order <- if (is.null(tp_order)) NULL else {
    z <- trimws(as.character(tp_order)); unique(z[nzchar(z) & !is.na(z)])
  }
  lv  <- if (!is.null(tp_order)) tp_order[tp_order %in% tpv] else unique(tpv[!is.na(tpv)])
  if (is.null(tp_order)) lv <- byod_order_timepoints(lv)
  meta$Timepoint <- factor(tpv, levels = lv)
  untimed <- is.na(meta$Timepoint)
  excluded <- NULL
  if (any(untimed)) {
    labs <- tpv[untimed]; labs[is.na(labs) | !nzchar(trimws(labs))] <- "(blank)"
    excluded <- list(n = sum(untimed), samples = rownames(meta)[untimed],
                     labels = sort(unique(labs)))
    meta <- meta[!untimed, , drop = FALSE]
    if (!nrow(meta))
      stop("No profile has a timepoint that is one of the declared levels. Check the timepoint column and its declared order.")
    meta$Timepoint <- droplevels(meta$Timepoint); lv <- levels(meta$Timepoint)
  }
  if (!is.null(subject_col) && !(subject_col %in% names(meta))) subject_col <- NULL
  others <- setdiff(names(meta), c(tp_col, "Timepoint", subject_col))
  cat_cols <- character(0); num_cols <- character(0)
  varying <- list()
  for (cc in others) {
    x <- meta[[cc]]
    const_ok <- TRUE
    if (!is.null(subject_col)) {
      by_subj <- tapply(as.character(x), as.character(meta[[subject_col]]),
                        function(v) length(unique(v[!is.na(v) & nzchar(v)])))
      bad <- names(by_subj)[!is.na(by_subj) & by_subj > 1]
      if (length(bad)) { const_ok <- FALSE; varying[[cc]] <- bad }
    }
    if (!const_ok) next
    if (is.numeric(x)) num_cols <- c(num_cols, cc)
    else { k <- length(unique(x[!is.na(x)])); if (k >= 2 && k <= 20) cat_cols <- c(cat_cols, cc) }
  }
  list(meta = meta, tp_levels = lv, subject = subject_col, excluded_untimed = excluded,
       clinical_cat = if (length(cat_cols)) cat_cols else NULL,
       clinical_num = if (length(num_cols)) num_cols else NULL,

       subject_varying = if (length(varying)) varying else NULL)
}

byod_build_layer <- function(expr_df, prep, unit = "expression", name = "layer",
                             norm = "Auto-detect", batch_col = NULL, max_genes_pca = 2000L,
                             id_map = "none", clean = NULL, omics = "other",
                             organism = NULL, annotation = NULL) {
  restore_rng <- trace_rng_scope(); on.exit(restore_rng(), add = TRUE)
  organism <- annotation_organism(organism)
  if (ncol(expr_df) < 2) stop("Expression file needs a feature column plus at least one sample column.")
  feat <- as.character(expr_df[[1]])
  X <- expr_df[, -1, drop = FALSE]

  M <- suppressWarnings(vapply(X, function(col)
    if (is.numeric(col)) as.numeric(col) else as.numeric(as.character(col)),
    numeric(nrow(X))))

  .na_mark <- function(z) {
    z <- toupper(trimws(as.character(z)))
    is.na(z) | !nzchar(z) |
      z %in% c("NA", "N/A", "#N/A", "NAN", "NULL", "NIL", "NONE", ".", "-", "--", "?",
               "INF", "-INF", "INFINITY", "-INFINITY") |
      startsWith(z, "#")            # #DIV/0!, #VALUE!, #REF!, #NAME?
  }

  bad <- names(X)[vapply(seq_along(X), function(j) {
    idx <- which(is.na(M[, j]))
    length(idx) > 0L && any(!.na_mark(X[[j]][idx]))
  }, logical(1))]
  if (length(bad)) stop(sprintf("Non-numeric text in sample column(s): %s. Expression matrices must contain only numbers (put feature IDs in the first column; use a blank or NA for missing values).",
                                paste(head(bad, 6), collapse = ", ")))
  if (is.matrix(M)) colnames(M) <- names(X) else M <- matrix(M, nrow = nrow(X), dimnames = list(NULL, names(X)))
  rownames(M) <- feat
  common <- intersect(colnames(M), rownames(prep$meta))
  if (length(common) < 4)
    stop(sprintf("Only %d sample IDs in layer '%s' match the metadata. Column names must equal metadata sample IDs.",
                 length(common), name))
  M <- M[, common, drop = FALSE]; meta <- prep$meta[common, , drop = FALSE]

  if (!id_map %in% c("none", "ensembl_hs", "entrez_hs")) stop("Unknown legacy ID conversion.")
  if (id_map != "none") {
    if (!identical(organism, "Human")) stop("Legacy human ID conversion requires an explicit Human organism declaration. Use layer annotation for other species.")
    annotation <- utils::modifyList(annotation %||% list(), list(mode = "human",
      id_type = if (id_map == "ensembl_hs") "ensembl" else "entrez"))
  }
  idinfo <- byod_detect_ids(rownames(M), omics)
  idinfo$coverage_before <- byod_id_coverage(rownames(M))
  idinfo$coverage_after <- idinfo$coverage_before
  idinfo$map <- list(applied = FALSE, msg = "Annotation is separate from the expression matrix.")

  diag0 <- byod_diagnose(M)
  cl <- byod_clean_matrix(M, clean, is_counts = identical(diag0$scale, "raw counts"))
  M <- cl$M
  if (nrow(M) < 10)
    stop(sprintf("Layer '%s': cleaning left only %d features. Relax the filters: %s",
                 name, nrow(M), paste(byod_clean_summary(cl), collapse = "; ")))
  if (ncol(M) < 4)
    stop(sprintf("Layer '%s': cleaning left only %d profiles.", name, ncol(M)))
  if (isTRUE(cl$any_applied)) {
    idinfo$ids_after_clean <- byod_detect_ids(rownames(M), omics)
    idinfo$coverage_after <- byod_id_coverage(rownames(M))
    idinfo$label <- idinfo$ids_after_clean$label
  }

  if (length(cl$dropped_samples)) {
    meta <- meta[setdiff(rownames(meta), cl$dropped_samples), , drop = FALSE]
    meta <- meta[colnames(M), , drop = FALSE]
  }
  raw_libsize <- colSums(M, na.rm = TRUE)                # input-scale check + optional normalization
  diagn <- byod_diagnose(M)
  qc <- NULL
  if (identical(diagn$scale, "raw counts") && ncol(M) > 3) {
    med <- stats::median(raw_libsize)
    if (is.finite(med) && med > 0) {
      ratio <- raw_libsize / med
      dead <- names(ratio)[is.finite(ratio) & ratio < 0.01]
      if (length(dead))
        qc <- list(kind = "library_size", samples = dead, ratio = ratio[dead], median = med,
                   msg = sprintf("%d profile%s a total count below 1%% of the median (%s): %s. %s deposited that way and %s been removed. This is usually a failed library, but in some designs it is real - a pathogen matrix before inoculation genuinely has almost no reads. Decide which applies here before drawing conclusions from %s.",
                                 length(dead), if (length(dead) == 1) " has" else "s have",
                                 format(round(med), big.mark = ","),
                                 paste(sprintf("%s = %s", dead,
                                               format(round(raw_libsize[dead]), big.mark = ",")),
                                       collapse = ", "),
                                 if (length(dead) == 1) "It was" else "They were",
                                 if (length(dead) == 1) "has not" else "have not",
                                 if (length(dead) == 1) "it" else "them"))
    }
  }
  nz <- byod_normalize(M, norm, diagn); M <- nz$M
  M <- M[byod_rowvars(M) > 0, , drop = FALSE]
  if (nrow(M) < 10) stop(sprintf("Layer '%s': fewer than 10 variable features after filtering.", name))
  batch_msg <- "none"                                    # optional ComBat batch correction
  if (!is.null(batch_col) && batch_col %in% names(meta)) {
    cb <- byod_combat(M, meta[[batch_col]], protect = meta$Timepoint); M <- cb$M
    batch_msg <- if (cb$applied) sprintf("%s — %s", batch_col, cb$msg) else sprintf("%s — skipped (%s)", batch_col, cb$msg)
    M <- M[byod_rowvars(M) > 0, , drop = FALSE]
  }
  Mi <- M
  if (anyNA(Mi)) { rm_ <- rowMeans(Mi, na.rm = TRUE); ix <- which(is.na(Mi), arr.ind = TRUE); Mi[ix] <- rm_[ix[, 1]] }
  rv2  <- byod_rowvars(Mi); keep <- head(order(rv2, decreasing = TRUE), min(max_genes_pca, sum(rv2 > 0)))
  pc   <- prcomp(t(Mi[keep, , drop = FALSE]), scale. = TRUE)
  npc  <- min(20L, ncol(pc$x)); pca <- pc$x[, seq_len(npc), drop = FALSE]
  varexp <- round(100 * (pc$sdev^2 / sum(pc$sdev^2))[seq_len(npc)], 1)
  umap_seed <- trace_seed()
  umap <- trace_compute_umap(pca, umap_seed)
  nlev <- function(v) length(unique(meta[[v]][!is.na(meta[[v]])]))
  cc <- prep$clinical_cat; if (!is.null(cc)) cc <- cc[vapply(cc, function(v) nlev(v) >= 2, logical(1))]
  cn <- prep$clinical_num; if (!is.null(cn)) cn <- cn[vapply(cn, function(v) nlev(v) >= 2, logical(1))]
  d <- list(expr = M, meta = meta, pca = pca, varexp = varexp, umap = umap, umap_seed = umap_seed,
       features = rownames(M), species = name, omic = name, tp_levels = prep$tp_levels, unit = unit,
       organism = if (is.null(organism) || !nzchar(organism)) NULL else organism,
       ids = idinfo, qc = qc, clean = cl[setdiff(names(cl), "M")], omics = omics,
       input_expression = expr_df,
       subject = prep$subject, diag = diagn, norm = nz$applied, batch = batch_msg, libsize = raw_libsize,
       clinical_cat = if (length(cc)) cc else NULL, clinical_num = if (length(cn)) cn else NULL,

       subject_varying = prep$subject_varying)
  d$annotation <- annotation_build(rownames(M), organism, annotation)
  ac <- annotation_collections(d)
  d$genesets_apply <- any(vapply(ac, function(k) byod_genesets_apply(d, k), logical(1)))
  hits <- intersect(rownames(M), unique(unlist(annotation_sets(d, "Hallmark"))))
  d$ids$coverage_after <- list(hit = length(hits), n = nrow(M))
  d
}

build_byod_study <- function(exprs, meta_df, tp_col = NULL, units = NULL, name = "Study",
                             tp_order = NULL, subject_col = NULL, norms = NULL, batch_col = NULL,
                             id_maps = NULL, cleans = NULL, omics = NULL, nest_col = NULL,
                             organism = NULL, annotations = NULL) {
  organism <- annotation_organism(organism)
  if (!length(exprs)) stop("No expression layers provided.")
  prep <- byod_prep_meta(meta_df, tp_col, tp_order, subject_col)
  layers <- list()

  if (!is.null(nest_col) && nzchar(nest_col)) {
    if (!nest_col %in% names(prep$meta))
      stop(sprintf("Higher-level unit '%s' is not a column in the metadata.", nest_col))
    if (is.null(subject_col))
      stop("A higher-level unit needs a subject column: declare which identifier nests inside it.")
    if (!byod_nest_ok(prep$meta, subject_col, nest_col))
      stop(sprintf(paste("'%s' does not nest '%s': every %s must belong to exactly one %s,",
                         "and there must be fewer %s than %s. Leave it blank to analyse at",
                         "subject level."),
                   nest_col, subject_col, subject_col, nest_col, nest_col, subject_col))
  }
  for (ln in names(exprs)) {
    u  <- if (!is.null(units) && !is.null(units[[ln]])) units[[ln]] else "expression"
    nm <- if (!is.null(norms) && !is.null(norms[[ln]])) norms[[ln]] else "Auto-detect"
    im <- if (!is.null(id_maps) && !is.null(id_maps[[ln]])) id_maps[[ln]] else "none"
    cle <- if (!is.null(cleans) && !is.null(cleans[[ln]])) cleans[[ln]] else NULL
    om <- if (!is.null(omics) && !is.null(omics[[ln]])) omics[[ln]] else "other"
    layers[[ln]] <- byod_build_layer(exprs[[ln]], prep, unit = u, name = paste0(name, " · ", ln),
                                     norm = nm, batch_col = batch_col, id_map = im,
                                     clean = cle, omics = om, organism = organism,
                                     annotation = annotations[[ln]])

    layers[[ln]]$nest <- if (!is.null(nest_col) && nzchar(nest_col)) nest_col else NULL
  }
  shared_s <- if (length(layers) >= 2) Reduce(intersect, lapply(layers, function(d) colnames(d$expr))) else character(0)
  shared_g <- if (length(layers) >= 2) Reduce(intersect, lapply(layers, function(d) rownames(d$expr))) else character(0)
  list(name = name, meta = prep$meta, input_metadata = meta_df, tp_levels = prep$tp_levels, subject = subject_col,
       organism = if (is.null(organism) || !nzchar(organism)) NULL else organism,
       nest = if (!is.null(nest_col) && nzchar(nest_col)) nest_col else NULL,
       excluded_untimed = prep$excluded_untimed,
       clinical_cat = prep$clinical_cat, clinical_num = prep$clinical_num,
       subject_varying = prep$subject_varying,
       layer_names = names(layers), layers = layers,
       matched_samples = shared_s, shared_genes = shared_g)
}

builtin_wrap_layer <- function(d, subject_col = "Kidney_ID", omic_name = NULL, nest_col = NULL) {
  d$subject <- if (subject_col %in% names(d$meta)) subject_col else NULL

  d$nest <- if (!is.null(nest_col) && byod_nest_ok(d$meta, d$subject, nest_col)) nest_col else NULL
  if (is.null(d$diag))    d$diag    <- byod_diagnose(d$expr)
  if (is.null(d$norm))    d$norm    <- "pre-computed (see Methods)"
  if (is.null(d$batch))   d$batch   <- "none (curated pipeline)"
  if (is.null(d$libsize)) d$libsize <- colSums(d$expr, na.rm = TRUE)
  if (!is.null(omic_name)) d$omic <- omic_name
  if (is.null(d$clinical_cat)) d$clinical_cat <- NULL
  if (is.null(d$clinical_num)) d$clinical_num <- NULL
  d
}

build_builtin_study <- function(mrna_key, prot_key, name, tp_levels,
                                subject_col = "Kidney_ID", nest_col = NULL,
                                species = c("Human", "Pig")) {
  species <- match.arg(species)
  L <- list()
  L[["mRNA"]]    <- builtin_wrap_layer(DS[[mrna_key]], subject_col, "mRNA", nest_col)
  L[["Protein"]] <- builtin_wrap_layer(DS[[prot_key]], subject_col, "Protein", nest_col)

  for (.ln in names(L)) {
    L[[.ln]]$organism <- species
    L[[.ln]]$annotation <- annotation_build(rownames(L[[.ln]]$expr), species,
      list(id_species = "human", mode = "human", source = "curated atlas; upstream human-symbol annotation"))
    L[[.ln]]$atlas_annotation <- TRUE
    if (species == "Human") {
      # Explicit contrasts used in the reference cohort's clinical pathway figure.
      L[[.ln]]$clinical_contrasts <- list(
        "CPR performed" = c("No", "Yes"), "Hypertension" = c("No", "Yes"),
        "Diabetes" = c("Yes", "No"), "Donor sex" = c("Female", "Male"),
        "Donor Type" = c("DCD", "DBD"), "Preservation technique" = c("HMP", "SCS"))
    }
  }

  m1 <- L[["mRNA"]]$meta; m2 <- L[["Protein"]]$meta
  common_cols <- intersect(names(m1), names(m2))

  extra <- setdiff(rownames(m2), rownames(m1))
  meta  <- rbind(m1[, common_cols, drop = FALSE], m2[extra, common_cols, drop = FALSE])
  meta$Timepoint <- factor(as.character(meta$Timepoint), levels = tp_levels)
  shared_s <- intersect(colnames(L[["mRNA"]]$expr), colnames(L[["Protein"]]$expr))
  shared_g <- intersect(rownames(L[["mRNA"]]$expr), rownames(L[["Protein"]]$expr))
  cc <- L[["mRNA"]]$clinical_cat %||% L[["Protein"]]$clinical_cat
  cn <- L[["mRNA"]]$clinical_num %||% L[["Protein"]]$clinical_num

  list(name = name, meta = meta, tp_levels = tp_levels, subject = subject_col, organism = species,
       nest = nest_col,
       clinical_cat = cc, clinical_num = cn,
       layer_names = names(L), layers = L,
       matched_samples = shared_s, shared_genes = shared_g,
       builtin = TRUE, species = species)
}

STUDIES <- list(

  nmp_pig = build_builtin_study("pig_mrna", "pig_protein",
              "Kidney Perfusion Reference — porcine", atlas$tp_pig, "Kidney_ID",
              nest_col = "Pig_Number", species = "Pig"),
  nmp_hum = build_builtin_study("hum_mrna", "hum_protein",
              "Kidney Perfusion Reference — human", atlas$tp_hum, "Kidney_ID",
              species = "Human")
)
for (.study_key in names(STUDIES))
  STUDIES[[.study_key]]$dataset_id <- paste0("builtin-", .study_key)
rm(.study_key)

byod_human_example <- function(dual = TRUE) {
  s <- STUDIES$nmp_hum
  s$reference_key <- "nmp_hum"
  s$dataset_id <- NULL
  s$builtin <- FALSE
  s$name <- if (dual) "Human kidney reference (mRNA + protein)" else "Human kidney reference (mRNA)"
  if (!dual) {
    s$layers <- s$layers["mRNA"]
    s$meta <- s$layers$mRNA$meta
    s$matched_samples <- s$shared_genes <- character(0)
  }
  s$layer_names <- names(s$layers)
  s$input_metadata <- data.frame(sample = rownames(s$meta), s$meta, check.names = FALSE)
  for (ln in s$layer_names) {
    d <- s$layers[[ln]]
    # Reuse the curated values and embeddings without another normalization pass.
    d$input_expression <- data.frame(feature = rownames(d$expr), d$expr, check.names = FALSE)
    d$genesets_apply <- any(vapply(annotation_collections(d), function(k) byod_genesets_apply(d, k), logical(1)))
    s$layers[[ln]] <- d
  }
  s
}

STUDY_LABELS <- c(
  nmp_pig = "Kidney Perfusion Reference — porcine (mRNA + protein, in vivo → 6 h)",
  nmp_hum = "Kidney Perfusion Reference — human (mRNA + protein + clinical, preNMP → 6 h)"
)

MODULE_HELP <- list(
  list(id="overview",   title="Overview & Normalization QC",
       body="Summary of the active study and per-sample value distributions for the active layer. Use the QC boxplots to confirm layers are on a comparable scale before comparing them; if a layer looks unnormalised, re-load with a different normalization method (uploads) — built-in studies are pre-normalised (mRNA = log2(TPM+1); protein = log2 intensity)."),
  list(id="landscape",  title="Landscape (PCA / UMAP)",
       body="Sample embedding coloured by any metadata column or by a gene's expression. UMAP uses the active random seed (default 1234), n_neighbors=15 (capped to n-1), min_dist=0.30 and single-threaded optimisation. It updates when a new seed is applied. PCA is deterministic. Use PCA to judge global structure/outliers; UMAP for local grouping."),
  list(id="trajectory", title="Trajectory (feature over time)",
       body="A single feature's trajectory across the time/condition axis. 'Mean ± SEM' shows group means; 'GAM smooth + band' fits a smooth trend (mgcv) with 95% band — use GAM when you have ≥4 timepoints. Split lines by any categorical variable (e.g. warm-ischemia status)."),
  list(id="differential", title="Differential expression (limma)",
       body="Moderated t-test between two groups (limma). Defaults: FDR=0.05, |log2FC|=1. 'Paired (subject block)' fits subject as a blocking factor — keep ON for repeated-measures designs (the NMP studies are paired by kidney). Label top N by combined significance. Live re-analysis; may differ slightly from the companion paper's full pipeline."),
  list(id="enrichment", title="Enrichment (GSEA / ORA)",
       body="GSEA ranks ALL genes by the DE t-statistic (fgsea; minSize=10, maxSize=500) — no arbitrary cutoff. ORA is a hypergeometric test on a gene list you paste, against the layer's measured features as the universe. Collections: Hallmark (recommended first pass), KEGG, Reactome, GO:BP."),
  list(id="integration", title="mRNA–Protein integration",
       body="For a chosen contrast, computes log2 fold-changes in both layers on genes measured in both, then correlates them (Pearson r). Points in the upper-right / lower-left are concordant (transcript and protein move together); off-diagonal points are discordant (post-transcriptional regulation)."),
  list(id="concordance", title="Concordance (genome-wide)",
       body="Per-gene correlation of mRNA and protein across matched samples (Spearman by default). The histogram shows the genome-wide distribution; positive tail = concordant genes, near-zero/negative = discordant. Click a gene to see its paired scatter. Requires two layers on matched samples."),
  list(id="pathways",   title="Pathways (GSVA / ssGSEA activity)",
       body="Sample-level pathway activity scores (GSVA or ssGSEA) for a collection, shown as a heatmap or as a single pathway's trajectory over time. Built-in studies use precomputed Hallmark GSVA when the GSVA package is unavailable. Good for asking 'when does this programme switch on during perfusion?'."),
  list(id="wgcna",      title="Co-expression (WGCNA)",
       body="Weighted gene co-expression network → modules → module–trait correlation heatmap. Defaults: top 2500 genes by variance, soft-power auto-picked, minModuleSize=20, deepSplit=2. Needs ≥15–20 samples to be meaningful. Use it to find gene modules tracking time or a clinical trait."),
  list(id="diablo",     title="DIABLO (multi-omics discrimination)",
       body="mixOmics DIABLO: a supervised multi-block model that finds correlated mRNA+protein signatures discriminating an outcome (e.g. timepoint or warm-ischemia). Defaults: nfeat≈25 per block, ≥8 matched samples required. Reports block loadings and the multi-omics sample plot. The FIT is descriptive and uses every matched sample, so its apparent separation is not evidence of generalisation — press 'Cross-validate' for that. The cross-validation holds out whole independent units (the declared nesting level, else the subject) and recomputes the variance filter, the scaling and the feature selection inside each training fold, so a sibling sample of a held-out animal never remains in training."),
  list(id="association", title="Metadata / clinical association",
       body="Associates a feature (gene or pathway) with a metadata variable, optionally restricted to one timepoint. Categorical variables → grouped boxplots (Kruskal–Wallis); numeric variables → scatter with Pearson r. For the human NMP study this exposes 33 donor characteristics."),
  list(id="crossspecies", title="Cross-species (pig vs human)",
       body="Shows a gene's perfusion trajectory in the porcine and human reference cohorts. The title identifies the human-pig mapping as one-to-one or ambiguous using the bundled, versioned orthology table. This is a descriptive comparison between cohorts."),
  list(id="clinical", title="Subject variables: features (one value per subject)",
       body="Genome-wide association of every feature with a donor/clinical variable, computed at DONOR resolution (one value per subject) so the correlation is across donors (n≤24), never across the ~5 timepoint biopsies per donor (which would pseudoreplicate). Continuous traits stay continuous (Spearman/Pearson); 2-level categoricals use a donor-level moderated-t contrast. BH FDR across features. Volcano + single-feature drill-down + enrichment of the top-correlated features. Human-study only. Low-coverage traits (effective donor n below the minimum) are declined rather than reported."),
  list(id="grouptime", title="Group×time interaction",
       body="Tests whether two groups have different time courses using limma. For the porcine paired design, animal-by-time and kidney terms account for the sampling structure. Where this fit has too few residual degrees of freedom, the module uses a subject-adjusted model. The result reports the fitted design and number of units. These are fixed-block screening models, not a fitted nested random-effects model."),
  list(id="temporal", title="Temporal programs",
       body="Genome-wide temporal-trend detection (limma over timepoints, subject-adjusted) followed by fuzzy c-means clustering of the FDR-significant features into temporal programs. Shows per-program mean profiles, a pattern heatmap, and per-program enrichment. Fixed seed for reproducible clusters. Needs ≥3 timepoints.")
)
module_help_ui <- function() {
  do.call(bslib::accordion, c(list(open = FALSE, class = "mb-2"),
    lapply(MODULE_HELP, function(h)
      bslib::accordion_panel(h$title, icon = shiny::icon("circle-question"),
        shiny::p(class = "statlab", style = "font-size:.82rem;", h$body)))))
}

# The effect is g2 minus g1. Resolve pairing at the declared independent unit.
byod_run_limma <- function(d, groupvar, g1, g2, subject_col = NULL, paired = TRUE, at_tp = "All",
                           filter_constant = FALSE) {
  gf <- grp_factor(d, groupvar)
  sel <- which(as.character(gf) %in% c(g1, g2) & !is.na(gf))
  if (!is.null(at_tp) && at_tp != "All") sel <- sel[d$meta$Timepoint[sel] == at_tp]
  empty <- data.frame(Feature = character(), logFC = numeric(), AveExpr = numeric(),
                      t = numeric(), P.Value = numeric(), adj.P.Val = numeric(), stringsAsFactors = FALSE)
  constant_features <- character()
  filter_record <- function(out) {
    attr(out, "filter_constant") <- isTRUE(filter_constant)
    attr(out, "constant_features") <- constant_features
    attr(out, "removed_features") <- if (isTRUE(filter_constant)) constant_features else character()
    out
  }
  declined <- function(msg) { out <- empty; attr(out, "declined") <- msg; filter_record(out) }
  if (length(sel) < 4)
    return(declined(sprintf("Only %d profile%s in '%s' and '%s' combined; a contrast needs at least 4.",
                            length(sel), if (length(sel) == 1) "" else "s", g1, g2)))
  mt <- d$meta[sel, , drop = FALSE]; E <- d$expr[, sel, drop = FALSE]
  grp <- factor(ifelse(as.character(gf[sel]) == g2, "g2", "g1"), levels = c("g1", "g2"))
  if (nlevels(droplevels(grp)) < 2 || min(table(grp)) < 2) {
    tb <- table(factor(ifelse(as.character(gf[sel]) == g2, g2, g1), levels = c(g1, g2)))
    return(declined(sprintf("Each side needs at least 2 profiles; this contrast has %s.",
                            paste(names(tb), as.integer(tb), sep = " = ", collapse = ", "))))
  }
  subj0 <- if (!is.null(subject_col) && subject_col %in% names(mt)) as.character(mt[[subject_col]]) else NULL
  if (isTRUE(paired) && (is.null(subj0) || anyNA(subj0) || any(!nzchar(subj0))))
    return(declined("Paired analysis requires a declared subject ID for every selected sample. Declare the pairing or turn Paired off."))
  paired_note <- FALSE; block_col <- NA_character_; collapsed <- FALSE
  excluded_subjects <- character(); excluded_profiles <- 0L
  constant_features <- rownames(E)[matrixStats::rowVars(E) == 0]
  if (isTRUE(filter_constant)) {
    E <- E[!rownames(E) %in% constant_features, , drop = FALSE]
    if (nrow(E) < 2L)
      return(declined("Fewer than two features remain after excluding features constant across the selected profiles."))
  }

  n_units <- NA_integer_; unit_lab <- NA_character_

  between_subj <- !is.null(subj0) && all(tapply(as.character(grp), subj0, function(z) length(unique(z))) == 1)
  if (between_subj) {
    collapsed <- TRUE
    su <- unique(subj0); idx <- lapply(su, function(u) which(subj0 == u))
    Ec <- vapply(idx, function(ii) rowMeans(E[, ii, drop = FALSE], na.rm = TRUE), numeric(nrow(E)))
    dimnames(Ec) <- list(rownames(E), su)
    grpc <- factor(vapply(idx, function(ii) as.character(grp[ii[1]]), character(1)), levels = c("g1", "g2"))
    metac <- mt[vapply(idx, `[`, integer(1), 1), , drop = FALSE]

    if (min(table(grpc)) < MIN_UNITS$group_contrast) {
      return(declined(byod_unit_floor_msg(table(grpc), MIN_UNITS$group_contrast,
                                          "A between-subject group contrast", ncol(E))))
    }
    if (isTRUE(paired)) {
      block_col <- d$nest %||% NULL
      if (is.null(block_col) || !byod_nest_ok(mt, subject_col, block_col))
        return(declined(paste("Each subject belongs to only one group. Paired analysis needs an explicitly declared higher-level pairing unit;",
                              "no metadata column is selected automatically. For independent groups, turn Paired off.")))
      blocks <- as.character(metac[[block_col]])
      complete <- names(which(vapply(split(as.character(grpc), blocks),
        function(z) all(c("g1", "g2") %in% z), logical(1))))
      if (length(complete) < MIN_UNITS$group_contrast)
        return(declined(sprintf("Only %d complete pairs in declared unit '%s'; at least %d are required.",
                                length(complete), block_col, MIN_UNITS$group_contrast)))
      excluded_subjects <- su[!blocks %in% complete]
      excluded_profiles <- sum(subj0 %in% excluded_subjects)
      # Each declared block contributes one group difference, regardless of its size.
      Ef <- vapply(complete, function(b) rowMeans(Ec[, blocks == b & grpc == "g2", drop = FALSE]) -
                       rowMeans(Ec[, blocks == b & grpc == "g1", drop = FALSE]), numeric(nrow(Ec)))
      dimnames(Ef) <- list(rownames(Ec), complete)
      design <- matrix(1, ncol(Ef), 1, dimnames = list(complete, "(Intercept)"))
      paired_note <- TRUE; n_units <- length(complete); unit_lab <- block_col
    } else { Ef <- Ec; grpf <- grpc; design <- model.matrix(~ grpf)
      n_units <- ncol(Ec); unit_lab <- subject_col %||% "subject" }  # collapsed table: one column per subject
  } else {

    use_pair <- isTRUE(paired) && !is.null(subj0) && !any(is.na(subj0)); Ef <- E; grpf <- grp; subjf <- NULL
    nest_col <- d$nest %||% NULL
    if (use_pair && !is.null(nest_col) && byod_nest_ok(mt, subject_col, nest_col)) {
      nf <- byod_nested_paired_fit(E, subj0, grp, as.character(mt[[nest_col]]))
      if (!is.null(nf) && !is.null(nf$declined)) return(declined(nf$declined))
      if (!is.null(nf)) {
        tt <- nf$tab
        out <- data.frame(Feature = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr,
                          t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
                          row.names = NULL, stringsAsFactors = FALSE)
        attr(out, "paired") <- TRUE; attr(out, "n") <- nf$n_profiles
        attr(out, "n_units") <- as.integer(nf$n_units); attr(out, "unit") <- nest_col
        attr(out, "n_subjects") <- as.integer(nf$n_subjects); attr(out, "subject_level") <- subject_col
        attr(out, "block") <- nest_col; attr(out, "collapsed") <- TRUE
        attr(out, "nested") <- TRUE
        attr(out, "excluded_subjects") <- nf$excluded_subjects
        attr(out, "excluded_profiles") <- ncol(E) - nf$n_profiles
        return(filter_record(out))
      }
    }
    if (use_pair) {
      sj <- factor(subj0); tab <- table(sj, grp)
      keepS <- rownames(tab)[tab[, "g1"] > 0 & tab[, "g2"] > 0]; keep <- sj %in% keepS
      if (length(keepS) >= 2 && sum(keep) >= 4) {
        Ef <- E[, keep, drop = FALSE]; grpf <- droplevels(grp[keep]); subjf <- droplevels(sj[keep])
        use_pair <- nlevels(subjf) >= 2 && nlevels(grpf) == 2 && ncol(Ef) > nlevels(subjf) + 1
        excluded_subjects <- setdiff(unique(subj0), keepS)
        excluded_profiles <- sum(!keep)
      } else use_pair <- FALSE
      if (!use_pair) return(declined("Too few complete subject pairs or residual degrees of freedom. No unpaired fallback was fitted."))
    }
    if (!use_pair) { Ef <- E; grpf <- grp }
    design <- if (use_pair) model.matrix(~ subjf + grpf) else model.matrix(~ grpf); paired_note <- use_pair
    if (use_pair) { n_units <- nlevels(subjf); unit_lab <- subject_col %||% "subject" }
    else { n_units <- ncol(Ef); unit_lab <- "profile" }   # no subject term: the profile is the unit

    if (n_units < MIN_UNITS$group_contrast)
      return(declined(byod_unit_floor_msg(
        stats::setNames(n_units, unit_lab %||% "unit"), MIN_UNITS$group_contrast,
        "A differential-expression contrast", ncol(Ef))))
  }
  fit <- tryCatch(eBayes(lmFit(Ef, design)), error = function(e) NULL)
  if (is.null(fit)) return(declined(sprintf(
    "The model could not be fitted: %d profiles against %d design columns leaves no residual degrees of freedom.",
    ncol(Ef), ncol(design))))
  coefn <- if (ncol(design) == 1L) 1L else grep("^grpf", colnames(design), value = TRUE)[1]
  tt <- topTable(fit, coef = coefn, number = Inf, sort.by = "P")

  out <- data.frame(Feature = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr,
             t = tt$t, P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val,
             row.names = NULL, stringsAsFactors = FALSE)
  attr(out, "paired") <- isTRUE(paired_note); attr(out, "n") <- ncol(E) - excluded_profiles
  attr(out, "n_units") <- as.integer(n_units); attr(out, "unit") <- unit_lab
  attr(out, "block") <- block_col; attr(out, "collapsed") <- collapsed
  attr(out, "excluded_subjects") <- excluded_subjects
  attr(out, "excluded_profiles") <- excluded_profiles
  filter_record(out)
}

byod_dimred <- function(d, type = "UMAP", colorby = "Timepoint", gene = NULL, cols = NULL, scale_stops = NULL) {
  co <- if (type == "UMAP") trace_umap(d) else as.data.frame(d$pca[, 1:2])
  ax <- if (type == "UMAP") c("UMAP1", "UMAP2")
        else c(paste0("PC1 (", d$varexp[1], "%)"), paste0("PC2 (", d$varexp[2], "%)"))
  df <- data.frame(sample = rownames(co), x = co[, 1], y = co[, 2],
                   d$meta[rownames(co), , drop = FALSE], check.names = FALSE, stringsAsFactors = FALSE)
  ttl <- paste0(d$species, " — ", type)
  grad <- identical(colorby, "Gene expression") ||
          (!is.null(colorby) && colorby %in% names(d$meta) && is_num_var(d, colorby))
  if (grad) {
    if (identical(colorby, "Gene expression")) {
      if (!has_feature(d, gene)) return(list(plot = emptyplot(paste(gene, "not detected")), summary = data.frame(),
        kind = "gradient", levels = character(0), ttl = ttl, ax = ax, scale_title = gene %||% "value"))
      cv <- as.numeric(d$expr[gene, rownames(co)]); sttl <- gene
    } else { cv <- as.numeric(df[[colorby]]); sttl <- colorby }
    stops <- scale_stops %||% heat_stops("YlOrRd")
    p <- plot_ly(df, x = ~x, y = ~y, type = "scatter", mode = "markers",
      marker = list(color = cv, colorscale = heat_colorscale(stops), showscale = TRUE,
        colorbar = list(title = sttl, len = .6), size = 9, line = list(width = .4, color = "white")),
      text = ~paste0(sample, "<br>", Timepoint, "<br>", sttl, " = ", signif(cv, 3)), hoverinfo = "text") %>%
      ply(xaxis = list(title = ax[1]), yaxis = list(title = ax[2]), title = list(text = ttl, x = .02, font = list(size = 15)))
    return(list(plot = p, summary = data.frame(sample = df$sample, x = df$x, y = df$y, value = cv, stringsAsFactors = FALSE),
      kind = "gradient", levels = character(0), ttl = ttl, ax = ax, scale_title = sttl))
  }
  g <- droplevels(grp_factor(d, colorby)[match(rownames(co), rownames(d$meta))])
  df$grp <- g; df <- df[!is.na(df$grp), , drop = FALSE]; g <- droplevels(df$grp); lv <- levels(g)
  cc <- if (is.null(cols)) setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = colorby)), lv) else cols[lv]
  p <- plot_ly(df, x = ~x, y = ~y, color = g, colors = unname(cc[lv]),
      type = "scatter", mode = "markers", marker = list(size = 9, line = list(width = .4, color = "white")),
      text = ~paste0(sample, "<br>", Timepoint), hoverinfo = "text") %>%
      ply(xaxis = list(title = ax[1]), yaxis = list(title = ax[2]), title = list(text = ttl, x = .02, font = list(size = 15)))
  list(plot = p, summary = data.frame(sample = df$sample, x = df$x, y = df$y, group = as.character(g), stringsAsFactors = FALSE),
       kind = "discrete", levels = lv, ttl = ttl, ax = ax, scale_title = colorby)
}

byod_empty_tc <- function(d, ttl, ylab) list(plot = emptyplot(ttl),
  summary = data.frame(Timepoint = factor(levels = d$tp_levels), Group = character(),
                       mean = numeric(), se = numeric(), n = numeric()),
  splits = character(0), ylab = ylab, ttl = ttl, cols = setNames(character(0), character(0)))
byod_gene_tc <- function(d, gene, splitby = "None", cols = NULL) {
  ttl <- paste0(gene, " over Timepoint"); v <- expr_vec(d, gene)
  if (is.null(v)) return(byod_empty_tc(d, paste(gene, "not found"), d$unit))
  td <- tc_data(d, setNames(v, colnames(d$expr)), splitby)
  cc <- if (is.null(cols)) resolve_cols(td$splits, PAL_NAMES[1], if (identical(splitby, "None")) NULL else splitby) else cols
  list(plot = tc_plotly(td, d$unit, ttl, cc), summary = td$summary, splits = td$splits,
       ylab = d$unit, ttl = ttl, cols = cc)
}

byod_heatmap <- function(d, genes, by = "Timepoint", colorscale = NULL) {
  genes <- intersect(genes, rownames(d$expr))
  if (length(genes) < 2) return(list(plot = emptyplot("Enter ≥ 2 valid genes"), mat = NULL, ttl = "", by = by))
  gf <- droplevels(grp_factor(d, by)); lv <- levels(gf)
  if (!length(lv)) return(list(plot = emptyplot("No groups available for this variable"), mat = NULL, ttl = "", by = by))
  M <- sapply(lv, function(L) rowMeans(d$expr[genes, which(gf == L), drop = FALSE], na.rm = TRUE))
  if (is.null(dim(M))) M <- matrix(M, nrow = length(genes), dimnames = list(genes, lv))
  Z <- t(scale(t(M))); Z[!is.finite(Z)] <- 0
  if (nrow(Z) > 2) Z <- Z[hclust(dist(Z))$order, , drop = FALSE]
  cs <- colorscale %||% list(c(0, HUMAN), c(.5, "#FFFFFF"), c(1, PIG))
  ttl <- paste(length(genes), "genes (row z-score)")
  p <- plot_ly(x = colnames(Z), y = rownames(Z), z = Z, type = "heatmap", colorscale = cs, zmid = 0,
    hovertemplate = "%{y} @ %{x}<br>z=%{z:.2f}<extra></extra>") %>%
    ply(xaxis = list(title = by), yaxis = list(title = ""), title = list(text = ttl, x = .02, font = list(size = 14)))
  list(plot = p, mat = Z, ttl = ttl, by = by)
}

byod_enrich_gsea <- function(d, groupvar, g1, g2, subject_col, paired, collection, n = 20) {
  tt <- byod_run_limma(d, groupvar, g1, g2, subject_col, paired)
  if (!nrow(tt)) {

    why <- attr(tt, "declined") %||% "Too few samples for this contrast"
    res <- data.frame(); attr(res, "declined") <- why
    return(list(plot = emptyplot(why), res = res))
  }
  stats <- setNames(tt$t, tt$Feature)
  res <- run_gsea(stats, collection, d = d)
  list(plot = plot_gsea(res, n, paste0(collection, ": ", g2, " vs ", g1))$plot, res = res)
}
byod_enrich_ora <- function(d, genes, collection, n = 20) {
  uni <- rownames(d$expr); res <- run_ora(intersect(genes, uni), collection, uni, d = d)
  list(plot = plot_ora(res, n, paste0(collection, " ORA (", length(intersect(genes, uni)), " genes)"))$plot, res = res)
}

byod_integrate <- function(d1, d2, groupvar, g1, g2, subject_col, paired,
                           lab1 = "layer 1", lab2 = "layer 2", cols = NULL, topn = 300) {
  QUAD <- c("both up", "both down", "discordant")
  fail <- function(m) list(plot = emptyplot(m), tab = data.frame(), r = NA_real_, p = NA_real_, n = 0,
    summary = data.frame(), splits = character(0), ylab = paste0(lab2, " mean difference"), ttl = m,
    cols = setNames(character(0), character(0)))
  rt <- byod_run_limma(d1, groupvar, g1, g2, subject_col, paired)
  pt <- byod_run_limma(d2, groupvar, g1, g2, subject_col, paired)
  if (!nrow(rt) || !nrow(pt)) return(fail("Too few samples in one layer for this contrast"))
  m <- merge(rt[, c("Feature", "logFC", "adj.P.Val")], pt[, c("Feature", "logFC", "adj.P.Val")],
             by = "Feature", suffixes = c("_1", "_2"))
  if (nrow(m) < 3) return(fail(sprintf("Only %d shared features between the two layers", nrow(m))))
  effect_r <- suppressWarnings(cor(m$logFC_1, m$logFC_2, use = "pairwise.complete.obs"))
  m$quad <- factor(ifelse(m$logFC_1 > 0 & m$logFC_2 > 0, "both up",
             ifelse(m$logFC_1 < 0 & m$logFC_2 < 0, "both down", "discordant")), levels = QUAD)
  ms <- head(m[order(-abs(m$logFC_1 + m$logFC_2)), ], topn)
  splits <- intersect(QUAD, as.character(unique(ms$quad)))
  cc <- if (is.null(cols)) fb_cols(splits, NULL, PAL_NAMES[2]) else cols
  ttl <- sprintf("%s vs %s · descriptive r = %.2f (%d shared features)", g2, g1,
                 effect_r, nrow(m))
  p <- plot_ly(ms, x = ~logFC_1, y = ~logFC_2, color = ~factor(quad, splits), colors = unname(cc[splits]),
      type = "scatter", mode = "markers", marker = list(size = 7, opacity = .75),
      text = ~paste0(Feature, "<br>", lab1, " logFC=", logFC_1, "<br>", lab2, " logFC=", logFC_2), hoverinfo = "text") %>%
    ply(xaxis = list(title = paste0(lab1, " mean difference"), zeroline = TRUE),
        yaxis = list(title = paste0(lab2, " mean difference"), zeroline = TRUE),
        legend = list(orientation = "h", y = 1.1), title = list(text = ttl, x = .02, font = list(size = 13)))
  list(plot = p, tab = m[order(m$adj.P.Val_1), ], r = effect_r, p = NA_real_, n = nrow(m),
       summary = ms, splits = splits, xlab = paste0(lab1, " mean difference"), ylab = paste0(lab2, " mean difference"),
       ttl = ttl, cols = cc)
}

byod_mp_plot <- function(o, cols = NULL) {
  ms <- o$summary; if (is.null(ms) || !nrow(ms)) return(emptyplot(gsub("<[^>]+>", "", o$ttl %||% "no data")))
  splits <- o$splits; cc <- if (is.null(cols)) o$cols else cols
  plot_ly(ms, x = ~logFC_1, y = ~logFC_2, color = ~factor(quad, splits), colors = unname(cc[splits]),
      type = "scatter", mode = "markers", marker = list(size = 7, opacity = .75),
      text = ~paste0(Feature, "<br>", logFC_1, " / ", logFC_2), hoverinfo = "text") %>%
    ply(xaxis = list(title = o$xlab, zeroline = TRUE), yaxis = list(title = o$ylab, zeroline = TRUE),
        legend = list(orientation = "h", y = 1.1), title = list(text = o$ttl, x = .02, font = list(size = 13)))
}

byod_assoc <- function(d, feature, var, at_tp = "All", cols = NULL) {
  ylab <- paste0(feature, " (", d$unit, ")")
  fail <- function(msg) list(plot = emptyplot(msg), summary = data.frame(Group = factor()), splits = character(0),
    ylab = ylab, ttl = msg, cols = setNames(character(0), character(0)), gg = gg_clinbox, csv = data.frame(), kind = "box")
  v <- expr_vec(d, feature); if (is.null(v)) return(fail(paste(feature, "not found")))
  mt <- d$meta; sel <- if (at_tp == "All") rep(TRUE, nrow(mt)) else as.character(mt$Timepoint) == at_tp
  x <- mt[[var]][sel]; y <- v[sel]; tp <- as.character(mt$Timepoint[sel])
  subj <- d$subject; sj <- if (!is.null(subj) && subj %in% names(mt)) as.character(mt[[subj]][sel]) else NULL
  ok <- !is.na(x) & !is.na(y); x <- x[ok]; y <- y[ok]; tp <- tp[ok]; if (!is.null(sj)) sj <- sj[ok]

  donor_lvl <- FALSE
  if (!is.null(sj) && length(sj) && all(tapply(as.character(x), sj, function(z) length(unique(z[!is.na(z)])) <= 1))) {
    du <- unique(sj); donor_lvl <- TRUE
    yn <- vapply(du, function(u) mean(y[sj == u], na.rm = TRUE), numeric(1))
    x  <- if (is.numeric(x)) vapply(du, function(u) { xx <- x[sj == u]; xx[!is.na(xx)][1] }, numeric(1))
          else vapply(du, function(u) { xx <- as.character(x[sj == u]); xx[!is.na(xx)][1] }, character(1))
    y <- yn; names(y) <- du; tp <- rep("subject", length(du)); sj <- du
  }
  unit_n <- if (donor_lvl) "subjects" else "samples"
  if (length(y) < 3) return(fail(paste("Too few non-missing", unit_n)))

  n_subj <- if (is.null(sj)) 0L else length(unique(sj))
  if (is.numeric(x)) {
    grp <- if (donor_lvl) factor(rep("subject-level", length(y))) else factor(tp, levels = intersect(d$tp_levels, unique(tp)))
    lv <- levels(grp)
    cc <- if (is.null(cols)) setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = if (donor_lvl) "cluster" else "Timepoint")), lv) else cols[lv]
    ttl <- if (donor_lvl) {
      ct <- suppressWarnings(cor.test(x, y))
      sprintf("r = %.2f, p = %s, n = %d subjects", unname(ct$estimate), signif(ct$p.value, 2), length(y))
    } else if (n_subj >= 3L) {

      sf  <- factor(sj)
      tf2 <- droplevels(factor(tp))
      adj <- if (nlevels(tf2) > 1L) ~ sf + tf2 else ~ sf
      rw <- suppressWarnings(stats::cor(
        residuals(stats::lm(update(adj, x ~ .))), residuals(stats::lm(update(adj, y ~ .)))))
      pv <- tryCatch(summary(stats::lm(update(adj, y ~ x + .)))$coefficients["x", "Pr(>|t|)"],
                     error = function(e) NA_real_)
      sprintf("within-subject%s r = %.2f, p = %s, %d profiles from %d subjects",
              if (nlevels(tf2) > 1L) ", time-adjusted" else "",
              rw, signif(pv, 2), length(y), n_subj)
    } else {
      sprintf("r = %.2f, descriptive only: no subject declared to adjust for (%d profiles)",
              suppressWarnings(stats::cor(x, y)), length(y))
    }
    p <- plot_ly(x = x, y = y, color = grp, colors = unname(cc[lv]), type = "scatter", mode = "markers",
        marker = list(size = 8), text = paste0(var, "=", signif(x, 3), "<br>", tp), hoverinfo = "text") %>%
      ply(xaxis = list(title = var), yaxis = list(title = ylab), legend = list(orientation = "h", y = 1.12),
          title = list(text = ttl, x = .02, font = list(size = 13)))
    return(list(plot = p, summary = data.frame(x = x, y = y, group = as.character(grp), stringsAsFactors = FALSE),
      splits = lv, ylab = ylab, ttl = ttl, cols = cc, kind = "scatter",
      gg = function(s, yl, tt, co) gg_clin_scatter(s, yl, tt, co, var),
      csv = data.frame(sample = names(y), value = unname(x), expression = unname(y), timepoint = tp, stringsAsFactors = FALSE)))
  }
  xf <- factor(as.character(x)); lv <- levels(xf)
  cc <- if (is.null(cols)) setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = var)), lv) else cols[lv]
  ttl <- if (donor_lvl) {
    pval <- tryCatch(summary(stats::aov(y ~ xf))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_)
    sprintf("ANOVA p = %s (n = %d subjects)", signif(pval, 2), length(y))
  } else if (n_subj >= 3L) {

    sf  <- factor(sj); pval <- NA_real_
    tf2 <- droplevels(factor(tp))
    aliased <- nlevels(tf2) < 2L ||
      all(tapply(as.character(xf), tf2, function(z) length(unique(z))) == 1L)
    dfa <- data.frame(y = y, xf = xf, sf = sf, tf2 = tf2)
    fit <- tryCatch(if (aliased) stats::aov(y ~ xf + Error(sf), data = dfa)
                    else        stats::aov(y ~ xf + tf2 + Error(sf), data = dfa),
                    error = function(e) NULL)
    if (!is.null(fit)) for (st in summary(fit)) {
      tb <- tryCatch(st[[1]], error = function(e) NULL)
      if (is.null(tb) || is.null(rownames(tb))) next
      hit <- trimws(rownames(tb)) == "xf"
      if (any(hit) && "Pr(>F)" %in% colnames(tb)) pval <- tb[hit, "Pr(>F)"][1]
    }
    if (is.finite(pval))
      sprintf("within-subject%s ANOVA p = %s (%d profiles from %d subjects)",
              if (aliased) "" else ", time-adjusted", signif(pval, 2), length(y), n_subj)
    else
      sprintf("descriptive only: the within-subject model is not estimable (%d profiles from %d subjects)",
              length(y), n_subj)
  } else {
    sprintf("descriptive only: no subject declared to adjust for (%d profiles)", length(y))
  }
  p <- plot_ly(x = xf, y = y, color = xf, colors = unname(cc[lv]), type = "box", boxpoints = "all", jitter = .4, pointpos = 0) %>%
    ply(xaxis = list(title = var), yaxis = list(title = ylab), showlegend = FALSE,
        title = list(text = ttl, x = .02, font = list(size = 13)))
  summ <- data.frame(Group = factor(lv, levels = lv)); attr(summ, "raw") <- data.frame(Group = as.character(xf), value = y, stringsAsFactors = FALSE)
  list(plot = p, summary = summ, splits = lv, ylab = ylab, ttl = ttl, cols = cc, kind = "box", gg = gg_clinbox,
       csv = data.frame(group = as.character(xf), expression = unname(y), stringsAsFactors = FALSE))
}

# Between-donor traits use one value per donor, not one value per biopsy.
byod_donor_collapse <- function(d, var, at_tp = "All") {
  subj_col <- d$subject %||% "Kidney_ID"
  subj <- d$meta[[subj_col]]
  if (is.null(subj) || !(var %in% names(d$meta)))
    return(list(ok = FALSE, msg = "Needs a subject/pairing column and the selected variable."))
  tp <- as.character(d$meta$Timepoint)
  keep <- if (identical(at_tp, "All")) seq_len(ncol(d$expr)) else which(tp == at_tp)
  if (!length(keep)) return(list(ok = FALSE, msg = paste0("No samples at timepoint ", at_tp, ".")))
  sj <- as.character(subj)[keep]; E <- d$expr[, keep, drop = FALSE]
  cvr <- d$meta[[var]][keep]; isnum <- is.numeric(cvr)
  donors <- unique(sj)
  Ed <- vapply(donors, function(dn) {
    cc <- which(sj == dn)
    if (length(cc) == 1L) E[, cc] else rowMeans(E[, cc, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(E)))
  dimnames(Ed) <- list(rownames(E), donors)
  clin <- vapply(donors, function(dn) {
    z <- cvr[sj == dn]; z <- z[!is.na(z)]
    if (!length(z)) return(if (isnum) NA_real_ else NA_character_)
    if (isnum) as.numeric(z[1]) else as.character(z[1])
  }, if (isnum) numeric(1) else character(1))
  list(ok = TRUE, Ed = Ed, clin = clin, isnum = isnum, donors = donors, subj_col = subj_col)
}

.cor_rp <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y); n <- sum(ok)
  if (n < 3L) return(c(rho = NA_real_, p = NA_real_, n = n))
  r <- suppressWarnings(stats::cor(x[ok], y[ok], method = method))
  if (!is.finite(r)) return(c(rho = NA_real_, p = NA_real_, n = n))
  if (abs(r) >= 1) return(c(rho = r, p = 0, n = n))
  tv <- r * sqrt((n - 2) / (1 - r^2)); c(rho = r, p = 2 * stats::pt(-abs(tv), n - 2), n = n)
}

byod_clin_corr <- function(d, var, at_tp = "All", method = "spearman", min_n = 8L, adjust = "BH") {
  cc <- byod_donor_collapse(d, var, at_tp)
  if (!isTRUE(cc$ok)) return(c(cc, list(isnum = NA)))
  Ed <- cc$Ed; clin <- cc$clin
  fail <- function(m) list(ok = FALSE, msg = m, isnum = cc$isnum)
  if (cc$isnum) {
    okc <- is.finite(clin)
    if (sum(okc) < min_n) return(fail(sprintf("Only %d donors have a non-missing %s (need >= %d).", sum(okc), var, min_n)))
    Ed <- Ed[, okc, drop = FALSE]; cl <- clin[okc]; n <- length(cl)
    r <- suppressWarnings(as.numeric(stats::cor(t(Ed), cl, method = method)))
    r[!is.finite(r)] <- NA
    tv <- r * sqrt((n - 2) / (1 - r^2))
    p <- 2 * stats::pt(-abs(tv), n - 2); p[is.finite(r) & abs(r) >= 1] <- 0
    df <- data.frame(Feature = rownames(Ed), logFC = r, AveExpr = rowMeans(Ed, na.rm = TRUE),
                     t = NA_real_, P.Value = p, adj.P.Val = p.adjust(p, adjust), stringsAsFactors = FALSE)
    df <- df[order(df$adj.P.Val, -abs(df$logFC)), ]
    ylab <- if (method == "spearman") "Spearman rho" else "Pearson r"
    return(list(ok = TRUE, isnum = TRUE, tab = df, summary = df, method = method, var = var, at_tp = at_tp,
                n_donors = n, ylab = ylab,
                ttl = sprintf("%s correlation with %s (n=%d donors, %s)", ylab, var, n,
                              if (identical(at_tp, "All")) "donor-averaged" else at_tp)))
  }

  okc <- !is.na(clin); lv <- levels(droplevels(factor(as.character(clin)[okc])))
  if (length(lv) != 2) return(fail(sprintf("'%s' must have exactly 2 non-missing groups (has %d: %s).",
                                           var, length(lv), paste(lv, collapse = "/"))))
  Ed <- Ed[, okc, drop = FALSE]; grp <- factor(as.character(clin)[okc], levels = lv)
  if (min(table(grp)) < 3) return(fail(sprintf("Too few donors per group (%s).",
                                    paste(names(table(grp)), as.integer(table(grp)), sep = "=", collapse = ", "))))
  design <- stats::model.matrix(~ grp)
  fit <- limma::eBayes(limma::lmFit(Ed, design), trend = TRUE)
  tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
  df <- data.frame(Feature = rownames(tt), logFC = tt$logFC, AveExpr = tt$AveExpr, t = tt$t,
                   P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  df <- df[order(df$adj.P.Val, -abs(df$logFC)), ]
  list(ok = TRUE, isnum = FALSE, tab = df, summary = df, method = "limma", var = var, at_tp = at_tp,
       n_donors = sum(okc), levels = lv, ylab = "Mean difference (analysis scale)",
       ttl = sprintf("%s: %s vs %s (n=%d donors, %s)", var, lv[2], lv[1], sum(okc),
                     if (identical(at_tp, "All")) "donor-averaged" else at_tp))
}

byod_clin_volcano <- function(o, fdr = 0.05, eff = 0, cols = NULL, n_label = 8) {
  if (!isTRUE(o$ok) || !nrow(o$summary)) return(list(plot = emptyplot(o$msg %||% "no result"),
    summary = data.frame(), splits = VOLCANO_LEVELS, ylab = "-log10 P", ttl = o$msg %||% "", cols = cols))
  shown <- o$summary[is.finite(o$summary$P.Value) & is.finite(o$summary$logFC), , drop = FALSE]
  v <- plot_volcano(shown, o$ttl, fdr = fdr, lfc = eff, cols = cols, n_label = n_label)
  if (isTRUE(o$isnum)) v$plot <- v$plot %>% plotly::layout(xaxis = list(title = o$ylab))
  v
}

byod_clin_scatter <- function(d, feature, var, at_tp = "All", method = "spearman", cols = NULL) {
  ylab <- paste0(feature, " (", d$unit, ")")
  fail <- function(msg) list(plot = emptyplot(msg), summary = data.frame(Group = factor()), splits = character(0),
    ylab = ylab, ttl = msg, cols = setNames(character(0), character(0)), gg = gg_clinbox, csv = data.frame(), kind = "box")
  cc <- byod_donor_collapse(d, var, at_tp); if (!isTRUE(cc$ok)) return(fail(cc$msg))
  if (!(feature %in% rownames(cc$Ed))) return(fail(paste(feature, "not found")))
  ex <- cc$Ed[feature, ]; cl <- cc$clin
  if (cc$isnum) {
    ok <- is.finite(ex) & is.finite(cl); if (sum(ok) < 3) return(fail("Too few donors with non-missing values"))
    rp <- .cor_rp(ex, cl, method); lab <- if (method == "spearman") "Spearman rho" else "Pearson r"
    ttl <- sprintf("%s vs %s: %s = %.2f, p = %s (n=%d donors)", feature, var, lab, rp["rho"], signif(rp["p"], 2), rp["n"])
    p <- plot_ly(x = cl[ok], y = ex[ok], type = "scatter", mode = "markers",
        marker = list(size = 9, color = "#2C8C99", opacity = .8),
        text = paste0(cc$donors[ok], "<br>", var, "=", signif(cl[ok], 3)), hoverinfo = "text") %>%
      ply(xaxis = list(title = var), yaxis = list(title = ylab), showlegend = FALSE,
          title = list(text = ttl, x = .02, font = list(size = 13)))
    return(list(plot = p, summary = data.frame(x = cl[ok], y = ex[ok], group = "donor", stringsAsFactors = FALSE),
      splits = "donor", ylab = ylab, ttl = ttl, cols = setNames("#2C8C99", "donor"), kind = "scatter",
      gg = function(s, yl, tt, co) gg_clin_scatter(s, yl, tt, co, var),
      csv = data.frame(donor = cc$donors[ok], clinical = unname(cl[ok]), expression = unname(ex[ok]), stringsAsFactors = FALSE)))
  }
  ok <- !is.na(cl); xf <- factor(as.character(cl)[ok]); lv <- levels(xf); yv <- ex[ok]
  if (length(lv) < 2) return(fail("Need >= 2 groups"))
  cols2 <- if (is.null(cols)) setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = var)), lv) else cols[lv]
  pval <- tryCatch(t.test(yv ~ xf)$p.value, error = function(e) NA)
  ttl <- sprintf("%s by %s: t-test p = %s (n=%d donors)", feature, var, signif(pval, 2), sum(ok))
  p <- plot_ly(x = xf, y = yv, color = xf, colors = unname(cols2[lv]), type = "box", boxpoints = "all", jitter = .5, pointpos = 0) %>%
    ply(xaxis = list(title = var), yaxis = list(title = ylab), showlegend = FALSE,
        title = list(text = ttl, x = .02, font = list(size = 13)))
  summ <- data.frame(Group = factor(lv, levels = lv)); attr(summ, "raw") <- data.frame(Group = as.character(xf), value = yv, stringsAsFactors = FALSE)
  list(plot = p, summary = summ, splits = lv, ylab = ylab, ttl = ttl, cols = cols2, kind = "box", gg = gg_clinbox,
       csv = data.frame(donor = cc$donors[ok], group = as.character(xf), expression = unname(yv), stringsAsFactors = FALSE))
}

byod_clin_enrich <- function(o, d, collection, mode = "gsea", fdr = 0.05, top_n = 200L, cols = NULL) {
  if (!isTRUE(o$ok) || !nrow(o$tab)) return(list(plot = emptyplot("run the correlation first"),
    summary = data.frame(), splits = character(0), ylab = "", ttl = "", cols = cols))
  tabb <- o$tab
  if (mode == "gsea") {
    rk <- if (isTRUE(o$isnum)) setNames(tabb$logFC, tabb$Feature) else setNames(tabb$t, tabb$Feature)
    rk <- rk[is.finite(rk)]; rk <- rk[!duplicated(names(rk))]
    return(plot_gsea(run_gsea(rk, collection, d = d), 20, paste0(collection, " - GSEA vs ", o$var), cols = cols))
  }
  genes <- head(tabb$Feature[is.finite(tabb$adj.P.Val) & tabb$adj.P.Val < fdr], top_n)
  plot_ora(run_ora(genes, collection, rownames(d$expr), d = d), 20,
           paste0(collection, " ORA (", length(genes), " genes, FDR<", fdr, ")"), cols = cols)
}

byod_fit_error <- function(msg) structure(list(msg = msg), class = "byod_unfittable")

byod_fixed_fit <- function(E, design, testcoef, trend = TRUE) {
  resid_df <- ncol(E) - qr(design)$rank
  if (resid_df < 1)
    stop(sprintf("This design is saturated: %d profiles against %d estimable design columns leaves %d residual degrees of freedom, so no variance can be estimated. Every combination of the declared factors has a single profile; the model needs replication within at least one of them.",
                 ncol(E), qr(design)$rank, resid_df), call. = FALSE)
  fit <- limma::eBayes(limma::lmFit(E, design), trend = trend)
  sb <- if (length(testcoef) >= 2) "F" else "p"
  limma::topTable(fit, coef = testcoef, number = Inf, sort.by = sb)
}

MIN_UNITS <- list(
  group_contrast   = 3L,   # independent units per group for any group comparison
  interaction      = 3L,   # units per group for a group x time interaction
  temporal         = 3L,   # units contributing to a temporal trajectory
  coexpression     = 8L    # minimum independent units permitted for network analysis
)

byod_count_units <- function(meta, fac, subject_col = NULL) {
  fac <- as.character(fac)
  if (!is.null(subject_col) && subject_col %in% names(meta)) {
    su <- as.character(meta[[subject_col]])
    ok <- !is.na(su) & !is.na(fac)
    if (any(ok)) {
      tab <- tapply(fac[ok], su[ok], function(z) unique(z)[1])
      return(table(factor(tab, levels = sort(unique(fac[ok])))))
    }
  }
  table(factor(fac[!is.na(fac)]))
}

byod_unit_floor_msg <- function(counts, floor_n, what, n_samples = NULL) {
  sprintf("%s needs at least %d independent units per group; this design has %s%s. Add units, or analyse at a coarser resolution.",
          what, floor_n,
          paste(sprintf("%s=%d", names(counts), as.integer(counts)), collapse = ", "),
          if (!is.null(n_samples)) sprintf(" (from %d profiles)", n_samples) else "")
}

byod_gti_cands <- function(d) {
  v <- setdiff(c(if (!is.null(d$meta$WI_status)) "WI_status", d$clinical_cat), c("Timepoint", "Timepoint_WI"))
  v[vapply(v, function(x) length(group_levels(d, x)) == 2, logical(1))]
}

byod_find_pair_block <- function(meta_kept, s, g, subj_col, group_var) {
  if (is.null(s)) return(NULL)
  ns <- nlevels(s); glv <- levels(droplevels(g))

  if (any(tapply(as.character(g), s, function(z) length(unique(z))) != 1)) return(NULL)
  subj_grp <- tapply(as.character(g), s, function(z) z[1])          # one group per subject
  cand <- setdiff(names(meta_kept), c("Timepoint", "Timepoint_WI", subj_col, group_var))
  best <- NULL; best_lv <- Inf
  for (bcol in cand) {
    x <- meta_kept[[bcol]]; if (is.numeric(x)) next
    blk <- droplevels(factor(as.character(x))); if (anyNA(blk)) next
    nb <- nlevels(blk); if (nb < 2 || nb >= ns) next                                  # coarser than subject
    if (any(tapply(as.character(blk), s, function(z) length(unique(z))) != 1)) next   # subject nests in blk

    subj_blk <- tapply(as.character(blk), s, function(z) z[1])       # one block per subject
    ok <- 0L; tot <- 0L
    for (bl in levels(blk)) {
      subj_in <- names(subj_blk)[subj_blk == bl]; if (!length(subj_in)) next
      tot <- tot + 1L
      cnt <- table(factor(subj_grp[subj_in], levels = glv))
      if (all(cnt >= 1L) && length(unique(as.integer(cnt))) == 1L && max(cnt) <= 2L) ok <- ok + 1L
    }
    if (tot < 4L || ok / tot < 0.8) next                            # enough, mostly-balanced matched sets
    if (nb < best_lv) { best <- bcol; best_lv <- nb }
  }
  best
}

byod_gti <- function(d, group_var, min_per = 3L) {
  scol <- d$subject; subj <- if (!is.null(scol) && scol %in% names(d$meta)) factor(as.character(d$meta[[scol]])) else NULL
  has_subj <- !is.null(subj) && !all(is.na(subj))                 # uploads may have no subject column
  gf <- droplevels(grp_factor(d, group_var)); lv <- levels(gf)
  if (length(lv) != 2) return(list(ok = FALSE, msg = sprintf("'%s' needs exactly 2 groups (has %d).", group_var, length(lv))))
  timef <- droplevels(factor(d$meta$Timepoint, levels = d$tp_levels))
  keep <- !is.na(gf) & !is.na(timef) & (if (has_subj) !is.na(subj) else TRUE)
  E <- d$expr[, keep, drop = FALSE]; g <- droplevels(gf[keep]); tf <- droplevels(timef[keep])
  sj <- if (has_subj) droplevels(subj[keep]) else NULL
  if (nlevels(tf) < 2) return(list(ok = FALSE, msg = "Need >= 2 timepoints for an interaction."))

  if (nlevels(g) != 2)
    return(list(ok = FALSE, msg = sprintf("'%s' has %d group(s) among the %d profiles with both a timepoint and a subject, so there is nothing to contrast over time.",
                                          group_var, nlevels(g), sum(keep))))

  gt <- table(g, tf) > 0
  if (any(rowSums(gt) < 2))
    return(list(ok = FALSE, msg = sprintf("Group '%s' is present at only one timepoint, so '%s' and time are confounded and no interaction is estimable.",
                                          rownames(gt)[which(rowSums(gt) < 2)[1]], group_var)))

  u_counts <- byod_count_units(d$meta[keep, , drop = FALSE], g, if (has_subj) scol else NULL)
  if (min(u_counts) < MIN_UNITS$interaction)
    return(list(ok = FALSE, msg = byod_unit_floor_msg(u_counts, MIN_UNITS$interaction,
                                                      "A group x time interaction", sum(keep))))
  if (min(table(g)) < min_per) return(list(ok = FALSE, msg = sprintf("Too few samples in a group (%s).",
                                     paste(names(table(g)), as.integer(table(g)), sep = "=", collapse = ", "))))

  blocked <- has_subj && nlevels(sj) >= 2 && nlevels(sj) < ncol(E) &&
             sum(table(sj) >= 2) >= 2
  blk <- if (blocked) byod_find_pair_block(d$meta[keep, , drop = FALSE], sj, g, scol, group_var) else NULL
  paired <- FALSE
  if (!is.null(blk)) {                                            # within-pair design D
    ptf <- droplevels(factor(paste(as.character(d$meta[[blk]])[keep], as.character(tf), sep = "_")))
    kid <- sj
    design <- stats::model.matrix(~ ptf + kid + g:tf); paired <- TRUE
    n_units <- nlevels(droplevels(factor(as.character(d$meta[[blk]])[keep])))
  } else {                                                        # subject-blocked (A) or unblocked
    design <- if (blocked) stats::model.matrix(~ sj + tf + tf:g) else stats::model.matrix(~ g * tf)
    n_units <- if (has_subj) nlevels(sj) else sum(keep)
  }
  design <- design[, colSums(design != 0) > 0, drop = FALSE]
  qrd <- qr(design)   # drop aliased columns (keep the nuisance blocks) rather than pseudoreplicate
  if (qrd$rank < ncol(design)) design <- design[, sort(qrd$pivot[seq_len(qrd$rank)]), drop = FALSE]
  int <- grep(":", colnames(design), fixed = TRUE)

  if (paired && (length(int) < nlevels(tf) - 1L || (ncol(E) - ncol(design)) < 30L)) {
    paired <- FALSE; n_units <- nlevels(sj)
    design <- stats::model.matrix(~ sj + tf + tf:g)
    design <- design[, colSums(design != 0) > 0, drop = FALSE]
    qrd <- qr(design); if (qrd$rank < ncol(design)) design <- design[, sort(qrd$pivot[seq_len(qrd$rank)]), drop = FALSE]
    int <- grep(":", colnames(design), fixed = TRUE)
  }
  if (!length(int)) return(list(ok = FALSE, msg = "No interaction terms estimable (unbalanced group x timepoint)."))

  tt <- tryCatch(byod_fixed_fit(E, design, int), error = function(e) e)
  if (inherits(tt, "error")) return(list(ok = FALSE, msg = conditionMessage(tt)))
  sc <- if ("F" %in% names(tt)) "F" else "t"
  df <- data.frame(Feature = rownames(tt), stat = tt[[sc]], P.Value = tt$P.Value, adj.P.Val = tt$adj.P.Val, stringsAsFactors = FALSE)
  names(df)[2] <- sc; df <- df[order(df$adj.P.Val, -abs(df[[2]])), ]

  lbl <- if (paired) sprintf("n=%d %s%s, within-%s paired", n_units, blk, if (n_units == 1L) "" else "s", blk)
         else if (blocked) sprintf("n=%d subjects, subject-adjusted", n_units) else sprintf("n=%d profiles", sum(keep))
  list(ok = TRUE, tab = df, group_var = group_var, levels = lv, stat = sc, blocked = blocked,
       paired = paired, pair_block = if (paired) blk else NA_character_, n_units = n_units, n = sum(keep),

       unit = if (paired) blk else (scol %||% "profile"),
       ttl = sprintf("%s x time interaction (%s vs %s, %s; %d profiles)", group_var, lv[2], lv[1], lbl, sum(keep)))
}
byod_gti_enrich <- function(o, d, collection, fdr = 0.05, top_n = 200L, cols = NULL, display_n = 10L) {
  if (!isTRUE(o$ok)) return(list(plot = emptyplot("run the interaction first"), summary = data.frame(), splits = character(0), ylab = "", ttl = "", cols = cols))
  genes <- head(o$tab$Feature[is.finite(o$tab$adj.P.Val) & o$tab$adj.P.Val < fdr], top_n)
  plot_ora(run_ora(genes, collection, rownames(d$expr), d = d), display_n,
           sprintf("%s ORA (%d interaction genes, FDR<%s)", collection, length(genes), fdr), cols = cols)
}

byod_outcome_plan <- function(module, o, pr = list()) {
  st  <- byod_outcome_of(o)
  eng <- pr$.engine; pr$.engine <- NULL
  if (identical(st$status, "declined"))
    return(list(list(module = module, params = c(pr, list(reason = st$reason)),
                     outcome = "declined", engine = eng)))
  if (identical(st$status, "stepped down"))
    pr <- c(pr, list(step_down = st$reason,
                     requested_method = o$requested_method %||% NA,
                     executed_method  = o$executed_method  %||% NA))
  out <- list(list(module = module, params = pr, outcome = st$status, engine = eng))

  if (is.list(o) && length(o$declined))
    for (nm in names(o$declined))
      out[[length(out) + 1L]] <- list(
        module = module,
        params = c(pr, list(sub = nm, reason = as.character(o$declined[[nm]])[1])),
        outcome = "declined", engine = eng)
  out
}

byod_outcome_of <- function(o) {
  no <- function(why) list(status = "declined", reason = why)
  if (is.null(o))                       return(no("the analysis returned no result"))
  if (inherits(o, c("error", "try-error")))
    return(no(conditionMessage(attr(o, "condition") %||% simpleError(as.character(o)))))
  dec <- attr(o, "declined")
  if (!is.null(dec) && length(dec) && nzchar(as.character(dec)[1]))
    return(no(as.character(dec)[1]))
  if (is.list(o) && !is.null(o$ok) && !isTRUE(o$ok))
    return(no(o$msg %||% "the resolved design does not support it"))
  if (is.data.frame(o) && !nrow(o))     return(no("the analysis produced no rows"))
  if (is.matrix(o) && !length(o))       return(no("the analysis produced an empty matrix"))
  if (is.list(o) && !is.data.frame(o)) {

    pay <- c("mat", "tab", "res", "z", "ME", "colors", "summary")
    have <- intersect(pay, names(o))
    if (length(have) && all(vapply(have, function(k) {
          v <- o[[k]]
          is.null(v) || (is.data.frame(v) && !nrow(v)) ||
            ((is.matrix(v) || is.array(v)) && !length(v)) ||
            (is.atomic(v) && !length(v))
        }, logical(1))))
      return(no("the analysis produced an empty result"))
  }
  if (isTRUE(o$stepped_down))
    return(list(status = "stepped down",
                reason = o$step_down_reason %||% "a weaker but estimable method was used"))
  list(status = "completed", reason = NULL)
}

collapse_layer_to_unit <- function(d, keep = NULL) {
  ucol <- resolved_unit_col(d, keep = keep)
  if (is.null(ucol) || identical(ucol, d$subject %||% "")) {
    n <- if (!is.null(ucol)) length(unique(as.character(d$meta[[ucol]]))) else ncol(d$expr)
    return(list(d = d, unit = ucol %||% "profile", n_units = n, collapsed = FALSE))
  }
  mt <- d$meta
  keep <- intersect(unique(c("Timepoint", keep)), names(mt))
  key  <- interaction(lapply(c(ucol, keep), function(k) as.character(mt[[k]])),
                      drop = TRUE, sep = "\r")
  idx  <- split(seq_len(ncol(d$expr)), key)
  E2   <- vapply(idx, function(i) rowMeans(d$expr[, i, drop = FALSE], na.rm = TRUE),
                 numeric(nrow(d$expr)))
  rownames(E2) <- rownames(d$expr)
  parts <- do.call(rbind, strsplit(names(idx), "\r", fixed = TRUE))
  m2 <- as.data.frame(parts, stringsAsFactors = FALSE); names(m2) <- c(ucol, keep)
  rownames(m2) <- colnames(E2) <- names(idx)
  m2$Timepoint <- factor(m2$Timepoint, levels = d$tp_levels)
  d2 <- d; d2$expr <- E2; d2$meta <- m2; d2$subject <- ucol; d2$nest <- NULL
  list(d = d2, unit = ucol, n_units = length(unique(m2[[ucol]])), collapsed = TRUE)
}

byod_temporal <- function(d, mode = "factor", fdr = 0.05, k = 6L, min_tp = 3L, max_heat = 300L, core_cut = 0,
                          seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  timef <- droplevels(factor(d$meta$Timepoint, levels = d$tp_levels))
  if (nlevels(timef) < min_tp) return(list(ok = FALSE, msg = sprintf("Need >= %d timepoints (has %d).", min_tp, nlevels(timef))))

  n_profiles_in <- ncol(d$expr)
  cl <- collapse_layer_to_unit(d)
  d  <- cl$d
  timef <- droplevels(factor(d$meta$Timepoint, levels = d$tp_levels))
  subj <- as.character(d$meta[[d$subject %||% "Kidney_ID"]])

  n_units_t <- cl$n_units
  if (n_units_t < MIN_UNITS$temporal)
    return(list(ok = FALSE, msg = sprintf("Temporal clustering needs at least %d independent units; this design has %d (from %d profiles).",
                                          MIN_UNITS$temporal, n_units_t, ncol(d$expr))))
  note <- NULL; use_spline <- identical(mode, "spline")

  subjf <- factor(subj)
  blocked <- nlevels(subjf) >= 2 && nlevels(subjf) < ncol(d$expr) &&
             sum(table(subjf) >= 2) >= 2
  if (use_spline) {
    tnum <- as.numeric(timef)   # ordinal (evenly-spaced) time, like the GAM trajectory

    dfree <- max(2, min(3, nlevels(timef) - 1))
    B <- NULL
    for (dtry in seq(dfree, 1)) {
      B <- tryCatch(splines::ns(tnum, df = dtry), error = function(e) NULL)
      if (!is.null(B)) { dfree <- dtry; break }
    }
    if (is.null(B))
      return(list(ok = FALSE, msg = sprintf(
        "A smooth trend cannot be fitted: the %d timepoints are too unevenly sampled for a spline basis (%s). Use the factor mode, which makes no assumption about spacing.",
        nlevels(timef),
        paste(sprintf("%s n=%d", levels(timef), as.integer(table(timef))), collapse = ", "))))
    colnames(B) <- paste0("ns", seq_len(ncol(B)))
    design <- if (blocked) stats::model.matrix(~ subjf + B) else stats::model.matrix(~ B)
    note <- "Smooth natural-spline trend over the ordered timepoints (like the GAM trajectory)."
  } else design <- if (blocked) stats::model.matrix(~ subjf + timef) else stats::model.matrix(~ timef)
  if (nrow(design) != ncol(d$expr))
    return(list(ok = FALSE, msg = sprintf("%d of %d profiles have no usable timepoint%s, so the model cannot be fitted. Fix those rows in the metadata and rebuild the study.",
                                          ncol(d$expr) - nrow(design), ncol(d$expr),
                                          if (blocked) " or subject label" else "")))
  design <- design[, colSums(design != 0) > 0, drop = FALSE]
  qrd <- qr(design)   # drop any aliased columns, keep the subject adjustment
  if (qrd$rank < ncol(design)) design <- design[, sort(qrd$pivot[seq_len(qrd$rank)]), drop = FALSE]
  testcoef <- which(!grepl("^\\(Intercept\\)$|^subjf", colnames(design)))
  if (!length(testcoef)) return(list(ok = FALSE, msg = "No temporal terms estimable."))
  tt <- tryCatch(byod_fixed_fit(d$expr, design, testcoef), error = function(e) e)
  if (inherits(tt, "error")) return(list(ok = FALSE, msg = conditionMessage(tt)))
  tt$Feature <- rownames(tt)
  sig <- tt$Feature[is.finite(tt$adj.P.Val) & tt$adj.P.Val < fdr]
  if (length(sig) < max(k + 1L, 10L)) return(list(ok = FALSE, msg = sprintf("Only %d features change over time at FDR<%s (too few for %d clusters).", length(sig), fdr, k)))
  lvl <- levels(timef)
  gm <- t(vapply(sig, function(f) tapply(d$expr[f, ], timef, mean, na.rm = TRUE)[lvl], numeric(length(lvl))))
  colnames(gm) <- lvl; z <- t(scale(t(gm))); z <- z[stats::complete.cases(z), , drop = FALSE]
  kk <- min(as.integer(k), nrow(z) - 1L)

  fz <- if (requireNamespace("e1071", quietly = TRUE))
          tryCatch(e1071::cmeans(z, centers = kk, m = 1.5, iter.max = 200), error = function(e) NULL) else NULL
  if (is.null(fz)) { km <- stats::kmeans(z, centers = kk, nstart = 10, iter.max = 50); hard <- km$cluster; member <- rep(1, nrow(z)); cmethod <- "k-means" }
  else { hard <- fz$cluster; member <- fz$membership[cbind(seq_len(nrow(z)), hard)]; cmethod <- "fuzzy c-means" }
  peak <- tapply(seq_len(nrow(z)), hard, function(i) mean(max.col(z[i, , drop = FALSE], ties.method = "first")))
  relab <- setNames(order(order(peak)), names(peak)); clu <- paste0("C", relab[as.character(hard)])
  clusters <- data.frame(Feature = rownames(z), cluster = clu, membership = round(member, 3), stringsAsFactors = FALSE)
  fcol <- if ("F" %in% names(tt)) "F" else "t"
  tab <- tt[tt$Feature %in% rownames(z), c("Feature", fcol, "P.Value", "adj.P.Val")]
  tab$cluster <- clusters$cluster[match(tab$Feature, clusters$Feature)]
  tab$membership <- clusters$membership[match(tab$Feature, clusters$Feature)]; tab <- tab[order(tab$adj.P.Val), ]
  n_core <- sum(member >= core_cut)
  list(ok = TRUE, tab = tab, z = z, clusters = clusters, timepoints = lvl, mode = if (use_spline) "spline" else "factor",
       note = note, n_sig = length(sig), n_clustered = nrow(z), k = length(unique(clu)), max_heat = max_heat, fdr = fdr,
       unit = cl$unit, n_units = cl$n_units, n_profiles = n_profiles_in, collapsed = cl$collapsed,
       blocked = blocked, core_cut = core_cut, n_core = n_core, cmethod = cmethod,
       ttl = sprintf("%d temporal programs from %d features (FDR<%s, %s%s%s)", length(unique(clu)), nrow(z), fdr, cmethod,
                     if (blocked) ", subject-adjusted" else "", if (core_cut > 0) sprintf(", core>=%.2f n=%d", core_cut, n_core) else ""))
}
.tp_members <- function(o) { m <- o$clusters; if ((o$core_cut %||% 0) > 0) m <- m[m$membership >= o$core_cut, , drop = FALSE]; m }
byod_temporal_profiles <- function(o, cols = NULL) {
  cm <- .tp_members(o); z <- o$z[cm$Feature, , drop = FALSE]; cl <- cm$cluster; tps <- o$timepoints; rows <- list()
  for (cc in sort(unique(cl))) { zz <- z[cl == cc, , drop = FALSE]
    rows[[cc]] <- data.frame(Timepoint = factor(tps, levels = tps), Group = paste0(cc, " (n=", nrow(zz), ")"),
      mean = colMeans(zz, na.rm = TRUE), se = apply(zz, 2, function(x) stats::sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))),
      n = nrow(zz), stringsAsFactors = FALSE) }
  st <- do.call(rbind, rows); rownames(st) <- NULL; splits <- unique(as.character(st$Group))
  cc2 <- if (is.null(cols)) resolve_cols(splits, PAL_NAMES[1], "cluster") else cols
  plt <- tc_plotly(list(summary = st, splits = splits), "mean z-score", o$ttl, cc2) %>%
    plotly::layout(legend = list(orientation = "h", y = -0.18, x = 0), margin = list(b = 70, t = 8))
  list(plot = plt, summary = st, splits = splits, ylab = "mean z-score", ttl = o$ttl, cols = cc2)
}
byod_temporal_heat <- function(o, colorscale = NULL) {
  cm <- .tp_members(o); z <- o$z[cm$Feature, , drop = FALSE]; clv <- cm$cluster; memb <- cm$membership
  if (nrow(z) > o$max_heat) { ord0 <- order(-memb)[seq_len(o$max_heat)]; z <- z[ord0, , drop = FALSE]; clv <- clv[ord0] }
  ord <- order(clv); M <- z[ord, , drop = FALSE]; showy <- nrow(M) <= 50
  p <- plot_ly(x = o$timepoints, y = rownames(M), z = M, type = "heatmap", colorscale = colorscale %||% "RdBu", reversescale = TRUE, zmid = 0) %>%
    ply(xaxis = list(title = "Timepoint"), yaxis = list(title = "", showticklabels = showy, autorange = "reversed"),
        title = list(text = sprintf("%d features%s, rows ordered by program", nrow(M), if ((o$core_cut %||% 0) > 0) " (core)" else " (top by membership)"), x = .02, font = list(size = 12)))
  list(plot = p, mat = M, ttl = o$ttl, by = "Timepoint")
}
byod_temporal_enrich <- function(o, d, cluster, collection, cols = NULL) {
  if (!isTRUE(o$ok)) return(list(plot = emptyplot("run first"), summary = data.frame(), splits = character(0), ylab = "", ttl = "", cols = cols))
  cm <- .tp_members(o); genes <- cm$Feature[cm$cluster == cluster]
  plot_ora(run_ora(genes, collection, rownames(d$expr), d = d), 20,
           sprintf("%s ORA (%s, %d genes%s)", collection, cluster, length(genes), if ((o$core_cut %||% 0) > 0) " core" else ""), cols = cols)
}

.n_phrase <- function(p) {
  nu <- p[["n_units"]]; ul <- p[["unit"]]; np <- p[["n_profiles"]] %||% p[["n"]]
  if (is.null(nu) || !length(nu) || is.na(nu[1])) return(sprintf("n=%s profiles", np %||% "?"))
  if (identical(as.character(ul)[1], "profile")) return(sprintf("effective n=%s profiles", nu[1]))
  ns <- p[["n_subjects"]]; sl <- p[["subject_level"]]
  if (!is.null(ns) && length(ns) && !is.na(ns[1]) && !is.null(sl) && length(sl))
    return(sprintf("effective n=%s independent %s%s, containing %s paired %s%s across %s profiles",
                   nu[1], as.character(ul)[1], if (identical(as.integer(nu[1]), 1L)) "" else "s",
                   ns[1], as.character(sl)[1], if (identical(as.integer(ns[1]), 1L)) "" else "s",
                   np %||% "?"))
  sprintf("effective n=%s independent %s%s, from %s profiles", nu[1], as.character(ul)[1],
          if (identical(as.integer(nu[1]), 1L)) "" else "s", np %||% "?")
}

byod_design_text <- function(d) {
  if (is.null(d)) return("no study loaded")
  n_prof <- if (!is.null(d$meta)) nrow(d$meta) else NA_integer_
  subj <- d$subject %||% NULL
  nest <- d$nest %||% NULL
  if (is.null(subj))
    return(sprintf("%s profiles; no repeated-measures subject declared", n_prof))
  n_subj <- length(unique(as.character(d$meta[[subj]])))
  if (is.null(nest) || !nest %in% names(d$meta))
    return(sprintf("%d %ss, repeatedly sampled; %s profiles", n_subj, subj, n_prof))
  n_nest <- length(unique(as.character(d$meta[[nest]])))
  sprintf("%d %ss; %d %ss nested within them; repeated measurements by %s; %s profiles",
          n_nest, nest, n_subj, subj, subj, n_prof)
}

byod_methods_text <- function(records, study_name = NULL) {
  if (!length(records)) return(paste0(
    "Run one or more analyses (Differential, Enrichment, Pathways, WGCNA, DIABLO, Clinical ",
    "correlates, Group×time, Temporal programs, ...). A paste-ready Methods paragraph ",
    "will assemble here from the exact parameters you use."))
  v <- function(pkg) tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "NA")
  last <- list(); for (r in records) if (!is.null(r$module)) last[[r$module]] <- r$params

  declined <- Filter(function(p) identical(p$outcome, "declined"), last)
  last <- last[setdiff(names(last), names(declined))]
  g <- function(p, k, d = "?") { x <- p[[k]]; if (is.null(x) || !length(x) || (is.atomic(x) && is.na(x[1]))) d else as.character(x)[1] }
  study_name <- study_name[!is.na(study_name) & nzchar(study_name)]
  s <- sprintf("Interactive analyses were performed in TRACE-Omics (R %s)%s.", getRversion(),
               if (length(study_name) == 1L) sprintf(" on the %s dataset", study_name)
               else if (length(study_name) > 1L)
                 sprintf(" on %d separately declared datasets (%s), which were analysed independently and never merged",
                         length(study_name), paste(study_name, collapse = "; "))
               else "")

  if (length(study_name) > 1L)
    s <- c(s, "Each description below is of the most recent run of that module; the session log records which dataset each run used, and the two datasets were never pooled into a single model.")
  if (!is.null(p <- last[["DE (limma)"]])) {
    blk <- g(p, "block", "-"); pr <- tolower(g(p, "paired", ""))
    pair_phrase <- if (!blk %in% c("-", "NA", "?")) sprintf(", paired by %s", blk) else if (pr %in% c("true", "yes")) ", paired by subject" else ", unpaired"
    s <- c(s, sprintf("Differential expression used limma (v%s) moderated t-tests on the declared analysis scale%s (%s); the %s vs %s contrast was called significant at FDR < %s and absolute mean difference ≥ %s.",
      v("limma"), pair_phrase, .n_phrase(p), g(p, "g2"), g(p, "g1"), g(p, "fdr", "0.05"), g(p, "lfc", "1")))
  }
  if (!is.null(p <- last[["Enrichment"]]))
    s <- c(s, sprintf("Functional enrichment used %s against the %s collection.",
      if (identical(g(p, "mode"), "gsea")) sprintf("gene-set enrichment analysis (fgsea v%s)", v("fgsea")) else "hypergeometric over-representation",
      g(p, "collection")))
  annotated <- Filter(function(r) r$module %in% c("Enrichment", "Pathways", "cross-dataset pathway concordance") &&
    !identical(r$params$outcome, "declined"), records)
  if (length(annotated)) {
    notes <- unique(vapply(annotated, function(r) r$params$annotation %||% "", character(1)))
    s <- c(s, notes[nzchar(notes)])
  }
  if (!is.null(p <- last[["Pathways"]]))
    s <- c(s, if (requireNamespace("GSVA", quietly = TRUE))
      sprintf("Per-sample pathway activity was scored with %s (GSVA v%s) over the %s collection.", g(p, "method", "GSVA"), v("GSVA"), g(p, "collection"))
      else sprintf("Per-sample pathway activity was scored with a mean z-score fallback (the GSVA package was not available in this instance) over the %s collection.", g(p, "collection")))
  if (!is.null(p <- last[["WGCNA"]]))
    s <- c(s, if (requireNamespace("WGCNA", quietly = TRUE))
      sprintf("Co-expression modules were constructed with WGCNA (v%s) on the top %s most-variable features (soft-thresholding power %s).", v("WGCNA"), g(p, "ngenes", "2000"), g(p, "power", "auto-selected"))
      else sprintf("Co-expression modules were built with a correlation-based fallback (the WGCNA package was not available in this instance) on the top %s most-variable features.", g(p, "ngenes", "2000")))
  if (!is.null(p <- last[["DIABLO"]]))
    s <- c(s, sprintf("Cross-omics integration used a DIABLO multiblock sPLS-DA (mixOmics v%s) on the outcome %s (%s features per block), fitted on all matched samples. The fit itself is descriptive: it is not cross-validated and must not be read as a predictor.",
      v("mixOmics"), g(p, "out"), g(p, "nfeat")))

  if (!is.null(p <- last[["DIABLO cross-validation"]]))
    s <- c(s, sprintf("That model was then cross-validated with %s-fold cross-validation repeated %s times, holding out whole %ss (n=%s) rather than samples; the variance pre-filter, the centring and scaling, and the sPLS-DA feature selection were all computed inside each training fold. Balanced accuracy was %s.",
      g(p, "folds"), g(p, "repeats"), g(p, "unit"), g(p, "n_units"), g(p, "balanced_accuracy")))
  if (!is.null(last[["Integration"]]) || !is.null(last[["Concordance"]]))
    s <- c(s, "Cross-layer integration described the correlation between separately estimated mean differences. Per-gene concordance described correlations over matched profiles; repeated profiles were not treated as independent observations for significance testing.")
  if (!is.null(p <- last[["Subject variables: features"]]))
    s <- c(s, sprintf("Molecular correlates of clinical variables were computed at donor resolution (one value per subject) as per-feature %s correlations against %s (%s), with Benjamini–Hochberg FDR control.",
      sub("^(\\w)", "\\U\\1", g(p, "method", "Spearman"), perl = TRUE), g(p, "var"),
      if (!identical(g(p, "at_tp", "All"), "All")) sprintf("at the %s timepoint", g(p, "at_tp")) else "donor-averaged across timepoints"))
  if (!is.null(p <- last[["Group x time interaction"]])) {
    des <- if (tolower(g(p, "paired", "")) %in% c("true", "yes"))
             sprintf("as a within-%s paired contrast (fixed-block limma approximation to a (1|%s)+(1|subject) mixed model; n=%s matched %ss as the effective replication)",
                     g(p, "block", "block"), g(p, "block", "block"), g(p, "n_units", "?"), tolower(g(p, "block", "block")))
           else "with the subject modelled as a fixed nuisance block"
    s <- c(s, sprintf("Group-by-time interaction was tested with limma %s, identifying features whose temporal trajectory differed by %s.", des, g(p, "var")))
  }
  if (!is.null(p <- last[["Temporal programs"]]))
    s <- c(s, sprintf("Temporal programs were derived by limma trend testing across timepoints (%s design, subject-adjusted) and clustering the FDR-significant features' z-scored mean profiles by %s into %s programs.",
      g(p, "mode", "factor"), if (requireNamespace("e1071", quietly = TRUE)) "fuzzy c-means (e1071)" else "k-means", g(p, "k", "6")))
  cross_modules <- c("cross-dataset differential concordance", "cross-dataset pathway concordance")
  for (mod in intersect(cross_modules, names(last))) {
    p <- last[[mod]]
    sides <- vapply(c("A", "B"), function(side) sprintf(
      "%s: %s, %s, n=%s %s", g(p, paste0("study_", side)),
      g(p, paste0(side, "_contrast")), g(p, paste0(side, "_pairing")),
      g(p, paste0(side, "_n_units")), g(p, paste0(side, "_unit"))), character(1))
    s <- c(s, paste0("Cross-study contrasts were fitted separately with limma (", paste(sides, collapse = "; "), ")."))
    if (mod == cross_modules[1])
      s <- c(s, sprintf("Effect estimates were compared over %s shared features (%s), with FDR < %s assessed within the shared feature set in each study.",
        g(p, "shared_features"), g(p, "feature_rule"), g(p, "fdr")))
    else s <- c(s, sprintf("Pathway concordance compared %s scores for %s sets in the %s collection.",
      g(p, "scored_by"), g(p, "sets"), g(p, "collection")))
    if (isTRUE(suppressWarnings(as.numeric(g(p, "permutations", "0"))) > 0))
      s <- c(s, sprintf("Agreement was tested with %s permutations (%s; seed %s).",
        g(p, "permutations"), g(p, "permutation_scheme"), g(p, "permutation_seed")))
  }
  used_seeds <- unique(unlist(lapply(records, function(r) r$params$random_seed)))
  if (length(used_seeds)) s <- c(s, sprintf(
    "Stochastic analyses used the recorded random seed%s %s; the session log identifies the seed for each run.",
    if (length(used_seeds) == 1L) "" else "s", paste(used_seeds, collapse = ", ")))
  s <- c(s, "Figures (PNG/PDF) and result tables (CSV) were exported directly from the application.")

  templated <- c("DE (limma)", "Enrichment", "Pathways", "WGCNA", "DIABLO",
                 "DIABLO cross-validation", "Integration", "Concordance",
                 "Subject variables: features", "Group x time interaction",
                 "Temporal programs", cross_modules)
  extra <- setdiff(names(last), c(templated, "load-study"))
  for (mod in extra) {
    p <- last[[mod]]
    eng  <- g(p, "engine", NA_character_)
    unit <- g(p, "independent_units", NA_character_)
    bits <- character(0)
    if (!is.na(eng) && nzchar(eng) && eng != "?") bits <- c(bits, sprintf("using %s", eng))
    if (!is.na(unit) && nzchar(unit) && unit != "?") bits <- c(bits, sprintf("on %s", unit))
    s <- c(s, sprintf("%s was run%s%s.", mod,
                      if (length(bits)) paste0(" ", paste(bits, collapse = " ")) else "",
                      if (!identical(g(p, "declared_nest", "?"), "(none declared)") &&
                          !identical(g(p, "declared_nest", "?"), "?"))
                        sprintf(", with %s nested inside %s",
                                g(p, "declared_subject", "the subject"), g(p, "declared_nest", "?"))
                      else ""))
  }
  for (mod in names(declined))
    s <- c(s, sprintf("%s was requested and declined: %s.", mod,
                      g(declined[[mod]], "reason", "the resolved design does not support it")))
  paste(s, collapse = " ")
}

byod_module <- function(d, genes, splitby = "None", cols = NULL) {
  ylab <- "module score (mean z)"
  ms <- module_score(d, intersect(genes, rownames(d$expr)))
  if (is.null(ms)) return(c(byod_empty_tc(d, "No matching genes", ylab), list(n = 0)))
  td <- tc_data(d, setNames(ms$score, colnames(d$expr)), splitby)
  cc <- if (is.null(cols)) resolve_cols(td$splits, PAL_NAMES[1], if (identical(splitby, "None")) NULL else splitby) else cols

  cu <- collapse_to_resolved_unit(d, setNames(ms$score, colnames(d$expr)),
                                  keep = if (identical(splitby, "None")) NULL else splitby)
  tt <- byod_time_test(cu, d, split_col = splitby)
  ttl <- sprintf("%d genes · %s", ms$n, tt$text)
  list(plot = tc_plotly(td, ylab, ttl, cc), summary = td$summary, splits = td$splits,
       ylab = ylab, ttl = ttl, cols = cc, n = ms$n)
}

byod_build_traits <- function(d, samples) {
  tr <- data.frame(row.names = samples)
  tr[["Timepoint(ord)"]] <- as.numeric(factor(d$meta[samples, "Timepoint"], levels = d$tp_levels))
  if (!is.null(d$clinical_num)) for (v in d$clinical_num)
    tr[[v]] <- suppressWarnings(as.numeric(d$meta[samples, v]))
  if (!is.null(d$clinical_cat)) for (v in d$clinical_cat) {
    f <- factor(as.character(d$meta[samples, v]))
    if (nlevels(f) == 2) tr[[paste0(v, "=", levels(f)[2])]] <- as.numeric(f) - 1
  }
  tr[, colSums(!is.na(tr)) > 2, drop = FALSE]
}
cor_pmat <- function(r, n) { t <- r * sqrt((n - 2) / pmax(1e-12, 1 - r^2)); 2 * pt(-abs(t), n - 2) }

byod_coexpr <- function(d, nGenes = 2500L, power = NULL, minModuleSize = 20L, deepSplit = 2L,
                       seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)

  ucol_c    <- resolved_unit_col(d)
  n_units_c <- if (!is.null(ucol_c)) length(unique(as.character(d$meta[[ucol_c]])))
               else ncol(d$expr)
  if (n_units_c < MIN_UNITS$coexpression)
    return(list(ok = FALSE, msg = sprintf("A co-expression network needs at least %d independent units; this design has %d (from %d profiles).",
                                          MIN_UNITS$coexpression, n_units_c, ncol(d$expr))))
  rv   <- byod_rowvars(d$expr); keep <- head(order(rv, decreasing = TRUE), min(nGenes, sum(rv > 0)))
  X    <- t(d$expr[keep, , drop = FALSE])                 # samples x genes
  if (anyNA(X)) X <- apply(X, 2, function(v) { v[is.na(v)] <- mean(v, na.rm = TRUE); v })
  useW <- requireNamespace("WGCNA", quietly = TRUE) && requireNamespace("dynamicTreeCut", quietly = TRUE)
  if (useW) {
    utils::capture.output(suppressWarnings({              # keep pickSoftThreshold noise out of the logs
      if (is.null(power)) {
        sft <- try(WGCNA::pickSoftThreshold(X, powerVector = c(1:10, seq(12, 20, 2)),
                   networkType = "signed", verbose = 0), silent = TRUE)
        power <- if (!inherits(sft, "try-error") && !is.null(sft$powerEstimate) && !is.na(sft$powerEstimate))
                   sft$powerEstimate else 6L
      }
      adj  <- WGCNA::adjacency(X, power = power, type = "signed")
      diss <- 1 - WGCNA::TOMsimilarity(adj, TOMType = "signed", verbose = 0)
    }))
    tree <- if (requireNamespace("fastcluster", quietly = TRUE)) fastcluster::hclust(as.dist(diss), method = "average")
            else hclust(as.dist(diss), method = "average")
    lab  <- dynamicTreeCut::cutreeDynamic(dendro = tree, distM = diss, deepSplit = deepSplit,
                                 pamRespectsDendro = FALSE, minClusterSize = minModuleSize, verbose = 0)
    colors <- WGCNA::labels2colors(lab)
    ME <- WGCNA::moduleEigengenes(X, colors)$eigengenes
    names(ME) <- sub("^ME", "", names(ME))
  } else {                                                # base fallback
    if (is.null(power)) power <- 6L
    Cor  <- cor(X, use = "pairwise.complete.obs"); adj <- ((1 + Cor) / 2)^power
    tree <- hclust(as.dist(1 - adj), method = "average")
    k    <- max(2L, min(14L, floor(ncol(X) / (minModuleSize * 2))))
    cl   <- cutree(tree, k = k); tb <- table(cl)
    small <- as.integer(names(tb)[tb < minModuleSize]); cl[cl %in% small] <- 0L
    colors <- ifelse(cl == 0L, "grey", paste0("M", match(cl, sort(unique(cl[cl > 0])))))
    mods <- setdiff(unique(colors), "grey")
    ME <- as.data.frame(lapply(mods, function(m) {
      g <- which(colors == m); pcx <- prcomp(X[, g, drop = FALSE], scale. = TRUE)$x[, 1]
      if (suppressWarnings(cor(pcx, rowMeans(X[, g, drop = FALSE]))) < 0) pcx <- -pcx; pcx }))
    names(ME) <- mods; rownames(ME) <- rownames(X)
  }
  keepM  <- setdiff(colnames(ME), c("grey", "0")); ME <- ME[, keepM, drop = FALSE]
  traits <- byod_build_traits(d, rownames(X))

  subj_col <- ucol_c %||% d$subject
  subjW <- if (!is.null(subj_col) && subj_col %in% names(d$meta)) as.character(d$meta[rownames(X), subj_col]) else NULL
  MT <- MTp <- MTn <- matrix(NA_real_, nrow = ncol(ME), ncol = ncol(traits), dimnames = list(colnames(ME), colnames(traits)))
  if (ncol(ME) && ncol(traits)) for (tj in colnames(traits)) {
    tvec <- traits[[tj]]
    bs <- !is.null(subjW) && all(tapply(tvec, subjW, function(z) length(unique(z[!is.na(z)])) <= 1))
    if (bs) { du <- unique(subjW)
      MEc <- vapply(colnames(ME), function(m) tapply(ME[[m]], subjW, mean, na.rm = TRUE)[du], numeric(length(du)))
      tc  <- vapply(du, function(u) { z <- tvec[subjW == u]; z <- z[!is.na(z)]; if (length(z)) z[1] else NA_real_ }, numeric(1))
    } else { MEc <- as.matrix(ME); tc <- tvec }
    for (m in colnames(ME)) {
      ok <- is.finite(MEc[, m]) & is.finite(tc); ne <- sum(ok); MTn[m, tj] <- ne
      r <- if (ne >= 4) suppressWarnings(cor(MEc[ok, m], tc[ok])) else NA_real_; MT[m, tj] <- r
      MTp[m, tj] <- if (is.finite(r) && abs(r) < 1) { tval <- r * sqrt((ne - 2) / (1 - r^2)); 2 * pt(-abs(tval), ne - 2) } else NA_real_
    }
  }
  MTq <- MTp; if (any(is.finite(MTp))) MTq[] <- p.adjust(MTp, "BH")   # BH over the full module×trait family
  kME <- suppressWarnings(cor(X, ME, use = "pairwise.complete.obs"))    # gene x module
  sizes <- table(factor(colors, levels = c(keepM, "grey")))
  hub <- sapply(colnames(ME), function(m) { gi <- which(colors == m)
    if (!length(gi)) return(NA_character_); colnames(X)[gi][which.max(kME[gi, m])] })
  list(colors = setNames(colors, colnames(X)), ME = ME, MT = MT, MTp = MTp, MTq = MTq, MTn = MTn, kME = kME,
       hub = hub, sizes = sizes, power = power,

       engine = if (useW) "WGCNA package" else "base correlation fallback",
       executed_method = if (useW) "WGCNA" else "base signed correlation",
       stepped_down = !useW,
       unit = ucol_c %||% "profile", n_units = n_units_c, n_profiles = ncol(d$expr),
       nGenes = ncol(X), genes = colnames(X))
}

byod_mt_height <- function(obj, per_row = 26L, min_h = 260L, max_h = 900L) {
  n <- if (length(obj$MT)) nrow(obj$MT) else 0L
  as.integer(min(max_h, max(min_h, 130L + per_row * n)))
}

byod_mt_heatmap_gg <- function(obj) {
  MT <- obj$MT
  if (is.null(MT) || !length(MT) || !nrow(MT))
    return(emptyplot_gg("No module-trait correlation to show."))
  q <- obj$MTq %||% obj$MTp
  sizes <- obj$sizes
  df <- expand.grid(module = rownames(MT), trait = colnames(MT),
                    stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  df$r <- as.vector(MT)
  df$q <- if (is.null(q)) NA_real_ else as.vector(q)
  df$n <- if (is.null(obj$MTn)) NA_integer_ else as.vector(obj$MTn)
  star <- function(z) ifelse(is.na(z), "", ifelse(z < .001, "***",
                      ifelse(z < .01, "**", ifelse(z < .05, "*", ""))))

  dense <- nrow(MT) > 20L || ncol(MT) > 18L
  df$lab <- if (dense) star(df$q)
            else ifelse(is.finite(df$r), paste0(sprintf("%.2f", df$r), star(df$q)), "")

  lab_of <- function(m) {
    n <- suppressWarnings(as.integer(sizes[m]))
    ifelse(is.na(n), m, sprintf("%s  (%s genes)", m, format(n, big.mark = ",")))
  }
  df$module <- factor(lab_of(df$module), levels = rev(lab_of(rownames(MT))))
  ggplot(df, aes(trait, module, fill = r)) +
    geom_tile(colour = "white", linewidth = .7) +
    geom_text(aes(label = lab), size = if (dense) 3.2 else 3, colour = "grey10",
              vjust = if (dense) 0.72 else 0.5) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-1, 1),
                         name = expression(italic(r))) +
    labs(title = sprintf("Module-trait correlation (%s, soft power %s)",
                         obj$engine %||% "WGCNA", obj$power %||% "?"),
         subtitle = paste(strwrap(sprintf(
           "%d modules from %s genes. Modules are named by colour, which is WGCNA's convention and carries no meaning; grey is reserved for genes in no module. Stars are BH-adjusted: * < 0.05, ** < 0.01, *** < 0.001, at n = %s.",
           ncol(obj$ME), format(obj$nGenes %||% NA, big.mark = ","),
           if (all(is.na(df$n))) "the resolved unit" else format(max(df$n, na.rm = TRUE), big.mark = ",")),
           width = 92), collapse = "\n"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 7.5, colour = "grey30", lineheight = 1.15),
          axis.text.x   = element_text(angle = 35, hjust = 1, size = 8.5),
          axis.text.y   = element_text(size = 8.5),
          panel.grid    = element_blank(),
          legend.key.height = unit(.5, "cm"))
}

byod_mt_heatmap <- function(obj, colorscale = NULL) {
  MT <- obj$MT; if (!length(MT) || !nrow(MT)) return(emptyplot("No modules/traits to correlate"))
  q <- obj$MTq %||% obj$MTp; q[is.na(q)] <- 1        # stars use BH-adjusted p at the correct n
  star <- ifelse(q < .001, "***", ifelse(q < .01, "**", ifelse(q < .05, "*", "")))

  dense <- ncol(MT) > 18L || nrow(MT) > 20L
  lab <- if (dense) star else matrix(paste0(sprintf("%.2f", MT), star), nrow = nrow(MT))
  fsz <- if (dense) 11 else 10

  vertical <- ncol(MT) > 12L
  trunc_trait <- function(x, n = 26L)
    ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1L), "\u2026"), x)
  wrap_trait <- function(x, width = 14L)
    vapply(x, function(z) if (nchar(z) <= width) z
           else paste(strwrap(z, width = width), collapse = "<br>"),
           character(1), USE.NAMES = FALSE)
  xlab <- if (vertical) trunc_trait(colnames(MT)) else wrap_trait(colnames(MT))

  cs <- colorscale %||% list(c(0, HUMAN), c(.5, "#FFFFFF"), c(1, PIG))
  p <- plot_ly(x = xlab, y = rownames(MT), z = MT, type = "heatmap", colorscale = cs, zmid = 0,
      zmin = -1, zmax = 1, xgap = 1, ygap = 1,
      text = matrix(sprintf("%.3f", MT), nrow = nrow(MT)),
      hovertemplate = "%{y} ~ %{x}<br>r = %{text}<extra></extra>",
      colorbar = list(title = list(text = "r", side = "top"), len = .7, thickness = 12))
  if (any(nzchar(lab)))
    p <- add_annotations(p, x = rep(xlab, each = nrow(MT)), y = rep(rownames(MT), times = ncol(MT)),
      text = as.vector(lab), showarrow = FALSE, font = list(size = fsz, color = "#222"))
  ply(p,
      margin = list(l = 90, r = 20, t = 62, b = if (vertical) 190 else 150),
      xaxis = list(title = "", tickangle = if (vertical) -90 else -45,
                   tickfont = list(size = if (vertical) 9.5 else 10),
                   automargin = TRUE, side = "bottom"),
      yaxis = list(title = "module", tickfont = list(size = 10), automargin = TRUE),
      title = list(text = paste0("Module-trait correlation · ", obj$engine, " · power ", obj$power,
                                 if (dense) " · * = BH<.05 (hover for r)" else " · * = BH<.05",
                                 " · clinical traits at donor level"),
                   x = .02, font = list(size = 12.5)))
}

byod_module_summary <- function(obj) {
  mods <- colnames(obj$ME); if (!length(mods)) return(data.frame())
  do.call(rbind, lapply(mods, function(m) {
    if (length(obj$MT) && m %in% rownames(obj$MT)) {
      j <- which.max(abs(obj$MT[m, ])); tr <- colnames(obj$MT)[j]
      r <- round(obj$MT[m, j], 2); p <- signif(obj$MTp[m, j], 2)
    } else { tr <- NA; r <- NA; p <- NA }
    data.frame(module = m, size = as.integer(obj$sizes[m]), hub_gene = obj$hub[[m]],
               top_trait = tr, r = r, p = p, stringsAsFactors = FALSE)
  }))
}

byod_dist_plot <- function(d) {
  q <- apply(d$expr, 2, quantile, probs = c(.25, .5, .75), na.rm = TRUE)
  samp <- colnames(d$expr); tp <- factor(d$meta[samp, "Timepoint"], levels = d$tp_levels)
  df <- data.frame(sample = samp, tp = tp, lo = q[1, ], mid = q[2, ], hi = q[3, ], stringsAsFactors = FALSE)
  df <- df[order(df$tp, df$sample), ]; df$idx <- seq_len(nrow(df))
  lv <- levels(droplevels(df$tp)); cols <- setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = "Timepoint")), lv)
  p <- plot_ly()
  for (L in lv) { s <- df[df$tp == L, ]
    p <- add_trace(p, data = s, x = ~idx, y = ~mid, type = "scatter", mode = "markers", name = L,
      marker = list(color = cols[[L]], size = 6),
      error_y = list(array = s$hi - s$mid, arrayminus = s$mid - s$lo, color = cols[[L]], width = 0),
      text = ~paste0(sample, "<br>median=", round(mid, 2)), hoverinfo = "text") }
  p %>% ply(xaxis = list(title = "sample (ordered by timepoint)"),
    yaxis = list(title = paste0("per-sample IQR (", d$unit, ")")), legend = list(orientation = "h", y = 1.12),
    title = list(text = "Sample distributions — aligned medians/IQR indicate good normalization", x = .02, font = list(size = 12)))
}

byod_pathway_activity <- function(d, collection = "Hallmark", method = "GSVA", by = "Timepoint", maxpath = 40L) {
  sets <- annotation_sets(d, collection)
  expr <- as.matrix(d$expr[intersect(rownames(d$expr), annotation_universe(d, collection)), , drop = FALSE])
  gs <- tryCatch({
    par <- if (identical(method, "ssGSEA")) GSVA::ssgseaParam(expr, sets, minSize = 5, maxSize = 500)
           else GSVA::gsvaParam(expr, sets, minSize = 5, maxSize = 500)
    suppressMessages(GSVA::gsva(par))
  }, error = function(e) NULL)
  used <- method
  gsva_err <- NULL
  if (is.null(gs)) {
    gsva_err <- "the GSVA package is unavailable or failed on this matrix"
    Z <- t(scale(t(expr))); Z[!is.finite(Z)] <- 0
    ok <- names(sets)[vapply(names(sets), function(p) length(intersect(sets[[p]], rownames(Z))) >= 3, logical(1))]

    if (!length(ok)) {
      hit <- sum(toupper(rownames(expr)) %in% toupper(unique(unlist(sets))))
      return(list(mat = matrix(numeric(0), 0, 0), group_mat = matrix(numeric(0), 0, 0),
                  by = by, method = used, requested_method = method, executed_method = used,
                  stepped_down = FALSE, collection = collection, maxpath = maxpath,
                  ok = FALSE,
                  msg = sprintf("No %s gene set has three or more eligible features: %s of %s match the collection. Check the organism, identifier type and mapping coverage for %s, or supply a native-ID GMT.",
                                collection, format(hit, big.mark = ","),
                                format(nrow(expr), big.mark = ","), collection)))
    }
    gs <- do.call(rbind, lapply(ok, function(p) colMeans(Z[intersect(sets[[p]], rownames(Z)), , drop = FALSE])))
    rownames(gs) <- ok; used <- "mean z-score"
  }
  if (is.null(dim(gs)) || !nrow(gs))
    return(list(mat = matrix(numeric(0), 0, 0), group_mat = matrix(numeric(0), 0, 0),
                by = by, method = used, requested_method = method, executed_method = used,
                stepped_down = FALSE, collection = collection, maxpath = maxpath, ok = FALSE,
                msg = sprintf("%s returned no scoreable gene set for this matrix.", method)))
  rownames(gs) <- gsub("_", " ", rownames(gs)); colnames(gs) <- colnames(expr)
  gf <- droplevels(grp_factor(d, by)); lv <- levels(gf)
  M <- sapply(lv, function(L) rowMeans(gs[, which(gf == L), drop = FALSE], na.rm = TRUE))
  if (is.null(dim(M))) M <- matrix(M, nrow = nrow(gs), dimnames = list(rownames(gs), lv))
  ord <- order(apply(M, 1, function(x) var(x, na.rm = TRUE)), decreasing = TRUE)
  list(mat = gs, group_mat = M[ord, , drop = FALSE], by = by, method = used,
       requested_method = method, executed_method = used,
       stepped_down = !identical(used, method), step_down_reason = gsva_err,
       collection = collection, maxpath = maxpath, ok = TRUE)
}
byod_pathway_heatmap <- function(obj, colorscale = NULL) {
  M <- utils::head(obj$group_mat, obj$maxpath)
  if (!nrow(M)) return(list(plot = emptyplot(obj$msg %||% "no pathways could be scored"),
                            mat = NULL, ttl = "", by = obj$by))
  Z <- t(scale(t(M))); Z[!is.finite(Z)] <- 0
  if (nrow(Z) > 2) Z <- Z[hclust(dist(Z))$order, , drop = FALSE]
  cs <- colorscale %||% list(c(0, HUMAN), c(.5, "#FFFFFF"), c(1, PIG))
  ttl <- paste0(obj$collection, " activity (", obj$method, ", row z)")
  p <- plot_ly(x = colnames(Z), y = rownames(Z), z = Z, type = "heatmap", colorscale = cs, zmid = 0,
      hovertemplate = "%{y} @ %{x}<br>z=%{z:.2f}<extra></extra>") %>%
    ply(xaxis = list(title = obj$by), yaxis = list(title = "", tickfont = list(size = 9)), margin = list(l = 290),
        title = list(text = ttl, x = .02, font = list(size = 13)))
  list(plot = p, mat = Z, ttl = ttl, by = obj$by)
}
byod_pathway_tc <- function(d, obj, pathway, splitby = "None", cols = NULL) {
  ri <- which(rownames(obj$mat) == pathway); ylab <- paste0(obj$method, " score")
  if (!length(ri)) return(byod_empty_tc(d, "Pick a pathway", ylab))
  val <- setNames(as.numeric(obj$mat[ri[1], ]), colnames(obj$mat))
  td <- tc_data(d, val, splitby)
  cc <- if (is.null(cols)) resolve_cols(td$splits, PAL_NAMES[1], if (identical(splitby, "None")) NULL else splitby) else cols
  list(plot = tc_plotly(td, ylab, pathway, cc), summary = td$summary, splits = td$splits,
       ylab = ylab, ttl = pathway, cols = cc)
}

byod_concordance <- function(d1, d2, method = "spearman", shared_samples = NULL, shared_genes = NULL) {
  S <- shared_samples %||% intersect(colnames(d1$expr), colnames(d2$expr))
  G <- shared_genes  %||% intersect(rownames(d1$expr), rownames(d2$expr))

  if (length(S) < 4 || length(G) < 10) {
    why <- if (length(S) < 4 && length(G) < 10)
             sprintf("only %d matched profile(s) and %d shared feature(s)", length(S), length(G))
           else if (length(S) < 4)
             sprintf("only %d profile(s) appear in both layers (needs 4); check that the two matrices use the same sample IDs", length(S))
           else
             sprintf("the two layers share %d feature name(s) (needs 10), although %d profiles match; the layers are keyed on different identifier namespaces, so no feature can be paired",
                     length(G), length(S))
    return(list(tab = data.frame(), n = 0, nS = length(S), nG = length(G), method = method,
                ok = FALSE, msg = paste0("Cross-layer concordance needs both layers on the same samples and the same feature names: ", why, ".")))
  }
  M <- d1$expr[G, S, drop = FALSE]; P <- d2$expr[G, S, drop = FALSE]
  rho <- vapply(seq_len(nrow(M)), function(i)
    suppressWarnings(cor(M[i, ], P[i, ], method = method, use = "pairwise.complete.obs")), numeric(1))
  cls <- ifelse(is.na(rho), "n/a", ifelse(rho >= 0.4, "concordant",
         ifelse(rho <= 0, "discordant", "weak positive")))
  tab <- data.frame(gene = G, rho = round(rho, 3),
                    mRNA_mean = round(rowMeans(M, na.rm = TRUE), 2),
                    protein_mean = round(rowMeans(P, na.rm = TRUE), 2),
                    class = cls, stringsAsFactors = FALSE)
  tab <- tab[order(tab$rho, decreasing = TRUE), ]
  list(tab = tab, n = length(G), nS = length(S), method = method, ok = TRUE,
       med_rho = median(rho, na.rm = TRUE), frac_disc = mean(cls == "discordant", na.rm = TRUE),
       discordant_genes = tab$gene[tab$class == "discordant"], universe = G)
}
BYOD_CONC_COLS <- c(concordant = "#1E3A5F", "weak positive" = "#E69F00",
                    discordant = "#B0413E", "n/a" = "#CCCCCC")
byod_conc_hist <- function(obj) {
  tab <- obj$tab; if (!nrow(tab)) return(emptyplot(obj$msg %||% "Need matched samples"))
  tab$class <- factor(tab$class, levels = names(BYOD_CONC_COLS))
  plot_ly(tab, x = ~rho, color = ~class, colors = BYOD_CONC_COLS, type = "histogram", nbinsx = 40) %>%
    ply(barmode = "stack", xaxis = list(title = paste0("per-gene mRNA-protein ", obj$method, " r")),
        yaxis = list(title = "genes"), legend = list(orientation = "h", y = 1.12),
        title = list(text = sprintf("%d shared genes across %d matched samples · median r = %.2f · %.0f%% discordant",
                     obj$n, obj$nS, obj$med_rho, 100 * obj$frac_disc), x = .02, font = list(size = 12)))
}
byod_conc_scatter <- function(d1, d2, gene, S = NULL, lab1 = "mRNA", lab2 = "Protein") {
  S <- S %||% intersect(colnames(d1$expr), colnames(d2$expr))
  if (!(gene %in% rownames(d1$expr)) || !(gene %in% rownames(d2$expr)))
    return(emptyplot(paste(gene, "not measured in both layers")))
  x <- as.numeric(d1$expr[gene, S]); y <- as.numeric(d2$expr[gene, S])
  tp <- factor(d1$meta[S, "Timepoint"], levels = d1$tp_levels)
  rho <- suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
  lv <- levels(droplevels(tp)); cc <- setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = "Timepoint")), lv)
  plot_ly(x = x, y = y, color = droplevels(tp), colors = unname(cc[lv]), type = "scatter", mode = "markers",
      marker = list(size = 8, line = list(width = .4, color = "white")),
      text = paste0(S, "<br>", tp), hoverinfo = "text") %>%
    ply(xaxis = list(title = paste0(gene, " — ", lab1)), yaxis = list(title = paste0(gene, " — ", lab2)),
        legend = list(orientation = "h", y = 1.12),
        title = list(text = sprintf("%s · descriptive Spearman rho = %.2f (%d matched profiles)",
                     gene, rho, sum(is.finite(x) & is.finite(y))), x = .02, font = list(size = 13)))
}

byod_gene_gam <- function(d, gene, splitby = "None", cols = NULL) {
  v <- expr_vec(d, gene); if (is.null(v)) return(emptyplot(paste(gene, "not found")))
  tp <- factor(d$meta$Timepoint, levels = d$tp_levels); tnum <- as.numeric(tp)
  use_split <- !identical(splitby, "None") && !is.null(d$meta[[splitby]])
  grp <- if (use_split) as.character(droplevels(grp_factor(d, splitby))) else rep("all", length(v))

  ucol <- resolved_unit_col(d, keep = if (use_split) splitby else NULL)
  scol <- d$subject %||% NULL
  uvec <- if (!is.null(ucol)) as.character(d$meta[[ucol]]) else rep(NA_character_, length(v))
  svec <- if (!is.null(scol) && !identical(scol, ucol) && scol %in% names(d$meta))
            as.character(d$meta[[scol]]) else rep(NA_character_, length(v))
  df <- data.frame(t = tnum, y = v, grp = grp, unit = uvec, sub = svec, stringsAsFactors = FALSE)
  df <- df[is.finite(df$t) & is.finite(df$y) & !is.na(df$grp), , drop = FALSE]
  if (!nrow(df)) return(emptyplot("No data"))
  grps <- unique(df$grp)
  cc <- if (is.null(cols)) resolve_cols(grps, PAL_NAMES[1], if (use_split) splitby else NULL) else cols
  rgba <- function(hex, a) { r <- grDevices::col2rgb(hex); sprintf("rgba(%d,%d,%d,%.2f)", r[1], r[2], r[3], a) }
  xs <- seq(min(df$t), max(df$t), length.out = 100); p <- plot_ly()
  n_re <- 0L; n_sub <- 0L                     # levels that entered a random effect
  for (g in grps) {
    s <- df[df$grp == g, ]; col <- if (g %in% names(cc)) cc[[g]] else MUTED
    nu <- length(unique(s$t))
    pts <- data.frame(t = s$t, y = s$y)
    p <- add_trace(p, data = pts, x = ~t, y = ~y, type = "scatter", mode = "markers", name = g, legendgroup = g,
                   marker = list(color = rgba(col, .45), size = 6), showlegend = FALSE, hoverinfo = "skip")
    s$unit <- factor(s$unit)
    has_re  <- nlevels(s$unit) >= 3L
    s$sub   <- factor(s$sub)
    has_sub <- has_re && !all(is.na(s$sub)) && nlevels(s$sub) > nlevels(s$unit)
    if (has_re) n_re <- max(n_re, nlevels(s$unit))
    if (has_sub) n_sub <- max(n_sub, nlevels(s$sub))
    fit <- if (nu >= 4 && nrow(s) >= 6)
      tryCatch(if (has_sub) mgcv::gam(y ~ s(t, k = min(nu, 5)) + s(unit, bs = "re") + s(sub, bs = "re"), data = s)
               else if (has_re) mgcv::gam(y ~ s(t, k = min(nu, 5)) + s(unit, bs = "re"), data = s)
               else        mgcv::gam(y ~ s(t, k = min(nu, 5)), data = s),
               error = function(e) NULL) else NULL
    if (!is.null(fit)) {

      nd <- data.frame(t = xs)
      if (has_re)  nd$unit <- s$unit[1]
      if (has_sub) nd$sub  <- s$sub[1]
      pr <- predict(fit, newdata = nd, se.fit = TRUE,
                    exclude = c(if (has_re) "s(unit)", if (has_sub) "s(sub)"))
      band <- data.frame(x = xs, fit = as.numeric(pr$fit),
                         lo = as.numeric(pr$fit - 1.96 * pr$se.fit), hi = as.numeric(pr$fit + 1.96 * pr$se.fit))
      if (has_re)
        p <- add_ribbons(p, data = band, x = ~x, ymin = ~lo, ymax = ~hi, line = list(width = 0),
                         fillcolor = rgba(col, .18), name = g, legendgroup = g, showlegend = FALSE, hoverinfo = "skip")
      p <- add_trace(p, data = band, x = ~x, y = ~fit, type = "scatter", mode = "lines",
                     line = list(color = col, width = 3), name = g, legendgroup = g)
    } else {                                     # too few timepoints for a smooth: connect means
      mu <- tapply(s$y, s$t, mean, na.rm = TRUE)
      md <- data.frame(t = as.numeric(names(mu)), y = as.numeric(mu))
      p <- add_trace(p, data = md, x = ~t, y = ~y, type = "scatter", mode = "lines+markers",
                     line = list(color = col, width = 3), marker = list(color = col, size = 7), name = g, legendgroup = g)
    }
  }
  p %>% ply(xaxis = list(tickvals = seq_along(d$tp_levels), ticktext = d$tp_levels, title = ""),
    yaxis = list(title = d$unit), showlegend = (length(grps) > 1), legend = list(orientation = "h", y = 1.12),
    title = list(text = paste0(gene, if (n_re >= 3L)
        sprintf(" \u00b7 GAM smooth, 95%% band adjusted for %d %s%s", n_re, ucol,
                if (n_sub >= 3L) sprintf(" and %d %s within them", n_sub, scol) else "")
      else " \u00b7 GAM smooth, descriptive (too few independent units for an interval)"),
      x = .02, font = list(size = 13)))
}

# Full-data DIABLO fit for descriptive projection; validation is a separate operation.
byod_diablo <- function(study, l1, l2, outcome = "Timepoint", nfeat = 500L, ncomp = 2L, keepx = 25L,
                       seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  S <- study$matched_samples
  if (length(S) < 8) stop("DIABLO needs >= 8 sample-matched observations.")
  d1 <- study$layers[[l1]]; d2 <- study$layers[[l2]]
  Y <- droplevels(grp_factor(d1, outcome)[match(S, rownames(d1$meta))])
  keep <- !is.na(Y); S <- S[keep]; Y <- droplevels(Y[keep])
  if (nlevels(Y) < 2) stop("Outcome needs >= 2 groups among matched samples.")
  if (min(table(Y)) < 3) stop("Each outcome group needs >= 3 matched samples.")
  topv <- function(d) { rv <- byod_rowvars(d$expr[, S, drop = FALSE])
    names(sort(rv, decreasing = TRUE))[seq_len(min(nfeat, sum(rv > 0)))] }
  X <- list(); X[[l1]] <- t(d1$expr[topv(d1), S, drop = FALSE]); X[[l2]] <- t(d2$expr[topv(d2), S, drop = FALSE])
  ncomp <- max(1L, min(as.integer(ncomp), nlevels(Y) - 1L, 2L))
  design <- matrix(0.5, 2, 2); diag(design) <- 0; dimnames(design) <- list(c(l1, l2), c(l1, l2))
  keepX <- list(); keepX[[l1]] <- rep(keepx, ncomp); keepX[[l2]] <- rep(keepx, ncomp)
  model <- mixOmics::block.splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design)
  list(model = model, blocks = c(l1, l2), outcome = outcome, Y = Y, ncomp = ncomp, n = length(S), seed = seed)
}

byod_cv_complete <- function(pred, samples, repeats) {
  expected <- expand.grid(sample = samples, repeat_i = seq_len(repeats), stringsAsFactors = FALSE)
  key <- function(d) paste(d$sample, d$repeat_i, sep = "\r")
  nrow(pred) == nrow(expected) && !anyDuplicated(key(pred)) && setequal(key(pred), key(expected))
}

# Keep each independent unit in one fold; preprocess and select features within training folds.
byod_diablo_cv <- function(study, l1, l2, outcome = "Timepoint", nfeat = 500L,
                           ncomp = 2L, keepx = 25L, group_col = NULL,
                           folds = 5L, repeats = 5L, seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  S  <- study$matched_samples
  d1 <- study$layers[[l1]]; d2 <- study$layers[[l2]]
  if (length(S) < 8) return(list(ok = FALSE, msg = "Cross-validation needs >= 8 sample-matched observations."))
  Y <- droplevels(grp_factor(d1, outcome)[match(S, rownames(d1$meta))])
  keep <- !is.na(Y); S <- S[keep]; Y <- droplevels(Y[keep])
  if (nlevels(Y) < 2) return(list(ok = FALSE, msg = "Outcome needs >= 2 groups among matched samples."))

  if (is.null(group_col)) group_col <- d1$nest %||% d1$subject
  if (is.null(group_col) || !group_col %in% names(d1$meta)) {
    grp <- S; group_col <- "sample"
  } else {
    grp <- as.character(d1$meta[[group_col]][match(S, rownames(d1$meta))])
  }
  units <- unique(grp)
  if (length(units) < folds)
    return(list(ok = FALSE, msg = sprintf("Only %d independent %ss; %d-fold grouped cross-validation needs at least %d.",
                                          length(units), group_col, folds, folds)))

  X1 <- d1$expr[, S, drop = FALSE]; X2 <- d2$expr[, S, drop = FALSE]
  design <- matrix(0.5, 2, 2); diag(design) <- 0
  dimnames(design) <- list(c(l1, l2), c(l1, l2))
  ncomp <- max(1L, min(as.integer(ncomp), nlevels(Y) - 1L, 2L))

  assignments <- lapply(seq_len(repeats), function(i)
    setNames(sample(rep_len(seq_len(folds), length(units))), units))
  rows <- list(); per_rep <- numeric(0)   # per-repeat BALANCED accuracy
  for (rep_i in seq_len(repeats)) {
    assign <- assignments[[rep_i]]
    rep_obs <- character(0); rep_pred <- character(0)
    for (k in seq_len(folds)) {
      te <- which(assign[grp] == k); tr <- setdiff(seq_along(S), te)
      if (!length(te) || !length(tr)) next
      ytr <- droplevels(Y[tr])
      if (nlevels(ytr) < 2 || min(table(ytr)) < 2) next
      te <- te[as.character(Y[te]) %in% levels(ytr)]
      if (!length(te)) next

      topv <- function(M) { rv <- byod_rowvars(M[, tr, drop = FALSE])
        names(sort(rv, decreasing = TRUE))[seq_len(min(nfeat, sum(rv > 0)))] }
      f1 <- topv(X1); f2 <- topv(X2)
      Xtr <- setNames(list(t(X1[f1, tr, drop = FALSE]), t(X2[f2, tr, drop = FALSE])), c(l1, l2))
      Xte <- setNames(list(t(X1[f1, te, drop = FALSE]), t(X2[f2, te, drop = FALSE])), c(l1, l2))
      kx  <- setNames(list(rep(min(keepx, length(f1)), ncomp),
                           rep(min(keepx, length(f2)), ncomp)), c(l1, l2))
      fit <- tryCatch(mixOmics::block.splsda(X = Xtr, Y = ytr, ncomp = ncomp,
                                             keepX = kx, design = design),
                      error = function(e) NULL)
      if (is.null(fit)) next
      pr <- tryCatch(stats::predict(fit, newdata = Xte), error = function(e) NULL)
      if (is.null(pr)) next
      cls <- as.character(pr$WeightedVote$max.dist[, ncomp])
      rows[[length(rows) + 1L]] <- data.frame(
        repeat_i = rep_i, fold = k, sample = S[te],
        unit = grp[te], observed = as.character(Y[te]), predicted = cls,
        stringsAsFactors = FALSE)
      rep_obs <- c(rep_obs, as.character(Y[te])); rep_pred <- c(rep_pred, cls)
    }
    if (length(rep_obs)) {
      cmr <- table(factor(rep_obs, levels(Y)), factor(rep_pred, levels(Y)))
      per_rep <- c(per_rep, mean(diag(cmr) / pmax(1, rowSums(cmr)), na.rm = TRUE))
    }
  }
  if (!length(rows))
    return(list(ok = FALSE, msg = "No cross-validation fold was estimable on this design."))

  pred <- do.call(rbind, rows)
  lv   <- levels(Y)
  cm   <- table(observed = factor(pred$observed, lv), predicted = factor(pred$predicted, lv))

  per_class <- diag(cm) / pmax(1, rowSums(cm))
  bal  <- mean(per_class, na.rm = TRUE)
  acc  <- mean(pred$observed == pred$predicted, na.rm = TRUE)
  complete <- byod_cv_complete(pred, S, repeats)
  list(ok = complete, msg = if (!complete) "Cross-validation is incomplete: at least one held-out sample lacks a prediction. Inspect the partial predictions and revise the design or model settings." else NULL,
       pred = pred, confusion = cm, per_class = per_class,
       accuracy = acc, balanced_accuracy = bal, balanced_error_rate = 1 - bal,
       chance = 1 / nlevels(Y), per_repeat_balanced_accuracy = per_rep,
       repeat_sd = if (length(per_rep) > 1) stats::sd(per_rep) else NA_real_,

       complete = byod_cv_complete(pred, S, repeats),
       n = length(S), n_units = length(units), unit = group_col,
       folds = folds, repeats = repeats, ncomp = ncomp, nfeat = nfeat, keepx = keepx, seed = seed,
       note = sprintf(paste("%d-fold cross-validation repeated %d times, holding out whole %ss:",
                            "balanced accuracy %.1f%% against %.1f%% for chance (%d %ss, %d matched samples).",
                            "The variance pre-filter, the scaling and the feature selection are all computed",
                            "inside each training fold."),
                      folds, repeats, group_col, 100 * bal, 100 / nlevels(Y),
                      length(units), group_col, length(S)))
}

byod_diablo_plot <- function(obj, which = "circos") {
  restore_rng <- trace_rng_scope(obj$seed %||% trace_seed()); on.exit(restore_rng(), add = TRUE)
  m <- obj$model
  if (which == "circos" && obj$ncomp < 1) which <- "sample"
  r <- tryCatch(switch(which,
    "sample"  = mixOmics::plotIndiv(m, ind.names = FALSE, legend = TRUE, ellipse = TRUE,
                                    title = paste0("DIABLO samples — ", obj$outcome)),
    "var"     = mixOmics::plotVar(m, var.names = FALSE, legend = TRUE, cutoff = 0.5,
                                  title = "DIABLO correlation circle"),
    "circos"  = mixOmics::circosPlot(m, cutoff = 0.7, size.variables = 0.6, size.labels = 1.1, line = TRUE),
    "cim"     = mixOmics::cimDiablo(m, margins = c(8, 14), size.legend = 0.7,
                                    color.blocks = c("#1E3A5F", "#2C8C99")),
    "network" = mixOmics::network(m, cutoff = 0.6, color.node = c("#1E3A5F", "#2C8C99")),
    mixOmics::plotIndiv(m)), error = function(e) e)
  if (inherits(r, "error")) { grid::grid.newpage(); grid::grid.text(paste("DIABLO plot error:", conditionMessage(r))); return(invisible()) }
  if (inherits(r, "ggplot") || (is.list(r) && inherits(r$graph, "ggplot"))) print(if (inherits(r, "ggplot")) r else r$graph)
  invisible()
}

byod_cheatmap <- function(d, genes, by = "Timepoint", split_k = 0L, seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  genes <- intersect(genes, rownames(d$expr))
  if (length(genes) < 2) { grid::grid.newpage(); grid::grid.text("Enter >= 2 valid genes"); return(invisible()) }
  M <- d$expr[genes, , drop = FALSE]; Z <- t(scale(t(M))); Z[!is.finite(Z)] <- 0
  samp <- colnames(Z); tp <- factor(d$meta[samp, "Timepoint"], levels = d$tp_levels)
  ord <- order(tp); Z <- Z[, ord, drop = FALSE]
  ann_df <- data.frame(Timepoint = droplevels(tp[ord]))
  tl <- levels(ann_df$Timepoint)
  ann_cols <- list(Timepoint = setNames(as.character(resolve_appearance(tl, PAL_NAMES[1], splitby = "Timepoint")), tl))
  if (!identical(by, "Timepoint") && by %in% names(d$meta)) {
    bv <- droplevels(grp_factor(d, by)[match(samp, rownames(d$meta))])[ord]
    if (nlevels(bv) >= 2) { ann_df[[by]] <- bv; bl <- levels(bv)
      ann_cols[[by]] <- setNames(as.character(resolve_appearance(bl, PAL_NAMES[2], splitby = by)), bl) } }
  ha <- ComplexHeatmap::HeatmapAnnotation(df = ann_df, col = ann_cols, annotation_name_side = "left",
          annotation_name_gp = grid::gpar(fontsize = 9))
  hm <- ComplexHeatmap::Heatmap(Z, name = "z-score",
    col = circlize::colorRamp2(c(-2, 0, 2), c(HUMAN, "#FFFFFF", PIG)),
    top_annotation = ha, cluster_columns = FALSE, column_order = seq_len(ncol(Z)),
    show_column_names = FALSE, show_row_names = (length(genes) <= 70), row_names_gp = grid::gpar(fontsize = 8),
    row_km = if (split_k >= 2) as.integer(split_k) else 1L, row_km_repeats = 20L,
    row_title = if (split_k >= 2) "cluster %s" else character(0), row_gap = grid::unit(1.5, "mm"),
    column_title = paste0(d$species, " — ", length(genes), " genes (row z-score)"),
    heatmap_legend_param = list(direction = "horizontal"))
  ComplexHeatmap::draw(hm, heatmap_legend_side = "bottom", annotation_legend_side = "bottom", merge_legend = TRUE)
  invisible()
}

cs_timepoints <- function(d)
  levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))

cs_align_timepoints <- function(dA, dB, align = c("label", "position")) {
  align <- match.arg(align)
  a <- cs_timepoints(dA); b <- cs_timepoints(dB)
  if (align == "label") {
    shared <- a[a %in% b]
    return(list(align = "label", A = shared, B = shared, label = shared,
                n = length(shared),
                note = if (!length(shared))
                  sprintf("The two studies share no timepoint label (%s vs %s). Align by position instead, if the k-th timepoint of one is comparable to the k-th of the other.",
                          paste(a, collapse = "/"), paste(b, collapse = "/")) else NULL))
  }
  k <- min(length(a), length(b))
  if (k < 2) return(list(align = "position", A = character(0), B = character(0),
                         label = character(0), n = 0,
                         note = "At least two ordered timepoints are needed in both studies."))
  list(align = "position", A = a[seq_len(k)], B = b[seq_len(k)],
       label = sprintf("%s / %s", a[seq_len(k)], b[seq_len(k)]), n = k,
       note = sprintf("Matched by position: timepoint %d of each study is treated as comparable. %s%s",
                      k,
                      if (length(a) != length(b))
                        sprintf("The studies have %d and %d timepoints, so only the first %d are used. ",
                                length(a), length(b), k) else "",
                      "Check that the elapsed times are actually comparable before interpreting."))
}

cs_shared_timepoints <- function(dA, dB) cs_align_timepoints(dA, dB, "label")$A

xs_organism <- function(d) {
  v <- d[["organism"]]
  if (is.null(v) || !length(v) || !nzchar(as.character(v)[1])) NA_character_
  else as.character(v)[1]
}

xs_ortholog_applicable <- function(dA, dB) {
  a <- xs_organism(dA); b <- xs_organism(dB)
  if (is.na(a) || is.na(b)) return(list(ok = FALSE, why = "declare both source organisms for cross-species comparison"))
  if (identical(a, b)) return(list(ok = FALSE, why = paste("both source organisms are", a)))
  ok <- annotation_get(dA)$config$mode == "human" && annotation_get(dB)$config$mode == "human"
  list(ok = ok, why = if (ok) sprintf("%s against %s, through human-symbol annotation", a, b)
       else "the studies have no common annotated gene namespace")
}

cs_feature_pairs <- function(dA, dB, one_to_one_only = TRUE) {
  a <- xs_organism(dA); b <- xs_organism(dB)
  cross <- !is.na(a) && !is.na(b) && !identical(a, b)
  none <- function(rule) list(features = character(), rule = rule, keyA = NULL, keyB = NULL)
  if (xor(is.na(a), is.na(b))) return(none("Declare both source organisms before matching features."))
  aa <- annotation_get(dA)$config; ab <- annotation_get(dB)$config
  same_native <- !cross && !is.na(a) && identical(aa$id_type, ab$id_type) && identical(aa$id_species, ab$id_species)
  if (same_native) return(list(features = intersect(rownames(dA$expr), rownames(dB$expr)),
    rule = paste("exact feature IDs in the same declared identifier namespace within", a), keyA = NULL, keyB = NULL))
  if (cross && !xs_ortholog_applicable(dA, dB)$ok)
    return(none("No compatible annotation for this cross-species comparison; native names are not ortholog evidence."))
  if (cross || (!is.na(a) && annotation_get(dA)$config$mode == "human" &&
                              annotation_get(dB)$config$mode == "human")) {
    keyA <- annotation_keys(dA); keyB <- annotation_keys(dB)
    shared <- intersect(names(keyA), names(keyB))
    rule <- "unique mapped human-symbol identifiers (see each layer's mapping report)"
    if (cross && isTRUE(dA$atlas_annotation) && isTRUE(dB$atlas_annotation) && one_to_one_only) {
      ok <- vapply(shared, function(g) {
        r <- ORTHOLOG_MAP[[g]]
        !is.null(r) && identical(as.character(r$class), "one_to_one")
      }, logical(1))
      shared <- shared[ok]
      rule <- "curated atlas human-pig 1:1 orthologs; no fallback when coverage is low"
    }
    return(list(features = shared, rule = rule, keyA = keyA[shared], keyB = keyB[shared]))
  }
  shared <- intersect(rownames(dA$expr), rownames(dB$expr))
  list(features = shared, rule = if (is.na(a)) "exact feature IDs; source organisms undeclared, not an orthology claim"
       else paste("exact native feature IDs within", a), keyA = NULL, keyB = NULL)
}

cs_scale_info <- function(d) {
  u <- tolower(trimws(d$unit %||% ""))
  fam <- NA_character_; off <- 0; src <- "declared unit"
  has_off <- function(z) grepl("[+]\\s*1|1p", z)
  if (grepl("log1p", u))                                { fam <- "loge";  off <- 1 }
  else if (grepl("log2|log\\s*2", u))                   { fam <- "log2";  off <- if (has_off(u)) 1 else 0 }
  else if (grepl("log10|log\\s*10", u))                 { fam <- "log10"; off <- if (has_off(u)) 1 else 0 }
  else if (grepl("\\bln\\b|natural log|loge", u))        { fam <- "loge";  off <- if (has_off(u)) 1 else 0 }

  else if (grepl("\\blog\\b", u)) {
    fam <- "loge"; off <- if (has_off(u)) 1 else 0
    src <- "declared unit (no base given, natural log assumed - state log2 or log10 if that is what you have)"
  }
  else if (grepl("z-?score|centred|centered|scaled", u)) fam <- "centered"

  else if (grepl("beta[- ]?value|\\bbeta\\b|proportion|fraction|methylation", u)) fam <- "bounded"
  else if (grepl("count|tpm|fpkm|rpkm|cpm|intensity|abundance|linear", u)) fam <- "linear"

  v <- suppressWarnings(as.numeric(d$expr)); v <- v[is.finite(v)]
  mx <- if (length(v)) max(v) else NA_real_; mn <- if (length(v)) min(v) else NA_real_
  neg <- is.finite(mn) && mn < -1e-6
  if (is.na(fam) || fam == "linear") {

    guess <- if (neg) "centered"
             else if (is.finite(mx) && is.finite(mn) && mx <= 1 + 1e-8 && mn >= -1e-8) "bounded"
             else if (is.finite(mx) && mx > 100) "linear"
             else if (is.finite(mx) && mx <= 30) "log2"
             else "unknown"
    if (is.na(fam)) { fam <- guess; src <- "inferred from the value range" }
    else if (fam == "linear" && guess %in% c("log2", "centered", "bounded")) {
      fam <- guess; src <- "inferred from the value range, overriding the declared unit"
    }
  }
  logged <- fam %in% c("log2", "log10", "loge")
  list(family = fam, offset = off, logged = logged, negative = neg, source = src,
       min = mn, max = mx, declared = if (nzchar(u)) d$unit else "not declared",
       label = switch(fam,
         log2 = sprintf("log2 scale%s", if (off) " with a +1 offset" else ""),
         log10 = sprintf("log10 scale%s", if (off) " with a +1 offset" else ""),
         loge = sprintf("natural-log scale%s", if (off) " with a +1 offset" else ""),
         linear = "linear scale (not log-transformed)",
         centered = "centred or z-scored (has negative values)",
         bounded = "bounded proportion between 0 and 1 (e.g. a methylation beta value)",
         "unrecognised scale"))
}

cs_to_linear <- function(M, si) {
  switch(si$family,
    log2  = 2^M - si$offset,
    log10 = 10^M - si$offset,
    loge  = exp(M) - si$offset,
    M)                                   # linear / centered / unknown: already there
}

CS_MEASURES <- c("difference", "ratio", "percent")

cs_measure_advice <- function(si) {

  if (identical(si$family, "bounded"))
    return(list(
      difference = list(status = "recommended",
        why = "Values are proportions between 0 and 1, so a difference is the delta beta: directly interpretable and comparable across features."),
      ratio = list(status = "discouraged",
        why = "A ratio of two proportions is defined but behaves badly near zero, where a change of 0.01 becomes a large fold change."),
      percent = list(status = "usable",
        why = "Percent change relative to the baseline proportion. Interpretable, but like the ratio it is unstable for features with a near-zero baseline.")))
  if (si$family == "centered")
    return(list(
      difference = list(status = "recommended",
        why = "Values are centred or z-scored, so a difference is the only interpretable change."),
      ratio = list(status = "not meaningful",
        why = "A ratio needs a meaningful zero. These values are centred and cross zero, so the ratio is unstable and can flip sign."),
      percent = list(status = "not meaningful",
        why = "Percent change needs a meaningful zero. These values are centred and cross zero.")))
  if (si$logged)
    return(list(
      difference = list(status = "recommended",
        why = sprintf("The data are on a %s. A difference of logs IS the log ratio, which is what a fold change means.", si$label)),
      ratio = list(status = "usable",
        why = "Computed by returning to the linear scale first, then taking log2 of the ratio. On log data this is close to the difference and adds nothing, but it is not wrong."),
      percent = list(status = "usable",
        why = "Computed on the back-transformed linear values. Interpretable, but sensitive to features with a near-zero baseline.")))
  if (si$family == "linear")
    return(list(
      difference = list(status = "discouraged",
        why = "The data are linear, so a difference is scale-dependent: a change of 100 means something quite different for an abundant and a rare feature."),
      ratio = list(status = "recommended",
        why = "On linear data the log2 ratio is the standard fold change and is comparable across features."),
      percent = list(status = "usable",
        why = "Interpretable on linear data, but unstable for features with a near-zero baseline.")))
  list(
    difference = list(status = "usable", why = "The scale could not be resolved from the declared unit or the value range. Confirm the unit before interpreting."),
    ratio      = list(status = "usable", why = "The scale could not be resolved. Confirm the unit before interpreting."),
    percent    = list(status = "usable", why = "The scale could not be resolved. Confirm the unit before interpreting."))
}

cs_method_advice <- function(o = NULL) {
  list(
    pearson = list(status = "recommended",
      why = "Measures how close the two studies' responses are to a straight line, so it rewards agreement in magnitude as well as direction; its square is the R-squared reported by the reference sensitivity analysis."),
    spearman = list(status = "usable",
      why = "Uses ranks only, so a handful of extreme features cannot dominate it. Report it alongside Pearson when the scatter has long tails; if the two disagree, outliers are driving the Pearson value."))
}

cs_subject_delta <- function(d, baseline, endpoint, measure = c("difference", "ratio", "percent"),
                             subject_col = NULL) {
  measure <- match.arg(measure)
  subject_col <- subject_col %||% d$subject %||% "Kidney_ID"
  tp  <- as.character(d$meta$Timepoint)
  si  <- cs_scale_info(d)
  lin <- function(M) cs_to_linear(M, si)

  rel <- function(b, e) switch(measure,
      difference = e - b,

      ratio      = { bl <- lin(b); el <- lin(e)
                     ok <- is.finite(bl) & is.finite(el) & bl > 0 & el > 0
                     out <- rep(NA_real_, length(b)); out[ok] <- log2(el[ok] / bl[ok]); out },
      percent    = { bl <- lin(b); el <- lin(e)
                     ok <- is.finite(bl) & is.finite(el) & abs(bl) > 1e-8
                     out <- rep(NA_real_, length(b)); out[ok] <- 100 * (el[ok] - bl[ok]) / abs(bl[ok]); out })

  if (subject_col %in% names(d$meta)) {
    subj <- as.character(d$meta[[subject_col]])
    keep <- unique(subj[!is.na(subj)])
    cols <- lapply(keep, function(s) {
      ib <- which(subj == s & tp == baseline); ie <- which(subj == s & tp == endpoint)
      if (!length(ib) || !length(ie)) return(NULL)
      rel(rowMeans(d$expr[, ib, drop = FALSE], na.rm = TRUE),
          rowMeans(d$expr[, ie, drop = FALSE], na.rm = TRUE))
    })
    ok <- !vapply(cols, is.null, logical(1))
    if (any(ok)) {
      M <- do.call(cbind, cols[ok]); rownames(M) <- rownames(d$expr); colnames(M) <- keep[ok]
      return(list(delta = M, subjects = keep[ok], measure = measure, pairing = "paired",
                  baseline = baseline, endpoint = endpoint, subject_col = subject_col))
    }
  }

  ib <- which(tp == baseline); ie <- which(tp == endpoint)
  if (!length(ib) || !length(ie)) return(NULL)
  b <- rowMeans(d$expr[, ib, drop = FALSE], na.rm = TRUE)
  M <- vapply(ie, function(j) rel(b, d$expr[, j]), numeric(nrow(d$expr)))
  if (!is.matrix(M)) M <- matrix(M, nrow = nrow(d$expr))
  rownames(M) <- rownames(d$expr); colnames(M) <- colnames(d$expr)[ie]
  list(delta = M, subjects = colnames(M), measure = measure, pairing = "unpaired",
       n_baseline = length(ib), baseline = baseline, endpoint = endpoint,
       subject_col = if (subject_col %in% names(d$meta)) subject_col else NA_character_)
}

CS_MEASURE_LAB <- c(
  difference = "Δ log2 (endpoint − baseline)",
  ratio      = "log2 ratio (endpoint / baseline)",
  percent    = "% change from baseline")

cs_trajectory_one <- function(d, baseline, endpoints, measure, subject_col,
                              aggregate_to = NULL, meta_for_agg = NULL) {
  mats <- lapply(endpoints, function(e) {
    dd <- cs_subject_delta(d, baseline, e, measure, subject_col)
    if (is.null(dd)) return(NULL)
    M <- dd$delta

    if (!is.null(aggregate_to) && !is.null(meta_for_agg)) {
      key <- setNames(as.character(meta_for_agg[[aggregate_to]])[
                        match(colnames(M), as.character(meta_for_agg[[subject_col]]))],
                      colnames(M))
      if (!anyNA(key)) {
        u <- unique(key)
        M <- vapply(u, function(a) rowMeans(M[, key == a, drop = FALSE], na.rm = TRUE),
                    numeric(nrow(M)))
        colnames(M) <- u
      }
    }
    list(mean = rowMeans(M, na.rm = TRUE), n_units = ncol(M), units = colnames(M),
         pairing = dd$pairing %||% "paired", n_baseline = dd$n_baseline %||% NA_integer_)
  })
  ok <- !vapply(mats, is.null, logical(1))
  if (!any(ok)) return(NULL)
  traj <- do.call(cbind, lapply(mats[ok], `[[`, "mean"))
  colnames(traj) <- endpoints[ok]
  rownames(traj) <- rownames(d$expr)
  pr <- unique(vapply(mats[ok], `[[`, character(1), "pairing"))
  list(traj = traj, endpoints = endpoints[ok],
       n_units = max(vapply(mats[ok], `[[`, integer(1), "n_units")),

       pairing = if (length(pr) == 1) pr else "mixed",
       n_baseline = mats[ok][[1]]$n_baseline,
       units = mats[ok][[1]]$units)
}

cs_trajectory <- function(dA, dB, baseline, measure = "difference",
                          one_to_one_only = TRUE, method = c("pearson", "spearman"),
                          labA = "A", labB = "B", align = c("label", "position", "manual"),
                          manual = NULL) {
  method <- match.arg(method); align <- match.arg(align)

  al <- xs_align_timepoints(dA, dB, align, manual = manual)
  if (al$n < 2) return(list(ok = FALSE, msg = al$note %||% "The two studies have too few comparable timepoints."))

  i <- match(baseline, al$label)
  if (is.na(i)) i <- match(baseline, al$A)
  if (is.na(i) || i >= al$n)
    return(list(ok = FALSE, msg = "Choose a baseline with at least one later comparable timepoint."))
  idx  <- (i + 1):al$n
  epsA <- al$A[idx]; epsB <- al$B[idx]; eps_lab <- al$label[idx]
  baseA <- al$A[i];  baseB <- al$B[i]

  fp <- cs_feature_pairs(dA, dB, one_to_one_only)
  if (length(fp$features) < 10)
    return(list(ok = FALSE, msg = sprintf("Only %d features are shared by the two studies; they must use the same feature identifiers.",
                                          length(fp$features))))

  roll <- function(d) if ("Pig_Number" %in% names(d$meta)) "Pig_Number" else NULL
  tA <- cs_trajectory_one(dA, baseA, epsA, measure, dA$subject %||% "Kidney_ID",
                          aggregate_to = roll(dA), meta_for_agg = dA$meta)
  tB <- cs_trajectory_one(dB, baseB, epsB, measure, dB$subject %||% "Kidney_ID",
                          aggregate_to = roll(dB), meta_for_agg = dB$meta)
  if (is.null(tA) || is.null(tB))
    return(list(ok = FALSE, msg = sprintf("%s has no profile at the chosen baseline, so nothing can be expressed relative to it.",
                                          if (is.null(tA)) labA else labB)))

  keep_pos <- which(epsA %in% tA$endpoints & epsB %in% tB$endpoints)
  if (!length(keep_pos)) return(list(ok = FALSE, msg = "No post-baseline timepoint survived in both studies."))
  epsA <- epsA[keep_pos]; epsB <- epsB[keep_pos]; eps <- eps_lab[keep_pos]
  P <- tA$traj[, epsA, drop = FALSE]; H <- tB$traj[, epsB, drop = FALSE]
  colnames(P) <- colnames(H) <- eps

  g  <- fp$features
  rA <- if (is.null(fp$keyA)) setNames(g, g) else fp$keyA[g]
  rB <- if (is.null(fp$keyB)) setNames(g, g) else fp$keyB[g]
  ok <- !is.na(rA) & !is.na(rB) & rA %in% rownames(P) & rB %in% rownames(H)
  g <- g[ok]; rA <- rA[ok]; rB <- rB[ok]
  P <- P[rA, , drop = FALSE]; H <- H[rB, , drop = FALSE]
  rownames(P) <- rownames(H) <- g

  keep <- apply(is.finite(P), 1, all) & apply(is.finite(H), 1, all)
  n_dropped <- sum(!keep)
  P <- P[keep, , drop = FALSE]; H <- H[keep, , drop = FALSE]; g <- g[keep]
  if (length(g) < 10) return(list(ok = FALSE, msg = sprintf(
    "Only %d features have a finite value at every timepoint in both studies (%d were dropped). With '%s' a feature is dropped wherever it is zero or missing at the baseline or at any later timepoint.",
    length(g), n_dropped, measure)))

  long <- data.frame(
    feature   = rep(g, times = length(eps)),
    timepoint = factor(rep(eps, each = length(g)), levels = eps),
    A         = as.vector(P),
    B         = as.vector(H), stringsAsFactors = FALSE)
  long$agreement <- ifelse(sign(long$A) == sign(long$B), "same direction", "opposite")

  by_tp <- do.call(rbind, lapply(eps, function(e) {
    z <- long[long$timepoint == e, ]
    ct <- suppressWarnings(stats::cor.test(z$A, z$B, method = method))
    data.frame(timepoint = e, r = round(unname(ct$estimate), 4),
               r2 = round(unname(ct$estimate)^2, 4), p = ct$p.value,
               n_features = nrow(z),
               same_direction_pct = round(100 * mean(z$agreement == "same direction"), 1),
               stringsAsFactors = FALSE)
  }))
  by_tp$timepoint <- factor(by_tp$timepoint, levels = eps)

  n_t <- length(eps)
  gene_r <- if (n_t >= 3) {
    Pc <- P - rowMeans(P); Hc <- H - rowMeans(H)
    den <- sqrt(rowSums(Pc^2) * rowSums(Hc^2))
    ifelse(den > 0, rowSums(Pc * Hc) / den, NA_real_)
  } else rep(NA_real_, length(g))
  auc <- function(M) {
    if (ncol(M) == 1L) return(as.vector(M))
    w <- rep(1, ncol(M)); w[1] <- w[ncol(M)] <- 0.5
    as.vector(M %*% w) / (ncol(M) - 1)
  }
  peak <- function(M) eps[max.col(abs(M), ties.method = "first")]
  per_gene <- data.frame(
    feature = g, trajectory_r = round(gene_r, 4),
    A_auc = round(auc(P), 4), B_auc = round(auc(H), 4),
    A_peak = peak(P), B_peak = peak(H),
    A_max = round(P[cbind(seq_len(nrow(P)), max.col(abs(P), ties.method = "first"))], 4),
    B_max = round(H[cbind(seq_len(nrow(H)), max.col(abs(H), ties.method = "first"))], 4),
    stringsAsFactors = FALSE)
  per_gene$same_peak <- per_gene$A_peak == per_gene$B_peak
  per_gene$auc_agreement <- ifelse(sign(per_gene$A_auc) == sign(per_gene$B_auc),
                                   "same direction", "opposite")
  per_gene$mean_abs_auc <- round((abs(per_gene$A_auc) + abs(per_gene$B_auc)) / 2, 4)
  per_gene <- per_gene[order(-per_gene$mean_abs_auc), ]
  names(per_gene) <- sub("^A_", paste0(labA, "_"), names(per_gene))
  names(per_gene) <- sub("^B_", paste0(labB, "_"), names(per_gene))

  auc_ct <- suppressWarnings(stats::cor.test(auc(P), auc(H), method = method))
  pooled <- suppressWarnings(stats::cor.test(long$A, long$B, method = method))

  list(ok = TRUE, long = long, by_tp = by_tp, per_gene = per_gene,
       P = P, H = H, endpoints = eps, features = g,
       labA = labA, labB = labB,
       scaleA = cs_scale_info(dA), scaleB = cs_scale_info(dB),
       method = method, measure = measure, baseline = baseline,
       align = al$align, align_note = al$note, baselineA = baseA, baselineB = baseB,
       rule = fp$rule, n_features = length(g), n_dropped = n_dropped,
       unitsA = tA$n_units, unitsB = tB$n_units,
       pairingA = tA$pairing, pairingB = tB$pairing,

       pairing_note = {
         un <- c(labA, labB)[c(!identical(tA$pairing, "paired"), !identical(tB$pairing, "paired"))]
         if (!length(un)) NULL
         else sprintf("%s referenced to the baseline group mean (%s no subject measured at both the baseline and a later timepoint), so %s trajectory is unpaired and carries between-subject variation.",
                      paste(un, collapse = " and "),
                      if (length(un) == 1) "it has" else "they have",
                      if (length(un) == 1) "that" else "each")
       },
       r_pooled = unname(pooled$estimate), r_auc = unname(auc_ct$estimate),
       median_trajectory_r = if (all(is.na(gene_r))) NA_real_ else round(stats::median(gene_r, na.rm = TRUE), 4),
       pct_same_peak = round(100 * mean(per_gene$same_peak), 1),
       ttl = sprintf("%s vs %s · baseline %s · %d timepoints (%s) · %s features · %d vs %d independent units%s",
                     labA, labB,
                     if (identical(baseA, baseB)) baseA else sprintf("%s / %s", baseA, baseB),
                     length(eps), paste(eps, collapse = ", "),
                     format(length(g), big.mark = ","), tA$n_units, tB$n_units,
                     if (identical(tA$pairing, "paired") && identical(tB$pairing, "paired")) ""
                     else sprintf(" · %s / %s baseline", tA$pairing, tB$pairing)))
}

cs_traj_scatter <- function(o, cols = NULL, max_points = 6000L, seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  if (!isTRUE(o$ok)) return(emptyplot(o$msg))
  d <- o$long
  if (nrow(d) > max_points) {                      # thin for the browser, keep every timepoint
    d <- do.call(rbind, lapply(split(d, d$timepoint), function(s)
      s[sample.int(nrow(s), min(nrow(s), floor(max_points / nlevels(d$timepoint)))), , drop = FALSE]))
  }
  cc <- if (is.null(cols)) c("same direction" = PIG, "opposite" = GOLD) else cols

  lim <- stats::quantile(abs(c(o$long$A, o$long$B)), .995, na.rm = TRUE) * 1.05
  n_out <- sum(abs(o$long$A) > lim | abs(o$long$B) > lim, na.rm = TRUE)
  ax <- CS_MEASURE_LAB[[o$measure]]
  eps <- levels(o$long$timepoint)
  plots <- lapply(seq_along(eps), function(k) {
    e <- eps[k]; s <- d[d$timepoint == e, ]
    st <- o$by_tp[o$by_tp$timepoint == e, ]
    plot_ly(s, x = ~B, y = ~A, type = "scatter", mode = "markers",
            color = ~factor(agreement, names(cc)), colors = unname(cc),
            marker = list(size = 3.4, opacity = .38), showlegend = (k == 1),
            text = ~feature, hoverinfo = "text") %>%
      add_segments(x = -lim, xend = lim, y = -lim, yend = lim, inherit = FALSE,
                   line = list(dash = "dash", color = "#999", width = 1), showlegend = FALSE) %>%
      add_annotations(x = -lim * .92, y = lim * .88, xref = "x", yref = "y",
                      text = sprintf("<b>%s</b><br>r = %.3f<br>R<sup>2</sup> = %.3f",
                                     e, st$r, st$r2),
                      showarrow = FALSE, align = "left",
                      font = list(size = 10.5, color = "#222"))
  })
  subplot(plots, nrows = 1, shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE) %>%
    ply(margin = list(l = 72, r = 20, t = 66, b = 76),
        legend = list(orientation = "h", x = .5, xanchor = "center", y = -.2, yanchor = "top"),
        annotations = list(
          list(x = .5, y = -.15, xref = "paper", yref = "paper", showarrow = FALSE,
               text = paste(o$labB, ax), font = list(size = 12)),
          list(x = -.055, y = .5, xref = "paper", yref = "paper", showarrow = FALSE,
               textangle = -90, text = paste(o$labA, ax), font = list(size = 12))),
        title = list(text = paste0(o$ttl, "<br><sub>every feature at every timepoint · matched by ",
                                   o$rule,
                                   if (n_out > 0) sprintf(" · axes clipped to the central 99.5%%, %s points outside",
                                                          format(n_out, big.mark = ",")) else "",
                                   "</sub>"),
                     x = .02, font = list(size = 12.5)))
}

cs_traj_hist <- function(o) {
  if (!isTRUE(o$ok)) return(emptyplot(o$msg))
  r <- o$per_gene$trajectory_r; r <- r[is.finite(r)]
  if (!length(r)) return(emptyplot("Need at least three post-baseline timepoints for a per-feature trajectory correlation"))
  med <- stats::median(r)
  plot_ly(x = r, type = "histogram", nbinsx = 60, marker = list(color = PIG, opacity = .8),
          showlegend = FALSE,
          hovertemplate = "r %{x:.2f}<br>%{y} features<extra></extra>") %>%
    add_segments(x = med, xend = med, y = 0, yend = 1, yref = "paper", inherit = FALSE,
                 line = list(color = GOLD, width = 2, dash = "dash"),
                 name = sprintf("median %.3f", med)) %>%
    ply(xaxis = list(title = sprintf("per-feature trajectory correlation (%s vs %s)", o$labA, o$labB), range = c(-1, 1),
                     zeroline = TRUE),
        yaxis = list(title = "features"),
        legend = list(orientation = "h", y = 1.06),
        title = list(text = sprintf("Do individual features follow the same time course? · median r = %.3f · %s features over %d timepoints",
                                    med, format(length(r), big.mark = ","), length(o$endpoints)),
                     x = .02, font = list(size = 12.5)))
}

cs_traj_gene <- function(o, gene, cols = NULL) {
  if (!isTRUE(o$ok)) return(emptyplot(o$msg))
  if (!gene %in% rownames(o$P)) return(emptyplot(paste(gene, "has no complete trajectory in both studies")))
  eps <- o$endpoints
  df <- rbind(
    data.frame(timepoint = factor(eps, levels = eps), value = as.numeric(o$P[gene, ]), study = o$labA),
    data.frame(timepoint = factor(eps, levels = eps), value = as.numeric(o$H[gene, ]), study = o$labB))
  df$study <- factor(df$study, levels = c(o$labA, o$labB))
  pg <- o$per_gene[o$per_gene$feature == gene, ]
  cc <- if (is.null(cols)) setNames(c(PIG, HUMAN), c(o$labA, o$labB)) else cols
  plot_ly(df, x = ~timepoint, y = ~value, color = ~study, colors = unname(cc[c(o$labA, o$labB)]),
          type = "scatter", mode = "lines+markers", line = list(width = 2.5),
          marker = list(size = 8)) %>%
    add_segments(x = eps[1], xend = eps[length(eps)], y = 0, yend = 0, inherit = FALSE,
                 line = list(color = "#aaa", width = 1), showlegend = FALSE) %>%
    ply(xaxis = list(title = paste0("timepoint (relative to ", o$baseline, ")")),
        yaxis = list(title = CS_MEASURE_LAB[[o$measure]], zeroline = TRUE),
        legend = list(orientation = "h", y = 1.06),
        title = list(text = sprintf("%s · trajectory r = %s · peak: %s %s, %s %s",
                                    gene,
                                    if (nrow(pg) && is.finite(pg$trajectory_r)) sprintf("%.3f", pg$trajectory_r) else "n/a",
                                    o$labA, if (nrow(pg)) pg[[paste0(o$labA, "_peak")]] else "?",
                                    o$labB, if (nrow(pg)) pg[[paste0(o$labB, "_peak")]] else "?"),
                     x = .02, font = list(size = 12.5)))
}

cs_noise_ceiling <- function(dA, dB, baseline, measure = "difference",
                             one_to_one_only = TRUE, labA = "A", labB = "B",
                             align = c("label", "position", "manual"), n_perm = 25L, seed = trace_seed(),
                             manual = NULL) {
  align <- match.arg(align)
  al <- xs_align_timepoints(dA, dB, align, manual = manual)
  i <- match(baseline, al$label); if (is.na(i)) i <- match(baseline, al$A)
  if (is.na(i) || i >= al$n) return(NULL)
  idx <- (i + 1):al$n
  fp <- cs_feature_pairs(dA, dB, one_to_one_only)
  one <- function(d, lab, base_tp, eps) {
    sc <- d$subject %||% "Kidney_ID"
    if (!sc %in% names(d$meta)) return(NULL)
    subs <- unique(as.character(d$meta[[sc]]))
    if (length(subs) < 6) return(NULL)
    restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
    rs <- vapply(seq_len(n_perm), function(b) {
      half <- sample(subs, floor(length(subs) / 2))
      sub_traj <- function(keep) {
        sel <- as.character(d$meta[[sc]]) %in% keep
        dd <- d; dd$expr <- d$expr[, sel, drop = FALSE]; dd$meta <- d$meta[sel, , drop = FALSE]
        cs_trajectory_one(dd, base_tp, eps, measure, sc)
      }
      a <- sub_traj(half); z <- sub_traj(setdiff(subs, half))
      if (is.null(a) || is.null(z)) return(NA_real_)
      ep <- intersect(a$endpoints, z$endpoints)
      gg <- intersect(rownames(a$traj), fp$features)
      if (!length(ep) || length(gg) < 10) return(NA_real_)
      A <- a$traj[gg, ep, drop = FALSE]; Z <- z$traj[gg, ep, drop = FALSE]
      ok <- stats::complete.cases(A) & stats::complete.cases(Z)
      if (sum(ok) < 10) return(NA_real_)
      suppressWarnings(stats::cor(as.vector(A[ok, ]), as.vector(Z[ok, ])))
    }, numeric(1))
    rs <- rs[is.finite(rs)]
    if (!length(rs)) return(NULL)
    data.frame(study = lab, n_perm = length(rs),
               median_r = round(stats::median(rs), 4),
               q25 = round(unname(stats::quantile(rs, .25)), 4),
               q75 = round(unname(stats::quantile(rs, .75)), 4),
               stringsAsFactors = FALSE)
  }
  out <- Filter(Negate(is.null),
                list(one(dA, labA, al$A[i], al$A[idx]), one(dB, labB, al$B[i], al$B[idx])))
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

MIN_UNITS$trajectory_consistency <- 3L   # subjects whose trajectories can be compared
MIN_UNITS$baseline_sensitivity   <- 3L   # subjects contributing to each baseline

byod_traj_consistency <- function(d, min_tp = 3L, max_feat = 4000L,
                                  min_units = MIN_UNITS$trajectory_consistency) {
  subj_col <- d$subject
  if (is.null(subj_col) || !subj_col %in% names(d$meta))
    return(list(ok = FALSE, msg = "Trajectory consistency needs a declared subject: it compares each subject's trajectory with the other subjects'. Declare the repeatedly sampled unit and re-run."))
  tf <- droplevels(factor(d$meta$Timepoint, levels = d$tp_levels))
  tps <- levels(tf)
  if (length(tps) < min_tp)
    return(list(ok = FALSE, msg = sprintf("Needs >= %d timepoints to correlate a trajectory (has %d).", min_tp, length(tps))))
  subj <- as.character(d$meta[[subj_col]])

  per_subj <- vapply(unique(subj), function(s)
    length(unique(as.character(tf)[subj == s & !is.na(tf)])), integer(1))
  usable <- names(per_subj)[per_subj >= min_tp]
  if (length(usable) < min_units)
    return(list(ok = FALSE, msg = sprintf(
      "Trajectory consistency needs at least %d subjects sampled at >= %d timepoints; this design has %d (of %d %ss, from %d profiles).",
      min_units, min_tp, length(usable), length(per_subj), subj_col, ncol(d$expr))))

  nest_col <- d$nest
  nest_map <- NULL; n_indep <- length(usable)
  if (!is.null(nest_col) && nest_col %in% names(d$meta)) {
    nm <- vapply(usable, function(s) {
      z <- unique(as.character(d$meta[[nest_col]][subj == s])); z[!is.na(z)][1] %||% NA_character_
    }, character(1))
    if (!anyNA(nm)) {
      nest_map <- setNames(nm, usable)
      n_indep  <- length(unique(nm))
      if (n_indep < min_units)
        return(list(ok = FALSE, msg = sprintf(
          "Trajectory consistency needs at least %d independent units; the %d usable %ss nest inside %d %ss.",
          min_units, length(usable), subj_col, n_indep, nest_col)))
    }
  }

  rv <- byod_rowvars(d$expr)
  keep <- head(order(rv, decreasing = TRUE), min(max_feat, sum(rv > 0)))
  E <- d$expr[keep, , drop = FALSE]

  A <- array(NA_real_, dim = c(nrow(E), length(usable), length(tps)),
             dimnames = list(rownames(E), usable, tps))
  for (s in usable) for (tp in tps) {
    j <- which(subj == s & as.character(tf) == tp)
    if (length(j)) A[, s, tp] <- if (length(j) == 1L) E[, j] else rowMeans(E[, j, drop = FALSE], na.rm = TRUE)
  }

  n_f <- dim(A)[1]
  mean_r <- sd_r <- min_r <- rep(NA_real_, n_f)
  n_used <- integer(n_f)
  for (i in seq_len(n_f)) {
    M <- t(A[i, , , drop = TRUE])                      # timepoints x subjects
    if (is.null(dim(M))) next
    ok_s <- colSums(!is.na(M)) >= min_tp
    if (sum(ok_s) < min_units) { n_used[i] <- sum(ok_s); next }
    M <- M[, ok_s, drop = FALSE]; n_used[i] <- ncol(M)

    nest_of_col <- nest_map[colnames(M)]
    rs <- vapply(seq_len(ncol(M)), function(j) {
      drop_j <- if (is.null(nest_map)) j else which(nest_of_col == nest_of_col[j])
      if (length(drop_j) >= ncol(M)) return(NA_real_)
      others <- rowMeans(M[, -drop_j, drop = FALSE], na.rm = TRUE)
      ok <- is.finite(M[, j]) & is.finite(others)
      if (sum(ok) < min_tp) return(NA_real_)
      a <- M[ok, j]; b <- others[ok]
      if (stats::sd(a) == 0 || stats::sd(b) == 0) return(NA_real_)
      suppressWarnings(stats::cor(a, b))
    }, numeric(1))
    rs <- rs[is.finite(rs)]
    if (length(rs)) { mean_r[i] <- mean(rs); sd_r[i] <- stats::sd(rs); min_r[i] <- min(rs) }
  }

  tab <- data.frame(Feature = rownames(E), mean_r = mean_r, sd_r = sd_r, min_r = min_r,
                    n_subjects = n_used, stringsAsFactors = FALSE)
  tab <- tab[is.finite(tab$mean_r), , drop = FALSE]
  if (!nrow(tab))
    return(list(ok = FALSE, msg = "No feature had enough overlapping timepoints across subjects to correlate."))
  tab <- tab[order(-tab$mean_r), ]
  unit_lab <- if (is.null(nest_map)) subj_col else nest_col
  list(ok = TRUE, tab = tab, subj_col = subj_col, nest_col = nest_col,
       n_subjects = length(usable), n_units = n_indep, unit = unit_lab,
       nested = !is.null(nest_map), timepoints = tps, n_features = nrow(tab),
       ttl = sprintf("Trajectory consistency across %d %s%s over %d timepoints (%d features)",
                     n_indep, unit_lab, if (is.null(nest_map)) "s"
                       else sprintf("s (%d %ss nested inside them)", length(usable), subj_col),
                     length(tps), nrow(tab)),
       note = sprintf("Leave-one-out: each %s's trajectory is correlated with the mean of the others. %s",
                      subj_col,
                      if (is.null(nest_map))
                        "No level above the subject was declared, so subjects are the units."
                      else sprintf("The whole %s is left out each time, so a sibling %s never stays in the comparison mean.",
                                   nest_col, subj_col)))
}

byod_traj_consistency_gg <- function(obj, fill = "#4C8FBF") {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  df  <- obj$tab
  med <- stats::median(df$mean_r, na.rm = TRUE)
  q1  <- stats::quantile(df$mean_r, .25, na.rm = TRUE)
  q3  <- stats::quantile(df$mean_r, .75, na.rm = TRUE)
  p <- ggplot(df, aes(x = "", y = mean_r)) +
    geom_hline(yintercept = 0,  colour = "grey86", linewidth = .8) +
    geom_hline(yintercept = q1, colour = "grey80", linetype = "dotted", linewidth = .7) +
    geom_hline(yintercept = q3, colour = "grey80", linetype = "dotted", linewidth = .7) +
    geom_hline(yintercept = med, colour = "#1F5C8B", linewidth = 1) +
    geom_violin(fill = fill, colour = "grey20", linewidth = .8,
                scale = "width", trim = FALSE, width = .68) +
    geom_boxplot(width = .10, outlier.shape = NA, fill = "white",
                 colour = "grey20", linewidth = .75) +

    annotate("text", x = 1.44, y = med, hjust = 1,
             vjust = if (med > 0.75) 1.6 else -0.8, size = 3.9,
             colour = "#1F5C8B", fontface = "bold",
             label = sprintf("median = %.3f", med)) +
    coord_cartesian(ylim = c(-1, 1), clip = "off") +
    scale_y_continuous(breaks = seq(-1, 1, .5), labels = sprintf("%.1f", seq(-1, 1, .5)),
                       expand = expansion(mult = c(0, 0))) +

    labs(title = paste(strwrap(obj$ttl, width = 62), collapse = "\n"),
         subtitle = paste(strwrap(sprintf("%s Computed on %s features.", obj$note,
                                          format(nrow(df), big.mark = ",")),
                                  width = 82), collapse = "\n"), x = NULL,
         y = "Correlation to the other units' mean trajectory") +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.15),
          axis.title.y  = element_text(size = 10.5, face = "bold"),
          axis.line     = element_line(linewidth = .8, colour = "black"),
          axis.ticks    = element_line(linewidth = .8, colour = "black"),
          axis.ticks.length = unit(.16, "cm"),
          axis.text.x   = element_blank(), axis.ticks.x = element_blank(),
          plot.margin   = margin(10, 18, 10, 10))
  p
}

byod_traj_consistency_plot <- function(obj, fill = "#4C8FBF") {
  g <- byod_traj_consistency_gg(obj, fill)
  if (inherits(g, "plotly")) return(g)
  ggplotly(g) %>% ply()
}

byod_baseline_sensitivity <- function(d, baseline_a, baseline_b, test_tps,
                                      fdr = 0.05, subject_col = NULL, paired = TRUE) {
  lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
  test_tps <- unique(as.character(test_tps))
  if (!length(test_tps)) return(list(ok = FALSE, msg = "Pick at least one timepoint to test."))
  for (x in c(baseline_a, baseline_b, test_tps))
    if (!x %in% lv) return(list(ok = FALSE, msg = sprintf("Timepoint '%s' is not in this study.", x)))
  if (baseline_a == baseline_b)
    return(list(ok = FALSE, msg = "Pick two DIFFERENT baselines; comparing a baseline with itself answers nothing."))
  bad <- intersect(test_tps, c(baseline_a, baseline_b))
  if (length(bad))
    return(list(ok = FALSE, msg = sprintf("A test timepoint cannot also be a baseline (%s).",
                                          paste(bad, collapse = ", "))))

  test_tps <- lv[lv %in% test_tps]
  sc <- subject_col %||% d$subject

  one <- function(b, tp) byod_run_limma(d, "Timepoint", b, tp,
                                        subject_col = sc, paired = paired && !is.null(sc))
  rows <- list(); stats <- list(); declined <- list()
  for (tp in test_tps) {
    A <- tryCatch(one(baseline_a, tp), error = function(e) e)
    B <- tryCatch(one(baseline_b, tp), error = function(e) e)
    why <- NULL
    for (nm in c("a", "b")) {
      r <- if (nm == "a") A else B; b <- if (nm == "a") baseline_a else baseline_b
      if (inherits(r, "error")) why <- sprintf("contrast against '%s' failed: %s", b, conditionMessage(r))
      else if (!is.null(attr(r, "declined"))) why <- sprintf("contrast against '%s' declined: %s", b, attr(r, "declined"))
      if (!is.null(why)) break
    }

    if (!is.null(why)) { declined[[tp]] <- why; next }

    fa <- A$Feature %||% rownames(A); fb <- B$Feature %||% rownames(B)
    common <- intersect(fa, fb)
    if (length(common) < 10) { declined[[tp]] <- "the two contrasts share too few features"; next }
    ia <- match(common, fa); ib <- match(common, fb)
    df <- data.frame(Timepoint = tp, Feature = common,
                     logFC_a = A$logFC[ia], logFC_b = B$logFC[ib],
                     q_a = A$adj.P.Val[ia], q_b = B$adj.P.Val[ib], stringsAsFactors = FALSE)
    fin <- is.finite(df$logFC_a) & is.finite(df$logFC_b)
    r   <- suppressWarnings(stats::cor(df$logFC_a[fin], df$logFC_b[fin]))
    rho <- suppressWarnings(stats::cor(df$logFC_a[fin], df$logFC_b[fin], method = "spearman"))
    sa <- df$Feature[is.finite(df$q_a) & df$q_a < fdr]
    sb <- df$Feature[is.finite(df$q_b) & df$q_b < fdr]
    inter <- length(intersect(sa, sb)); uni <- length(union(sa, sb))
    df$agreement <- ifelse(df$Feature %in% sa & df$Feature %in% sb, "significant from both",
                    ifelse(df$Feature %in% sa, sprintf("only from %s", baseline_a),
                    ifelse(df$Feature %in% sb, sprintf("only from %s", baseline_b),
                           "significant from neither")))
    rows[[tp]]  <- df
    stats[[tp]] <- data.frame(Timepoint = tp, r = r, r2 = r^2, rho = rho,
                              n_a = length(sa), n_b = length(sb), n_both = inter,
                              only_a = length(setdiff(sa, sb)), only_b = length(setdiff(sb, sa)),
                              n_union = uni,
                              jaccard = if (uni) inter / uni else NA_real_,
                              pct_overlap = if (uni) 100 * inter / uni else NA_real_,
                              n_units = attr(A, "n_units"), unit = attr(A, "unit"),
                              stringsAsFactors = FALSE)
  }
  if (!length(rows))
    return(list(ok = FALSE,
                msg = paste("No timepoint could be compared.",
                            paste(sprintf("%s: %s", names(declined), unlist(declined)), collapse = "; "))))

  tab <- do.call(rbind, rows); st <- do.call(rbind, stats)
  tab$Timepoint <- factor(tab$Timepoint, levels = test_tps)
  st$Timepoint  <- factor(st$Timepoint,  levels = test_tps)
  rownames(tab) <- NULL; rownames(st) <- NULL
  ran <- as.character(st$Timepoint)
  list(ok = TRUE, tab = tab, stats = st, declined = declined,
       baselines = c(baseline_a, baseline_b), test_tps = ran, fdr = fdr,
       n_units = st$n_units[1], unit = st$unit[1],
       ttl = sprintf("Two baselines (%s vs %s) across %d timepoint%s (n=%s %s)",
                     baseline_a, baseline_b, length(ran), if (length(ran) == 1) "" else "s",
                     st$n_units[1], st$unit[1]),
       note = local({
         ok_pct <- is.finite(st$pct_overlap)
         tail_txt <- if (length(declined))
           sprintf(" %d timepoint(s) declined and are not shown.", length(declined)) else ""
         if (!any(ok_pct))
           return(sprintf("Neither baseline called any feature significant at FDR < %.2g at any timepoint, so there is no overlap to report.%s",
                          fdr, tail_txt))
         f <- which(ok_pct)[1]; l <- rev(which(ok_pct))[1]
         if (f == l)
           return(sprintf("Overlap is %.0f%% at %s; the other timepoint(s) had nothing significant from either baseline.%s",
                          st$pct_overlap[f], ran[f], tail_txt))
         sprintf("Agreement %s from %.0f%% at %s to %.0f%% at %s.%s",
                 if (st$pct_overlap[l] >= st$pct_overlap[f]) "rises" else "falls",
                 st$pct_overlap[f], ran[f], st$pct_overlap[l], ran[l], tail_txt)
       }))
}

byod_baseline_gg <- function(obj, cols = NULL, ncol = NULL) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  df <- obj$tab
  df <- df[is.finite(df$logFC_a) & is.finite(df$logFC_b), , drop = FALSE]
  lim <- max(abs(c(df$logFC_a, df$logFC_b)), na.rm = TRUE) * 1.02
  only_a <- sprintf("only from %s", obj$baselines[1])
  only_b <- sprintf("only from %s", obj$baselines[2])
  df$agreement <- factor(df$agreement,
                         levels = c("significant from both", only_a, only_b,
                                    "significant from neither"))
  pal <- setNames(c("#D55E00", "#0072B2", "#009E73", "grey85"),
                  c("significant from both", only_a, only_b, "significant from neither"))
  lab <- obj$stats
  lab$x <- -lim * 0.94; lab$y <- lim * 0.92
  lab$txt <- sprintf("italic(R)^2 == %.4f", lab$r2)
  bg <- df[df$agreement == "significant from neither", , drop = FALSE]
  fg <- df[df$agreement != "significant from neither", , drop = FALSE]
  ggplot(mapping = aes(logFC_b, logFC_a)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                linewidth = .4, colour = "grey45") +
    geom_point(data = bg, aes(colour = agreement), size = .4, alpha = .3) +
    geom_point(data = fg, aes(colour = agreement), size = .5, alpha = .65) +
    geom_smooth(data = df, method = "lm", formula = y ~ x, se = FALSE,
                colour = "black", linewidth = .5) +
    geom_text(data = lab, aes(x, y, label = txt), parse = TRUE,
              hjust = 0, vjust = 1, size = 3.4, inherit.aes = FALSE) +
    facet_wrap(~ Timepoint, ncol = ncol %||% min(4L, nlevels(df$Timepoint))) +
    scale_colour_manual(values = pal, drop = FALSE, name = NULL) +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(title = paste(strwrap(obj$ttl, width = 66), collapse = "\n"),
         subtitle = paste(strwrap(sprintf("%s Both baselines test the same features; \"only from\" means only that contrast called it significant. Dashed line y = x, solid line the fit.",
                                          obj$note), width = 100), collapse = "\n"),
         x = bquote("Mean log"[2] * "FC (" * .(obj$baselines[2]) * " baseline)"),
         y = bquote("Mean log"[2] * "FC (" * .(obj$baselines[1]) * " baseline)")) +
    byod_theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.15),
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank()) +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1)))
}

byod_baseline_overlap_gg <- function(obj) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  st <- obj$stats
  only_a <- sprintf("only from %s", obj$baselines[1])
  only_b <- sprintf("only from %s", obj$baselines[2])
  long <- rbind(
    data.frame(Timepoint = st$Timepoint, part = only_a, n = st$only_a),
    data.frame(Timepoint = st$Timepoint, part = "shared", n = st$n_both),
    data.frame(Timepoint = st$Timepoint, part = only_b, n = st$only_b))
  long$part <- factor(long$part, levels = c(only_b, "shared", only_a))
  pal <- setNames(c("#009E73", "#4FB0A5", "#0072B2"), c(only_b, "shared", only_a))
  ggplot(long, aes(Timepoint, n, fill = part)) +
    geom_col(width = .68) +
    geom_text(data = st, aes(Timepoint, n_union, label = sprintf("%.1f%%", pct_overlap)),
              inherit.aes = FALSE, vjust = -0.6, size = 3.3, fontface = "bold") +
    scale_fill_manual(values = pal, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, .12))) +
    labs(title = paste(strwrap(sprintf("Features called significant from each baseline (FDR < %.2f)",
                                       obj$fdr), width = 66), collapse = "\n"),
         subtitle = paste(strwrap(obj$note, width = 96), collapse = "\n"),
         x = NULL, y = sprintf("features at FDR < %.2f", obj$fdr)) +
    byod_theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.15),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank())
}

byod_baseline_plot <- function(obj, cols = NULL) {
  g <- byod_baseline_gg(obj, cols)
  if (inherits(g, "plotly")) return(g)
  # ggplotly() cannot convert the plotmath axis titles used for static export.
  g <- g + labs(x = sprintf("Mean log2FC (%s baseline)", obj$baselines[2]),
                y = sprintf("Mean log2FC (%s baseline)", obj$baselines[1]))
  ggplotly(g) %>% ply()
}

byod_baseline_overlap_plot <- function(obj) {
  g <- byod_baseline_overlap_gg(obj)
  if (inherits(g, "plotly")) return(g)
  ggplotly(g) %>% ply()
}

byod_pathway_clinical <- function(d, gsva_scores, vars, method = "spearman",
                                  min_n = 8L, adjust = "BH") {
  if (is.null(gsva_scores) || !nrow(gsva_scores))
    return(list(ok = FALSE, msg = "Compute pathway activity first, then run this."))
  subj_col <- d$subject
  if (is.null(subj_col) || !subj_col %in% names(d$meta))
    return(list(ok = FALSE, msg = "Needs a declared subject, so each donor's characteristic is counted once and not once per profile."))
  vars <- intersect(vars, names(d$meta))
  if (!length(vars)) return(list(ok = FALSE, msg = "None of the selected variables is in the metadata."))
  keep <- intersect(colnames(gsva_scores), rownames(d$meta))
  if (length(keep) < 4) return(list(ok = FALSE, msg = "Pathway scores and metadata share too few samples."))
  S <- gsva_scores[, keep, drop = FALSE]
  subj <- as.character(d$meta[keep, subj_col])
  donors <- unique(subj)

  Sd <- vapply(donors, function(dn) {
    j <- which(subj == dn)
    if (length(j) == 1L) S[, j] else rowMeans(S[, j, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(S)))
  dimnames(Sd) <- list(rownames(S), donors)

  out <- list(); skipped <- character()
  for (v in vars) {
    cv <- d$meta[keep, v]
    isnum <- is.numeric(cv)
    consistent <- vapply(donors, function(dn) length(unique(cv[subj == dn & !is.na(cv)])) <= 1L, logical(1))
    if (!all(consistent)) { skipped <- c(skipped, v); next }
    val <- vapply(donors, function(dn) {
      z <- cv[subj == dn]; z <- z[!is.na(z)]
      if (!length(z)) return(if (isnum) NA_real_ else NA_character_)
      if (isnum) as.numeric(z[1]) else as.character(z[1])
    }, if (isnum) numeric(1) else character(1))
    if (isnum) {
      ok <- is.finite(val)
      if (sum(ok) < min_n) next
      x <- val[ok]; M <- Sd[, ok, drop = FALSE]; n <- length(x)
      r <- suppressWarnings(as.numeric(stats::cor(t(M), x, method = method)))
      tv <- r * sqrt((n - 2) / (1 - r^2))
      p <- 2 * stats::pt(-abs(tv), n - 2)
      out[[v]] <- data.frame(pathway = rownames(M), variable = v, type = "continuous",
                             estimate = r, p = p, n = n, tested = NA_character_, reference = NA_character_,
                             stringsAsFactors = FALSE)
    } else {
      observed <- unique(val[!is.na(val)])
      preset <- d$clinical_contrasts[[v]]
      lv <- if (!is.null(preset)) preset else if (is.factor(cv)) rev(levels(cv)[levels(cv) %in% observed])
            else rev(sort(observed, method = "radix"))
      if (length(lv) != 2L) next                       # k-level factors need a different test
      if (!setequal(lv, observed)) next
      g <- factor(val, levels = lv)
      ok <- !is.na(g)
      if (min(table(g[ok])) < 3L) next
      M <- Sd[, ok, drop = FALSE]; gg <- droplevels(g[ok])
      res <- t(apply(M, 1, function(y) {
        tt <- tryCatch(stats::t.test(y ~ gg), error = function(e) NULL)
        if (is.null(tt)) c(NA_real_, NA_real_) else c(unname(diff(rev(tt$estimate))), tt$p.value)
      }))
      out[[v]] <- data.frame(pathway = rownames(M), variable = v, type = "two-group",
                             estimate = res[, 1], p = res[, 2], n = sum(ok),
                             tested = lv[1], reference = lv[2],
                             stringsAsFactors = FALSE)
    }
  }
  if (!length(out))
    return(list(ok = FALSE, msg = sprintf("No selected variable had at least %d subjects with a value (continuous) or two groups of at least three (categorical).", min_n)))
  tab <- do.call(rbind, out)
  tab$fdr <- stats::p.adjust(tab$p, adjust)            # one family across the whole grid
  tab <- tab[order(tab$fdr, -abs(tab$estimate)), ]
  rownames(tab) <- NULL
  list(ok = TRUE, tab = tab, n_units = length(donors), subj_col = subj_col,
       skipped_inconsistent_variables = skipped,
       n_sig = sum(tab$fdr < 0.05, na.rm = TRUE),
       ttl = sprintf("Pathway activity vs subject characteristics (n=%d %ss)", length(donors), subj_col),
       note = sprintf("%d of %d pathway x variable tests at FDR < 0.05. Scores are collapsed to one value per %s before testing, and the FDR is taken across the whole grid.",
                      sum(tab$fdr < 0.05, na.rm = TRUE), nrow(tab), subj_col))
}

byod_pathway_clinical_gg <- function(obj, max_path = 40L) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  tab <- obj$tab

  ord <- names(sort(tapply(abs(tab$estimate), tab$pathway, max, na.rm = TRUE), decreasing = TRUE))
  keep <- head(ord, max_path)
  tab <- tab[tab$pathway %in% keep, , drop = FALSE]
  tab$pathway  <- factor(gsub("_", " ", tab$pathway), levels = rev(gsub("_", " ", keep)))
  tab$variable <- factor(tab$variable, levels = sort(unique(tab$variable)))
  tab$size <- -log10(pmax(tab$fdr, 1e-12))
  lim <- max(abs(tab$estimate), na.rm = TRUE)

  sig <- tab[is.finite(tab$fdr) & tab$fdr < 0.05, , drop = FALSE]
  ggplot(tab, aes(variable, pathway)) +
    geom_point(aes(colour = estimate, size = size)) +
    { if (nrow(sig)) geom_point(data = sig, aes(size = size), shape = 21,
                                colour = "grey15", stroke = .45, fill = NA) } +
    scale_colour_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                           midpoint = 0, limits = c(-lim, lim), name = "effect") +
    scale_size_continuous(range = c(.6, 5.5), name = expression(-log[10]~FDR)) +
    labs(title = paste(strwrap(obj$ttl, width = 58), collapse = "\n"),
         subtitle = paste(strwrap(sprintf("%s Ringed cells pass FDR < 0.05.", obj$note),
                                  width = 88), collapse = "\n"),
         x = NULL, y = NULL) +
    byod_theme_bw(base_size = 10) +
    theme(axis.text.x   = element_text(angle = 40, hjust = 1, size = 8.5),
          axis.text.y   = element_text(size = 7.5),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.15),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(linewidth = .2, colour = "grey92"),
          legend.key.height = unit(.5, "cm"))
}

byod_pathway_clinical_plot <- function(obj, max_path = 40L) {
  g <- byod_pathway_clinical_gg(obj, max_path)
  if (inherits(g, "plotly")) return(g)
  ggplotly(g) %>% ply()
}

byod_multi_contrast <- function(d, reference, test_tps, fdr = 0.05, lfc = 0,
                                subject_col = NULL, paired = TRUE) {
  lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
  if (!reference %in% lv)
    return(list(ok = FALSE, msg = sprintf("Control timepoint '%s' is not in this study.", reference)))
  test_tps <- setdiff(unique(as.character(test_tps)), reference)
  miss <- setdiff(test_tps, lv)
  if (length(miss))
    return(list(ok = FALSE, msg = sprintf("Not in this study: %s.", paste(miss, collapse = ", "))))
  if (!length(test_tps))
    return(list(ok = FALSE, msg = "Pick at least one timepoint to compare against the control."))
  test_tps <- lv[lv %in% test_tps]                    # keep the declared order
  sc <- subject_col %||% d$subject

  res <- list(); sets <- list(); stats <- list(); declined <- list()
  for (tp in test_tps) {
    de <- tryCatch(byod_run_limma(d, "Timepoint", reference, tp,
                                  subject_col = sc, paired = paired && !is.null(sc)),
                   error = function(e) e)
    if (inherits(de, "error")) { declined[[tp]] <- conditionMessage(de); next }
    if (!is.null(attr(de, "declined"))) { declined[[tp]] <- attr(de, "declined"); next }
    de$Feature <- de$Feature %||% rownames(de)
    keep <- is.finite(de$adj.P.Val) & de$adj.P.Val < fdr & abs(de$logFC) >= lfc
    res[[tp]]   <- de
    sets[[tp]]  <- de$Feature[keep]
    stats[[tp]] <- data.frame(Timepoint = tp, n_tested = nrow(de), n_sig = sum(keep),
                              n_up = sum(keep & de$logFC > 0), n_down = sum(keep & de$logFC < 0),
                              n_units = attr(de, "n_units"), unit = attr(de, "unit"),
                              nested = isTRUE(attr(de, "nested")), stringsAsFactors = FALSE)
  }
  if (!length(res))
    return(list(ok = FALSE,
                msg = paste("No timepoint could be contrasted against the control.",
                            paste(sprintf("%s: %s", names(declined), unlist(declined)),
                                  collapse = "; "))))
  st <- do.call(rbind, stats); rownames(st) <- NULL
  ran <- st$Timepoint
  st$Timepoint <- factor(st$Timepoint, levels = ran)

  all_sig <- unique(unlist(sets))
  core <- if (length(sets) > 1) Reduce(intersect, sets) else all_sig
  only <- vapply(names(sets), function(k)
    length(setdiff(sets[[k]], unique(unlist(sets[names(sets) != k])))), integer(1))

  list(ok = TRUE, results = res, sets = sets, stats = st, declined = declined,
       reference = reference, test_tps = ran, fdr = fdr, lfc = lfc,
       n_any = length(all_sig), n_core = length(core), core = core, n_only = only,
       n_units = st$n_units[1], unit = st$unit[1], nested = st$nested[1],
       ttl = sprintf("Every timepoint against %s (n=%s %s%s)", reference,
                     st$n_units[1], st$unit[1], if (st$nested[1]) ", nested" else ""),
       note = sprintf("%s features are significant at FDR < %.2f in at least one timepoint; %s in all %d.%s",
                      format(length(all_sig), big.mark = ","), fdr,
                      format(length(core), big.mark = ","), length(sets),
                      if (length(declined)) sprintf(" %d timepoint(s) declined and are not shown.",
                                                    length(declined)) else ""))
}

byod_multi_volcano_gg <- function(obj, ncol = NULL, label_n = 0L) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  df <- do.call(rbind, lapply(names(obj$results), function(tp) {
    de <- obj$results[[tp]]
    data.frame(Timepoint = tp, Feature = de$Feature, logFC = de$logFC,
               q = de$adj.P.Val, stringsAsFactors = FALSE)
  }))
  df <- df[is.finite(df$logFC) & is.finite(df$q), , drop = FALSE]
  df$Timepoint <- factor(df$Timepoint, levels = obj$test_tps)
  df$sig <- df$q < obj$fdr & abs(df$logFC) >= obj$lfc
  lim <- max(abs(df$logFC), na.rm = TRUE) * 1.02
  lab <- obj$stats
  lab$txt <- sprintf("%s up / %s down", format(lab$n_up, big.mark = ","),
                     format(lab$n_down, big.mark = ","))
  g <- ggplot(df, aes(logFC, -log10(pmax(q, 1e-300)))) +
    geom_hline(yintercept = -log10(obj$fdr), linetype = "dashed",
               linewidth = .4, colour = "grey45") +
    geom_vline(xintercept = 0, linewidth = .25, colour = "grey80") +
    geom_point(aes(colour = sig), size = .4, alpha = .45) +
    geom_text(data = lab, aes(x = 0, y = Inf, label = txt), inherit.aes = FALSE,
              hjust = .5, vjust = 1.5, size = 3.1, colour = "grey20") +

    facet_wrap(~ Timepoint,
               ncol = ncol %||% (if (length(obj$test_tps) <= 6L) length(obj$test_tps)
                                 else ceiling(length(obj$test_tps) / 2))) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey78"), guide = "none") +
    scale_x_continuous(limits = c(-lim, lim)) +
    labs(title = paste(strwrap(obj$ttl, width = 66), collapse = "\n"),
         subtitle = paste(strwrap(obj$note, width = 96), collapse = "\n"),
         x = bquote("log"[2] * " fold change vs " * .(obj$reference)),
         y = bquote(-log[10] * italic(q))) +
    byod_theme_bw(base_size = 11) +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.15),
          strip.background = element_rect(fill = "grey94", colour = NA),
          strip.text = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank())
  if (label_n > 0L) {
    top <- do.call(rbind, lapply(split(df[df$sig, , drop = FALSE], df$Timepoint[df$sig]),
                                 function(x) head(x[order(x$q), ], label_n)))
    if (!is.null(top) && nrow(top))
      g <- g + ggrepel::geom_text_repel(data = top, aes(label = Feature), seed = trace_seed(), size = 2.6,
                                        max.overlaps = Inf, min.segment.length = 0,
                                        segment.size = .25, colour = "grey15")
  }
  g
}

byod_multi_upset_gg <- function(obj, max_sets = 12L) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  sets <- obj$sets
  ks <- names(sets)
  if (length(ks) < 2)
    return(emptyplot_gg("An UpSet needs at least two timepoints with a result."))
  all_f <- unique(unlist(sets))
  if (!length(all_f))
    return(emptyplot_gg(sprintf("No feature is significant at FDR < %.2f in any timepoint.", obj$fdr)))

  memb <- vapply(ks, function(k) all_f %in% sets[[k]], logical(length(all_f)))
  if (is.null(dim(memb))) memb <- matrix(memb, nrow = length(all_f), dimnames = list(NULL, ks))
  key <- apply(memb, 1, function(z) paste(as.integer(z), collapse = ""))
  tb  <- sort(table(key), decreasing = TRUE)
  tb  <- head(tb[names(tb) != paste(rep(0, length(ks)), collapse = "")], max_sets)
  if (!length(tb)) return(emptyplot_gg("No intersection to show."))

  pat <- names(tb)
  bars <- data.frame(i = seq_along(pat), n = as.integer(tb), stringsAsFactors = FALSE)

  band <- max(bars$n) * 0.42
  step <- band / (length(ks) + 1)
  dots <- do.call(rbind, lapply(seq_along(pat), function(i) {
    on <- strsplit(pat[i], "")[[1]] == "1"
    data.frame(i = i, set = ks, y = -step * seq_along(ks), on = on, stringsAsFactors = FALSE)
  }))
  seg <- do.call(rbind, lapply(seq_along(pat), function(i) {
    d <- dots[dots$i == i & dots$on, , drop = FALSE]
    if (nrow(d) < 2) return(NULL)
    data.frame(i = i, ymin = min(d$y), ymax = max(d$y))
  }))
  set_lab <- data.frame(set = ks, y = -step * seq_along(ks),
                        size = vapply(ks, function(k) length(sets[[k]]), integer(1)))

  ggplot() +
    geom_col(data = bars, aes(i, n), fill = "#2C6E9B", width = .68) +
    geom_text(data = head(bars, 3), aes(i, n, label = format(n, big.mark = ",")),
              vjust = -0.45, size = 2.8, colour = "grey20") +
    geom_hline(yintercept = 0, linewidth = .3, colour = "grey70") +
    { if (!is.null(seg)) geom_segment(data = seg, aes(x = i, xend = i, y = ymin, yend = ymax),
                                      linewidth = .7, colour = "#2C6E9B") } +
    geom_point(data = dots, aes(i, y, colour = on), size = 2.6) +
    geom_text(data = set_lab, aes(x = 0.2, y = y, label = sprintf("%s  (%s)", set,
                                  format(size, big.mark = ","))),
              hjust = 1, size = 3, colour = "grey20") +
    scale_colour_manual(values = c(`TRUE` = "#2C6E9B", `FALSE` = "grey86"), guide = "none") +
    scale_x_continuous(limits = c(-max(7, max(nchar(ks)) * 0.7 + 2), length(pat) + .6),
                       expand = expansion(mult = c(0, .01))) +

    scale_y_continuous(expand = expansion(mult = c(.06, .12)),
                       breaks = function(z) { b <- scales::extended_breaks()(c(0, max(z)))
                                              b[b >= 0] }) +
    labs(title = paste(strwrap(sprintf("Shared and timepoint-specific features vs %s",
                                       obj$reference), width = 66), collapse = "\n"),
         subtitle = paste(strwrap(sprintf("%s Bars are intersection sizes; a filled column shows which timepoints that intersection belongs to, and the number after each timepoint is its own set size.",
                                          obj$note), width = 100), collapse = "\n"),
         x = NULL, y = "features in the intersection") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.15),
          axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y  = element_text(size = 8.5),
          panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
}

byod_multi_volcano_one_gg <- function(obj, tp, label_n = 0L) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  if (!tp %in% names(obj$results))
    return(emptyplot_gg(sprintf("'%s' has no result to draw.", tp)))
  all_fc <- unlist(lapply(obj$results, function(x) x$logFC))
  lim <- max(abs(all_fc[is.finite(all_fc)]), na.rm = TRUE) * 1.02
  de <- obj$results[[tp]]
  df <- data.frame(Feature = de$Feature, logFC = de$logFC, q = de$adj.P.Val,
                   stringsAsFactors = FALSE)
  df <- df[is.finite(df$logFC) & is.finite(df$q), , drop = FALSE]
  df$sig <- df$q < obj$fdr & abs(df$logFC) >= obj$lfc
  s <- obj$stats[obj$stats$Timepoint == tp, ]
  g <- ggplot(df, aes(logFC, -log10(pmax(q, 1e-300)))) +
    geom_hline(yintercept = -log10(obj$fdr), linetype = "dashed",
               linewidth = .4, colour = "grey45") +
    geom_vline(xintercept = 0, linewidth = .25, colour = "grey80") +
    geom_point(aes(colour = sig), size = .5, alpha = .5) +
    scale_colour_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey78"), guide = "none") +
    scale_x_continuous(limits = c(-lim, lim)) +
    labs(title = sprintf("%s vs %s", tp, obj$reference),
         subtitle = sprintf("%s up / %s down at FDR < %.2f; n = %s %s%s",
                            format(s$n_up, big.mark = ","), format(s$n_down, big.mark = ","),
                            obj$fdr, s$n_units, s$unit,
                            if (isTRUE(s$nested)) ", nested" else ""),
         x = bquote("log"[2] * " fold change vs " * .(obj$reference)),
         y = bquote(-log[10] * italic(q))) +
    byod_theme_bw(base_size = 12) +
    theme(plot.title    = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey30"),
          panel.grid.minor = element_blank())
  if (label_n > 0L) {
    top <- head(df[df$sig, ][order(df$q[df$sig]), ], label_n)
    if (nrow(top))
      g <- g + ggrepel::geom_text_repel(data = top, aes(label = Feature), seed = trace_seed(), size = 2.8,
                                        max.overlaps = Inf, min.segment.length = 0,
                                        segment.size = .25, colour = "grey15")
  }
  g
}

byod_multi_volcano_export <- function(obj, dir, fmts = c("png", "pdf"),
                                      w_one = 5.2, h_one = 4.6, dpi = 300) {
  if (!isTRUE(obj$ok)) return(character(0))
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  safe <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
  out <- character(0)
  for (tp in names(obj$results)) for (fmt in fmts) {
    f <- file.path(dir, sprintf("volcano_%s_vs_%s.%s", safe(tp), safe(obj$reference), fmt))
    ok <- tryCatch({ byod_save_figure(byod_multi_volcano_one_gg(obj, tp), f, fmt,
                                    w_one, h_one, dpi); TRUE },
                   error = function(e) FALSE)
    if (ok) out <- c(out, f)
  }
  n <- length(obj$test_tps)
  ncol <- if (n <= 6L) n else ceiling(n / 2)
  for (fmt in fmts) {
    f <- file.path(dir, sprintf("volcano_ALL_vs_%s.%s", safe(obj$reference), fmt))
    ok <- tryCatch({ byod_save_figure(byod_multi_volcano_gg(obj), f, fmt,
                                    min(24, 2.9 * ncol + 1.2),
                                    3.4 * ceiling(n / ncol) + 1.2, dpi); TRUE },
                   error = function(e) FALSE)
    if (ok) out <- c(out, f)
  }
  out
}

byod_multi_volcano_plot <- function(obj, ncol = NULL) {
  g <- byod_multi_volcano_gg(obj, ncol)
  if (inherits(g, "plotly")) return(g)
  # ggplotly() cannot convert the plotmath axis titles used for static export.
  g <- g + labs(x = paste("log2 fold change vs", obj$reference), y = "-log10 q")
  ggplotly(g) %>% ply()
}
byod_multi_upset_plot <- function(obj, max_sets = 25L) {
  g <- byod_multi_upset_gg(obj, max_sets)
  if (inherits(g, "plotly")) return(g)
  ggplotly(g) %>% ply()
}

byod_response_composition <- function(d, reference, test_tps, group_var = NULL,
                                      fdr = 0.05, lfc = 0, subject_col = NULL, paired = TRUE) {
  lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
  if (!reference %in% lv)
    return(list(ok = FALSE, msg = sprintf("Baseline '%s' is not in this study.", reference)))
  test_tps <- lv[lv %in% setdiff(unique(as.character(test_tps)), reference)]
  if (!length(test_tps))
    return(list(ok = FALSE, msg = "Pick at least one timepoint to compare against the baseline."))
  sc <- subject_col %||% d$subject
  n_feat <- nrow(d$expr)

  arms <- list(all = seq_len(ncol(d$expr)))
  if (!is.null(group_var) && nzchar(group_var) && group_var %in% names(d$meta)) {
    g <- as.character(d$meta[[group_var]]); ug <- unique(g[!is.na(g)])
    if (length(ug) != 2L)
      return(list(ok = FALSE, msg = sprintf("'%s' has %d level(s); this view compares exactly two arms, or none at all.",
                                            group_var, length(ug))))
    arms <- setNames(lapply(ug, function(k) which(g == k)), ug)
  }

  sets <- list(); declined <- list()
  for (a in names(arms)) {
    da <- d; da$expr <- d$expr[, arms[[a]], drop = FALSE]; da$meta <- d$meta[arms[[a]], , drop = FALSE]
    for (tp in test_tps) {
      de <- tryCatch(byod_run_limma(da, "Timepoint", reference, tp,
                                    subject_col = sc, paired = paired && !is.null(sc)),
                     error = function(e) e)
      k <- paste(a, tp, sep = "|")
      if (inherits(de, "error")) { declined[[k]] <- conditionMessage(de); next }
      if (!is.null(attr(de, "declined"))) { declined[[k]] <- attr(de, "declined"); next }
      f <- de$Feature %||% rownames(de)
      sets[[k]] <- f[is.finite(de$adj.P.Val) & de$adj.P.Val < fdr & abs(de$logFC) >= lfc]
    }
  }
  ok_tps <- test_tps[vapply(test_tps, function(tp)
    all(paste(names(arms), tp, sep = "|") %in% names(sets)), logical(1))]
  if (!length(ok_tps))
    return(list(ok = FALSE, msg = paste("No timepoint could be summarised.",
                paste(sprintf("%s: %s", names(declined), unlist(declined)), collapse = "; "))))

  two <- length(arms) == 2L
  rows <- lapply(ok_tps, function(tp) {
    if (!two) {
      s <- sets[[paste("all", tp, sep = "|")]]
      data.frame(Timepoint = tp,
                 category = c("changed", "unchanged"),
                 n = c(length(s), n_feat - length(s)), stringsAsFactors = FALSE)
    } else {
      a <- sets[[paste(names(arms)[1], tp, sep = "|")]]
      b <- sets[[paste(names(arms)[2], tp, sep = "|")]]
      both <- length(intersect(a, b))
      data.frame(Timepoint = tp,
                 category = c("changed in both",
                              sprintf("only %s", names(arms)[1]),
                              sprintf("only %s", names(arms)[2]),
                              "unchanged"),
                 n = c(both, length(setdiff(a, b)), length(setdiff(b, a)),
                       n_feat - length(union(a, b))), stringsAsFactors = FALSE)
    }
  })
  tab <- do.call(rbind, rows)
  tab$Timepoint <- factor(tab$Timepoint, levels = ok_tps)
  tab$pct <- 100 * tab$n / n_feat
  ord <- if (two) c(sprintf("only %s", names(arms)[1]), sprintf("only %s", names(arms)[2]),
                    "changed in both", "unchanged")
         else c("changed", "unchanged")
  tab$category <- factor(tab$category, levels = ord)

  moved <- tapply(tab$n[tab$category != "unchanged"],
                  droplevels(tab$Timepoint[tab$category != "unchanged"]), sum)
  list(ok = TRUE, tab = tab, n_features = n_feat, reference = reference,
       test_tps = ok_tps, arms = names(arms), two_arm = two,
       group_var = if (two) group_var else NULL, fdr = fdr, lfc = lfc, declined = declined,
       moved_pct = 100 * moved / n_feat,
       ttl = sprintf("Share of the layer that has moved from %s%s", reference,
                     if (two) sprintf(", by %s", group_var) else ""),
       note = local({
         tail_txt <- if (length(declined)) sprintf(" %d arm-timepoint(s) declined.",
                                                   length(declined)) else ""
         if (length(ok_tps) == 1L)
           sprintf("%.1f%% of %s features have moved by %s.%s",
                   100 * moved[1] / n_feat, format(n_feat, big.mark = ","), ok_tps[1], tail_txt)
         else
           sprintf("%.1f%% of %s features at %s, %s to %.1f%% at %s.%s",
                   100 * moved[1] / n_feat, format(n_feat, big.mark = ","), ok_tps[1],
                   if (moved[length(moved)] >= moved[1]) "rising" else "falling",
                   100 * moved[length(moved)] / n_feat, ok_tps[length(ok_tps)], tail_txt)
       }))
}

byod_response_composition_gg <- function(obj, as_percent = TRUE) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  tab <- obj$tab
  tab$y <- if (as_percent) tab$pct else tab$n
  pal <- if (obj$two_arm)
    setNames(c("#B59FCF", "#8CD3C3", "#D9BE6C", "grey88"), levels(tab$category))
  else setNames(c("#D55E00", "grey88"), levels(tab$category))

  single <- nlevels(droplevels(tab$Timepoint)) < 2L
  ggplot(tab, aes(as.integer(Timepoint), y, fill = category)) +
    { if (single) geom_col(width = .55, colour = "white", linewidth = .3)
      else geom_area(colour = "white", linewidth = .3) } +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_continuous(breaks = seq_along(levels(tab$Timepoint)),
                       labels = levels(tab$Timepoint),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(expand = expansion(mult = c(0, .02)),
                       labels = if (as_percent) function(z) paste0(z, "%") else
                                function(z) format(z, big.mark = ",")) +
    labs(title = paste(strwrap(obj$ttl, width = 62), collapse = "\n"),
         subtitle = paste(strwrap(obj$note, width = 92), collapse = "\n"),
         x = NULL, y = if (as_percent) "share of measured features" else "features") +
    byod_theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.15),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())
}

byod_split_variable <- function(values, how = "median", probs = c(.25, .75)) {
  v <- suppressWarnings(as.numeric(as.character(values)))
  if (all(is.na(v))) {
    f <- factor(as.character(values))
    return(list(ok = nlevels(f) >= 2, f = f, label = "as given"))
  }
  cut_at <- switch(how,
    median = stats::median(v, na.rm = TRUE),
    mean   = mean(v, na.rm = TRUE),
    q25    = stats::quantile(v, .25, na.rm = TRUE),
    q75    = stats::quantile(v, .75, na.rm = TRUE),
    stats::median(v, na.rm = TRUE))
  lab <- switch(how, median = "median", mean = "mean", q25 = "25th percentile",
                q75 = "75th percentile", "median")
  f <- factor(ifelse(v > cut_at, "High", "Low"), levels = c("Low", "High"))
  f[is.na(v)] <- NA
  list(ok = nlevels(droplevels(f)) == 2L, f = f,
       label = sprintf("split at the %s (%.4g)", lab, cut_at))
}

spearman_p_rowwise <- function(M, x) {
  n <- length(x)
  if (n < 3L || !nrow(M)) return(rep(NA_real_, nrow(M)))
  rx <- rank(x); rx <- rx - mean(rx)
  sx <- sqrt(sum(rx^2))
  RM <- if (requireNamespace("matrixStats", quietly = TRUE))
          matrixStats::rowRanks(M, ties.method = "average")
        else t(apply(M, 1, rank))
  RM <- RM - rowMeans(RM)
  sy <- sqrt(rowSums(RM^2))
  r  <- as.numeric(RM %*% rx) / pmax(sx * sy, .Machine$double.eps)
  r[!is.finite(r)] <- 0
  r  <- pmax(pmin(r, 1 - 1e-12), -1 + 1e-12)
  tt <- r * sqrt((n - 2) / (1 - r^2))
  p  <- 2 * stats::pt(-abs(tt), df = n - 2)
  p[sy == 0] <- NA_real_          # a feature that is constant has no correlation
  p
}

# Adjust across features within each variable-timepoint test, not across the whole screen.
byod_variable_screen <- function(d, vars, test_tps = NULL, splits = NULL,
                                 fdr = 0.05, min_per_group = 3L, min_n_cont = 8L) {
  lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
  test_tps <- if (is.null(test_tps)) lv else lv[lv %in% as.character(test_tps)]
  if (!length(test_tps)) return(list(ok = FALSE, msg = "No timepoint selected."))
  vars <- intersect(vars, names(d$meta))
  if (!length(vars)) return(list(ok = FALSE, msg = "None of the selected variables is in the metadata."))
  n_feat <- nrow(d$expr)

  .dm_cache <- new.env(parent = emptyenv())
  donor_matrix <- function(tp) {
    key <- paste0("tp_", as.character(tp))
    hit <- .dm_cache[[key]]
    if (!is.null(hit)) return(if (identical(hit, NA)) NULL else hit)
    keep <- which(as.character(d$meta$Timepoint) == tp)
    if (!length(keep)) { assign(key, NA, envir = .dm_cache); return(NULL) }
    X <- d$expr[, keep, drop = FALSE]
    sc <- d$subject
    if (is.null(sc) || !sc %in% names(d$meta)) {
      colnames(X) <- rownames(d$meta)[keep]; assign(key, X, envir = .dm_cache); return(X)
    }
    g <- as.character(d$meta[[sc]][keep])

    S <- t(rowsum(t(X), g, na.rm = TRUE))
    N <- t(rowsum(t(!is.na(X)) + 0, g))
    M <- S / pmax(N, 1); M[N == 0] <- NA_real_
    rownames(M) <- rownames(X)
    assign(key, M, envir = .dm_cache); M
  }
  one_per_subject <- function(v) {
    sc <- d$subject
    if (is.null(sc) || !sc %in% names(d$meta)) return(setNames(d$meta[[v]], rownames(d$meta)))
    i <- !duplicated(as.character(d$meta[[sc]]))
    setNames(d$meta[[v]][i], as.character(d$meta[[sc]][i]))
  }
  is_continuous <- function(v) {
    if (!is.null(splits) && !is.null(splits[[v]])) return(FALSE)     # user asked for a split
    if (!is.null(d$clinical_num) && v %in% d$clinical_num) return(TRUE)
    if (!is.null(d$clinical_cat) && v %in% d$clinical_cat) return(FALSE)
    x <- d$meta[[v]]
    is.numeric(x) && length(unique(x[!is.na(x)])) > 2L
  }

  rows <- list(); skipped <- list(); order_i <- 0L
  for (v in vars) {
    order_i <- order_i + 1L
    if (is_continuous(v)) {
      cmp <- sprintf("%s: Spearman (continuous)", v)
      xs  <- suppressWarnings(as.numeric(one_per_subject(v)))
      names(xs) <- names(one_per_subject(v))
      for (tp in test_tps) {
        M <- donor_matrix(tp)
        x <- if (is.null(M)) numeric(0) else xs[colnames(M)]
        good <- is.finite(x)
        base <- data.frame(order = order_i, variable = v, route = "continuous",
                           comparison = cmp, Timepoint = tp, n_sig = NA_integer_,
                           pct = NA_real_, n1 = sum(good), n2 = NA_integer_,
                           split = "not split - correlated on its own scale",
                           note = "", stringsAsFactors = FALSE)
        if (is.null(M) || sum(good) < min_n_cont) {
          base$note <- sprintf("needs >= %d subjects with a value", min_n_cont)
          rows[[length(rows) + 1L]] <- base; next
        }
        pv <- spearman_p_rowwise(M[, good, drop = FALSE], x[good])
        ns <- sum(stats::p.adjust(pv, "BH") < fdr, na.rm = TRUE)
        base$n_sig <- ns; base$pct <- 100 * ns / n_feat
        rows[[length(rows) + 1L]] <- base
      }
      next
    }

    how <- (splits %||% list())[[v]] %||% "levels"
    sp <- if (identical(how, "levels")) list(f = droplevels(factor(as.character(d$meta[[v]]))),
                                             label = "levels as recorded")
          else byod_split_variable(d$meta[[v]], how)
    f <- sp$f
    lvls <- levels(droplevels(f))
    if (length(lvls) < 2L) { skipped[[v]] <- "fewer than two groups"; next }
    pairs <- if (length(lvls) == 2L) list(lvls) else lapply(lvls[-1], function(z) c(lvls[1], z))
    for (pr in pairs) {
      cmp <- sprintf("%s: %s vs %s", v, pr[2], pr[1])
      for (tp in test_tps) {
        keep <- which(as.character(d$meta$Timepoint) == tp & !is.na(f) & f %in% pr)
        gg <- droplevels(f[keep])
        n1 <- sum(gg == pr[1]); n2 <- sum(gg == pr[2])
        base <- data.frame(order = order_i, variable = v, route = "categorical",
                           comparison = cmp, Timepoint = tp, n_sig = NA_integer_,
                           pct = NA_real_, n1 = n1, n2 = n2, split = sp$label,
                           note = "", stringsAsFactors = FALSE)
        if (min(n1, n2) < min_per_group) {
          base$note <- sprintf("needs >= %d per group", min_per_group)
          rows[[length(rows) + 1L]] <- base; next
        }
        dsub <- d; dsub$expr <- d$expr[, keep, drop = FALSE]
        dsub$meta <- d$meta[keep, , drop = FALSE]; dsub$meta[[v]] <- gg
        de <- tryCatch(byod_run_limma(dsub, v, pr[1], pr[2], subject_col = NULL, paired = FALSE),
                       error = function(e) e)
        bad <- if (inherits(de, "error")) conditionMessage(de)
               else if (!is.null(attr(de, "declined"))) attr(de, "declined") else NA_character_
        if (is.na(bad)) {
          ns <- sum(de$adj.P.Val < fdr, na.rm = TRUE)
          base$n_sig <- ns; base$pct <- 100 * ns / n_feat
        } else base$note <- substr(bad, 1, 90)
        rows[[length(rows) + 1L]] <- base
      }
    }
  }
  if (!length(rows))
    return(list(ok = FALSE, msg = "No variable produced an estimable comparison at any timepoint."))
  tab <- do.call(rbind, rows)
  tab$Timepoint <- factor(tab$Timepoint, levels = test_tps)
  ran <- sum(!is.na(tab$n_sig))
  n_cont <- length(unique(tab$variable[tab$route == "continuous"]))
  n_cat  <- length(unique(tab$variable[tab$route == "categorical"]))
  list(ok = TRUE, tab = tab, skipped = skipped, fdr = fdr, n_features = n_feat,
       test_tps = test_tps, n_comparisons = length(unique(tab$comparison)),
       n_continuous = n_cont, n_categorical = n_cat,
       unit = d$subject %||% "profile",
       ttl = sprintf("%d comparison%s across %d timepoint%s",
                     length(unique(tab$comparison)),
                     if (length(unique(tab$comparison)) == 1) "" else "s",
                     length(test_tps), if (length(test_tps) == 1) "" else "s"),
       note = sprintf("%d of %d cells were estimable. %d continuous variable%s correlated on its own scale and %d categorical variable%s contrasted between levels, each inside one timepoint at %s resolution.%s",
                      ran, nrow(tab), n_cont, if (n_cont == 1) "" else "s",
                      n_cat, if (n_cat == 1) "" else "s", d$subject %||% "profile",
                      if (length(skipped)) sprintf(" Skipped: %s.",
                        paste(sprintf("%s (%s)", names(skipped), unlist(skipped)), collapse = "; ")) else ""))
}

byod_variable_screen_gg <- function(obj, value = c("count", "percent"), max_rows = 40L) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  value <- match.arg(value)
  tab <- obj$tab
  tab$val <- if (value == "percent") tab$pct else tab$n_sig

  ord  <- unique(tab$comparison[order(tab$order)])
  keep <- head(ord, max_rows)
  tab <- tab[tab$comparison %in% keep, , drop = FALSE]
  tab$comparison <- factor(tab$comparison, levels = rev(keep))

  ggplot(tab, aes(Timepoint, comparison, fill = val)) +
    geom_tile(colour = "white", linewidth = .6) +
    geom_text(aes(label = ifelse(is.na(val), "N/A",
                                 if (value == "percent") sprintf("%.1f", val)
                                 else format(val, big.mark = ","))),
              size = 2.6, colour = "grey15") +
    scale_fill_gradient(low = "#EAF2F8", high = "#B2182B", na.value = "grey93",
                        name = if (value == "percent") "% of features" else "features") +
    labs(title = paste(strwrap(obj$ttl, width = 62), collapse = "\n"),
         subtitle = paste(strwrap(obj$note, width = 96), collapse = "\n"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey30", lineheight = 1.15),
          axis.text.y   = element_text(size = 8),
          panel.grid    = element_blank(),
          legend.key.height = unit(.5, "cm"))
}

byod_response_composition_plot <- function(obj, as_percent = TRUE) {
  g <- byod_response_composition_gg(obj, as_percent)
  if (inherits(g, "plotly")) return(g); ggplotly(g) %>% ply()
}
byod_variable_screen_plot <- function(obj, value = "count") {
  g <- byod_variable_screen_gg(obj, value)
  if (inherits(g, "plotly")) return(g); ggplotly(g) %>% ply()
}

byod_order_candidates <- function(d, max_levels = 40L) {
  if (is.null(d$meta)) return(character(0))
  v <- vapply(names(d$meta), function(k) {
    z <- d$meta[[k]]; u <- unique(as.character(z[!is.na(z)]))
    length(u) > 1L && length(u) <= max_levels
  }, logical(1))
  out <- names(d$meta)[v]

  c(intersect("Timepoint", out), setdiff(out, c("Timepoint", "sample_id")))
}

byod_heat_order <- function(M, meta, by = "cluster", cluster_within = TRUE) {
  cn <- colnames(M)
  none <- list(order = seq_along(cn), group = NULL, by = "as given",
               label = "columns in the order they were supplied")
  if (is.null(cn) || !length(cn)) return(none)

  if (identical(by, "none")) return(none)

  if (identical(by, "cluster")) {
    if (ncol(M) < 3) return(none)
    Z <- M; Z[!is.finite(Z)] <- 0
    d1 <- tryCatch(stats::dist(t(Z)), error = function(e) NULL)
    if (is.null(d1)) return(none)
    hc <- tryCatch(stats::hclust(d1, method = "ward.D2"), error = function(e) NULL)
    if (is.null(hc)) return(none)
    return(list(order = hc$order, group = NULL, by = "cluster", hclust = hc,
                label = "columns clustered on their own profiles (Ward, Euclidean)"))
  }

  if (is.null(meta) || !by %in% names(meta)) return(none)
  g <- meta[cn, by]
  if (is.factor(g)) g <- droplevels(g) else g <- factor(as.character(g))
  if (nlevels(g) < 2L) return(none)
  ord <- order(g, seq_along(cn))

  if (isTRUE(cluster_within) && ncol(M) >= 3) {
    ord <- unlist(lapply(levels(g), function(k) {
      idx <- which(g == k)
      if (length(idx) < 3) return(idx)
      Z <- M[, idx, drop = FALSE]; Z[!is.finite(Z)] <- 0
      hc <- tryCatch(stats::hclust(stats::dist(t(Z)), method = "ward.D2"), error = function(e) NULL)
      if (is.null(hc)) idx else idx[hc$order]
    }), use.names = FALSE)
  }
  list(order = ord, group = g[ord], by = by,
       label = sprintf("columns grouped by %s%s", by,
                       if (isTRUE(cluster_within)) ", clustered within each group" else ""))
}

byod_sample_heatmap_gg <- function(M, meta = NULL, by = "cluster", title = "",
                                   value_name = "score", subtitle = NULL,
                                   max_rows = 60L, cluster_rows = TRUE,
                                   show_col_labels = NULL, stops = NULL) {
  if (is.null(M) || !length(M) || !nrow(M) || !ncol(M))
    return(emptyplot_gg("Nothing to show."))
  if (nrow(M) > max_rows) {
    keep <- head(order(byod_rowvars(M), decreasing = TRUE), max_rows)
    M <- M[sort(keep), , drop = FALSE]
  }
  ord <- byod_heat_order(M, meta, by)
  M <- M[, ord$order, drop = FALSE]
  if (isTRUE(cluster_rows) && nrow(M) >= 3) {
    Z <- M; Z[!is.finite(Z)] <- 0
    hr <- tryCatch(stats::hclust(stats::dist(Z), method = "ward.D2"), error = function(e) NULL)
    if (!is.null(hr)) M <- M[hr$order, , drop = FALSE]
  }
  cn <- colnames(M)
  df <- data.frame(
    col = rep(factor(cn, levels = cn), each = nrow(M)),
    row = rep(factor(rownames(M), levels = rev(rownames(M))), times = ncol(M)),
    value = as.vector(M), stringsAsFactors = FALSE)
  lim <- max(abs(df$value), na.rm = TRUE)

  show_lab <- show_col_labels %||% (ncol(M) <= 24L)
  g <- ggplot(df, aes(col, row, fill = value)) +
    geom_raster() +

    (if (length(stops))
       scale_fill_gradientn(colours = stops, limits = c(-lim, lim), name = value_name)
     else
       scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                            midpoint = 0, limits = c(-lim, lim), name = value_name)) +
    scale_x_discrete(expand = expansion(mult = c(0, 0))) +
    scale_y_discrete(expand = expansion(mult = c(0, 0))) +
    labs(title = paste(strwrap(title, width = 62), collapse = "\n"),
         subtitle = paste(strwrap(paste(c(subtitle, ord$label,
                                          sprintf("%d rows x %d columns.", nrow(M), ncol(M))),
                                        collapse = " "), width = 96), collapse = "\n"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 7.5, colour = "grey30", lineheight = 1.15),
          axis.text.y   = element_text(size = 6.5),
          axis.text.x   = if (show_lab) element_text(angle = 90, hjust = 1, vjust = .5, size = 5.5)
                          else element_blank(),
          axis.ticks    = element_blank(),
          panel.grid    = element_blank(),
          legend.key.height = unit(.55, "cm"))

  if (!is.null(ord$group)) {
    gp <- ord$group
    runs <- rle(as.character(gp))
    ends <- cumsum(runs$lengths); starts <- ends - runs$lengths + 1
    mids <- (starts + ends) / 2
    seps <- head(ends, -1) + .5
    g <- g +
      geom_vline(xintercept = seps, colour = "grey25", linewidth = .35) +
      annotate("text", x = mids, y = nrow(M) + .8, label = runs$values,
               size = 2.7, vjust = 0, colour = "grey15", fontface = "bold") +
      coord_cartesian(ylim = c(.5, nrow(M) + 2), clip = "off") +
      theme(plot.margin = margin(16, 10, 6, 10))
  }
  g
}

byod_dist_gg <- function(d, cols = NULL) {
  q <- apply(d$expr, 2, stats::quantile, probs = c(.25, .5, .75), na.rm = TRUE)
  samp <- colnames(d$expr)
  tp <- factor(d$meta[samp, "Timepoint"], levels = d$tp_levels)
  df <- data.frame(sample = samp, tp = tp, lo = q[1, ], mid = q[2, ], hi = q[3, ],
                   stringsAsFactors = FALSE)
  df <- df[order(df$tp, df$sample), ]; df$idx <- seq_len(nrow(df))
  df$tp <- droplevels(df$tp); lv <- levels(df$tp)
  cc <- cols %||% setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = "Timepoint")), lv)
  ggplot(df, aes(idx, mid, colour = tp)) +
    geom_linerange(aes(ymin = lo, ymax = hi), linewidth = .35, alpha = .8) +
    geom_point(size = 1.1) +
    scale_colour_manual(values = cc, name = NULL) +
    labs(x = "sample (ordered by timepoint)", y = paste0(d$unit, " (median, IQR)"),
         title = sprintf("%s - %s profiles", d$omic %||% "layer", format(ncol(d$expr), big.mark = ","))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(size = 12),
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

byod_gti_overview <- function(o, d = NULL, fdr = 0.05, n_label = 8L) {
  if (!isTRUE(o$ok)) return(list(ok = FALSE, msg = o$msg %||% "not available"))
  tb <- o$tab
  s <- data.frame(Feature = tb$Feature, P.Value = tb$P.Value,
                  adj.P.Val = tb$adj.P.Val, stringsAsFactors = FALSE)
  s$y <- -log10(pmax(s$P.Value, .Machine$double.xmin))

  eff <- rep(NA_real_, nrow(s)); names(eff) <- s$Feature
  gv <- o$group_var
  if (!is.null(d) && !is.null(gv) && gv %in% names(d$meta)) {
    tf <- droplevels(factor(d$meta$Timepoint, levels = d$tp_levels))
    gf <- droplevels(factor(as.character(d$meta[[gv]])))
    ok <- !is.na(tf) & !is.na(gf)
    lv2 <- levels(gf)
    if (length(lv2) == 2L && nlevels(droplevels(tf[ok])) >= 2L) {
      feats <- intersect(s$Feature, rownames(d$expr))
      E <- d$expr[feats, ok, drop = FALSE]
      key <- interaction(gf[ok], tf[ok], drop = TRUE, sep = "\r")
      M <- vapply(split(seq_len(ncol(E)), key), function(i)
             rowMeans(E[, i, drop = FALSE], na.rm = TRUE), numeric(nrow(E)))
      parts <- do.call(rbind, strsplit(colnames(M), "\r", fixed = TRUE))
      tps <- levels(droplevels(tf[ok]))
      gA <- M[, parts[, 1] == lv2[1], drop = FALSE]
      gB <- M[, parts[, 1] == lv2[2], drop = FALSE]
      colnames(gA) <- parts[parts[, 1] == lv2[1], 2]
      colnames(gB) <- parts[parts[, 1] == lv2[2], 2]
      shared <- intersect(intersect(colnames(gA), colnames(gB)), tps)
      if (length(shared) >= 2L) {
        gA <- gA[, shared, drop = FALSE]; gB <- gB[, shared, drop = FALSE]
        dA <- gA - gA[, 1]; dB <- gB - gB[, 1]
        eff[feats] <- apply(abs(dA - dB), 1, max, na.rm = TRUE)
      }
    }
  }
  s$effect <- unname(eff[s$Feature])
  have_eff <- any(is.finite(s$effect))
  if (!have_eff) {                       # no layer supplied: fall back to the statistic
    stat <- if ("F" %in% names(tb)) "F" else "t"
    s$effect <- tb[[stat]]
  }
  s$sig <- ifelse(is.finite(s$adj.P.Val) & s$adj.P.Val < fdr,
                  sprintf("FDR < %g", fdr), "n.s.")
  s <- s[order(s$adj.P.Val), ]
  lab <- head(s$Feature[s$sig != "n.s."], n_label)
  list(ok = TRUE, summary = s, labelled = lab, fdr = fdr, has_effect = have_eff,
       xlab = if (have_eff) sprintf("largest divergence in change from baseline (%s)", d$unit)
              else "test statistic",
       ylab = "-log10 P",
       n_sig = sum(s$sig != "n.s.", na.rm = TRUE),
       ttl = sprintf("%s x time - %s, n = %s %s",
                     o$group_var %||% "group",
                     if (isTRUE(o$paired)) "within-unit paired" else "subject-blocked",
                     o$n_units %||% NA, o$unit %||% "units"),
       sub = sprintf("%s of %s features at FDR < %g",
                     format(sum(s$sig != "n.s.", na.rm = TRUE), big.mark = ","),
                     format(nrow(s), big.mark = ","), fdr))
}
byod_gti_overview_gg <- function(ov, cols = NULL) {
  if (!isTRUE(ov$ok)) return(emptyplot_gg(ov$msg %||% "not available"))
  s <- ov$summary
  cc <- c("n.s." = "grey72"); cc[sprintf("FDR < %g", ov$fdr)] <- unname(VOLCANO_DEFAULT[["Up"]])
  if (!is.null(cols)) { ok <- intersect(names(cols), names(cc)); cc[ok] <- cols[ok] }
  g <- ggplot(s, aes(effect, y, colour = sig)) +
    geom_point(size = 1.1, alpha = .7) +
    scale_colour_manual(values = cc, name = NULL) +
    labs(x = ov$xlab, y = ov$ylab, title = gsub("·", "-", ov$ttl),
         subtitle = ov$sub) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 10, colour = "grey30"),
          panel.grid.minor = element_blank())
  if (length(ov$labelled) && requireNamespace("ggrepel", quietly = TRUE))
    g <- g + ggrepel::geom_text_repel(data = s[s$Feature %in% ov$labelled, , drop = FALSE],
      aes(label = Feature), seed = trace_seed(), size = 3, colour = "grey15", max.overlaps = Inf,
      min.segment.length = 0, box.padding = .3, segment.size = .3, show.legend = FALSE)
  g
}

byod_conc_hist_gg <- function(obj) {
  if (is.null(obj$tab) || !nrow(obj$tab)) return(emptyplot_gg(obj$msg %||% "Need matched samples"))
  tab <- obj$tab; tab$class <- factor(tab$class, levels = names(BYOD_CONC_COLS))
  ggplot(tab, aes(rho, fill = class)) +
    geom_histogram(bins = 40, colour = NA) +
    scale_fill_manual(values = BYOD_CONC_COLS, name = NULL) +
    labs(x = paste0("per-gene mRNA-protein ", obj$method, " r"), y = "genes",
         title = sprintf("%d shared genes across %d matched samples - median r = %.2f, %.0f%% discordant",
                         obj$n, obj$nS, obj$med_rho, 100 * obj$frac_disc)) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(size = 11),
          panel.grid.minor = element_blank())
}

byod_pathway_concordance <- function(pa1, pa2, lab1 = "mRNA", lab2 = "Protein",
                                     method = "spearman") {
  if (is.null(pa1$mat) || is.null(pa2$mat) || !nrow(pa1$mat) || !nrow(pa2$mat))
    return(list(ok = FALSE, msg = "Compute pathway activity on both layers first."))
  P <- intersect(rownames(pa1$mat), rownames(pa2$mat))
  S <- intersect(colnames(pa1$mat), colnames(pa2$mat))
  if (length(P) < 3 || length(S) < 4)
    return(list(ok = FALSE, msg = sprintf(
      "Needs pathways and samples measured in both layers: %d shared pathways, %d shared samples.",
      length(P), length(S))))
  r <- vapply(P, function(p) suppressWarnings(stats::cor(
    as.numeric(pa1$mat[p, S]), as.numeric(pa2$mat[p, S]), method = method,
    use = "complete.obs")), numeric(1))
  tab <- data.frame(pathway = P, r = unname(r), stringsAsFactors = FALSE)
  tab <- tab[is.finite(tab$r), ]; tab <- tab[order(-tab$r), ]
  list(ok = TRUE, tab = tab, n_pathways = nrow(tab), n_samples = length(S),
       method = method, lab1 = lab1, lab2 = lab2, median_r = stats::median(tab$r),
       ttl = sprintf("%d pathways across %d matched samples - median %s r = %.2f",
                     nrow(tab), length(S), method, stats::median(tab$r)))
}
byod_pathway_concordance_gg <- function(obj, top_n = 25L) {
  if (!isTRUE(obj$ok)) return(emptyplot_gg(obj$msg %||% "not available"))
  tb <- obj$tab
  tb <- rbind(utils::head(tb, ceiling(top_n / 2)), utils::tail(tb, floor(top_n / 2)))
  tb <- tb[!duplicated(tb$pathway), ]
  tb$pathway <- factor(tb$pathway, levels = rev(tb$pathway))
  ggplot(tb, aes(r, pathway, fill = r > 0)) +
    geom_col(width = .7) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = .35) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), guide = "none") +
    labs(x = sprintf("%s r between %s and %s", obj$method, obj$lab1, obj$lab2),
         y = NULL, title = gsub("·", "-", obj$ttl)) +
    byod_theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 11), panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 8))
}

byod_diablo_sample_gg <- function(obj, comps = c(1L, 2L)) {
  m <- obj$model
  if (is.null(m)) return(emptyplot_gg("run the multiblock model first"))
  V <- tryCatch(m$variates, error = function(e) NULL)
  if (is.null(V) || !length(V)) return(emptyplot_gg("the model carries no sample scores"))

  blocks <- names(V)[vapply(V, function(z) is.matrix(z) && ncol(z) >= max(comps), logical(1))]
  if (!length(blocks)) return(emptyplot_gg("too few components for a projection"))
  X <- Reduce(`+`, lapply(V[blocks], function(z) z[, comps, drop = FALSE])) / length(blocks)
  df <- data.frame(x = X[, 1], y = X[, 2],
                   group = as.character(obj$Y %||% m$Y), stringsAsFactors = FALSE)
  df <- df[is.finite(df$x) & is.finite(df$y), , drop = FALSE]
  lv <- if (is.factor(obj$Y)) levels(droplevels(obj$Y)) else unique(df$group)
  df$group <- factor(df$group, levels = lv)
  cc <- setNames(as.character(resolve_appearance(lv, PAL_NAMES[1], splitby = obj$outcome)), lv)
  ggplot(df, aes(x, y, colour = group)) +
    geom_point(size = 2.2, alpha = .85) +
    stat_ellipse(level = .95, linewidth = .4, show.legend = FALSE) +
    scale_colour_manual(values = cc, name = NULL) +
    labs(x = sprintf("consensus component %d", comps[1]),
         y = sprintf("consensus component %d", comps[2]),
         title = sprintf("Multiblock projection - descriptive fit on all %d matched observations, NOT cross-validated",
                         nrow(df))) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(size = 10.5),
          panel.grid.minor = element_blank())
}

byod_diablo_loadings <- function(obj, comp = 1L, top_n = 20L) {
  m <- obj$model
  if (is.null(m) || is.null(m$loadings)) return(list(ok = FALSE, msg = "run the multiblock model first"))
  blk <- setdiff(names(m$loadings), "Y")
  rows <- do.call(rbind, lapply(blk, function(b) {
    L <- m$loadings[[b]]
    if (is.null(L) || ncol(L) < comp) return(NULL)
    v <- L[, comp]; v <- v[v != 0]
    if (!length(v)) return(NULL)
    v <- v[order(-abs(v))][seq_len(min(top_n, length(v)))]
    data.frame(block = b, feature = names(v), loading = unname(v), stringsAsFactors = FALSE)
  }))
  if (is.null(rows) || !nrow(rows))
    return(list(ok = FALSE, msg = sprintf("no feature has a non-zero loading on component %d", comp)))
  list(ok = TRUE, tab = rows, comp = comp,
       ttl = sprintf("Multiblock loadings, component %d - %d selected features across %d blocks",
                     comp, nrow(rows), length(unique(rows$block))))
}
byod_diablo_loadings_gg <- function(o) {
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg %||% "not available"))
  tb <- o$tab
  tb$feature <- factor(tb$feature, levels = rev(unique(tb$feature[order(tb$block, -abs(tb$loading))])))
  ggplot(tb, aes(loading, feature, fill = loading > 0)) +
    geom_col(width = .72) +
    geom_vline(xintercept = 0, colour = "grey40", linewidth = .35) +
    facet_wrap(~ block, scales = "free_y") +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), guide = "none") +
    labs(x = sprintf("loading on component %d", o$comp), y = NULL,
         title = gsub("·", "-", o$ttl)) +
    byod_theme_bw(base_size = 12) +
    theme(plot.title = element_text(size = 11), panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 7.5))
}

byod_cv_confusion_gg <- function(cv) {
  if (!isTRUE(cv$ok) || is.null(cv$confusion)) return(emptyplot_gg(cv$msg %||% "run the cross-validation first"))
  M <- as.matrix(cv$confusion)
  df <- expand.grid(observed = rownames(M), predicted = colnames(M),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df$n <- as.vector(M)
  df$recall <- df$n / rowSums(M)[df$observed]
  df$observed  <- factor(df$observed,  levels = rev(rownames(M)))
  df$predicted <- factor(df$predicted, levels = colnames(M))
  ggplot(df, aes(predicted, observed, fill = recall)) +
    geom_tile(colour = "white", linewidth = .6) +
    geom_text(aes(label = ifelse(n > 0, format(n, big.mark = ","), ""), colour = recall > .55), size = 3) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey12"), guide = "none") +
    scale_fill_gradient(low = "#F2F6F9", high = "#1E3A5F", limits = c(0, 1),
                        labels = scales::percent_format(accuracy = 1), name = "recall") +
    labs(x = "predicted", y = "observed",
         title = sprintf("%d-fold x %d repeats, whole %s held out - balanced accuracy %.1f%% against %.1f%% for chance",
                         cv$folds, cv$repeats, cv$unit %||% "unit",
                         100 * cv$balanced_accuracy, 100 * cv$chance)) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(size = 10.5), panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1))
}

cs_traj_scatter_gg <- function(o, cols = NULL, max_points = 20000L, seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg %||% "not available"))
  d <- o$long
  if (nrow(d) > max_points) d <- d[sample.int(nrow(d), max_points), , drop = FALSE]
  cc <- cols %||% c("same direction" = PIG, "opposite" = GOLD)
  lim <- stats::quantile(abs(c(o$long$A, o$long$B)), .995, na.rm = TRUE) * 1.05
  ggplot(d, aes(A, B, colour = agreement)) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = .3) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = .3) +
    geom_point(size = .7, alpha = .45) +
    scale_colour_manual(values = cc, name = NULL) +
    coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    facet_wrap(~ timepoint) +
    labs(x = o$labA, y = o$labB, title = gsub("·", "-", o$ttl)) +
    byod_theme_bw(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(size = 10.5),
          panel.grid.minor = element_blank())
}
cs_traj_hist_gg <- function(o) {
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg %||% "not available"))
  r <- o$per_gene$trajectory_r; r <- r[is.finite(r)]
  if (!length(r)) return(emptyplot_gg(
    "Needs at least three post-baseline timepoints for a per-feature trajectory correlation"))
  med <- stats::median(r)
  ggplot(data.frame(r = r), aes(r)) +
    geom_histogram(bins = 60, fill = PIG, colour = NA, alpha = .85) +
    geom_vline(xintercept = med, colour = GOLD, linewidth = .8, linetype = "dashed") +
    annotate("text", x = med, y = Inf, vjust = 1.6, hjust = -0.08,
             label = sprintf("median %.3f", med), size = 3.2, colour = "grey20") +
    coord_cartesian(xlim = c(-1, 1)) +
    labs(x = sprintf("per-feature trajectory correlation (%s vs %s)", o$labA, o$labB),
         y = "features",
         title = sprintf("%s features with a complete trajectory in both studies",
                         format(length(r), big.mark = ","))) +
    byod_theme_bw(base_size = 13) +
    theme(plot.title = element_text(size = 11), panel.grid.minor = element_blank())
}
cs_traj_gene_gg <- function(o, gene, cols = NULL) {
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg %||% "not available"))
  if (!gene %in% rownames(o$P)) return(emptyplot_gg(paste(gene, "has no complete trajectory in both studies")))
  eps <- o$endpoints
  df <- rbind(
    data.frame(timepoint = factor(eps, levels = eps), value = as.numeric(o$P[gene, ]), study = o$labA),
    data.frame(timepoint = factor(eps, levels = eps), value = as.numeric(o$H[gene, ]), study = o$labB))
  df$study <- factor(df$study, levels = c(o$labA, o$labB))
  cc <- cols %||% setNames(c(PIG, HUMAN), c(o$labA, o$labB))
  pg <- o$per_gene[o$per_gene$feature == gene, ]
  ggplot(df, aes(timepoint, value, colour = study, group = study)) +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = .3) +
    geom_line(linewidth = 1) + geom_point(size = 2.6) +
    scale_colour_manual(values = cc, name = NULL) +
    labs(x = NULL, y = CS_MEASURE_LAB[[o$measure]] %||% "effect",
         title = if (nrow(pg)) sprintf("%s - trajectory r = %.2f", gene, pg$trajectory_r[1]) else gene) +
    byod_theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(size = 11),
          panel.grid.minor = element_blank())
}

byod_gsea_curve <- function(stats_vec, collection, pathway, d = NULL) {
  if (length(pathway) != 1L || is.na(pathway) || !nzchar(pathway))
    return(list(ok = FALSE, msg = "Select one pathway for the enrichment curve."))
  sets <- annotation_sets(d, collection)
  if (!is.null(d)) stats_vec <- stats_vec[names(stats_vec) %in% annotation_universe(d, collection)]
  stats_vec <- stats_vec[is.finite(stats_vec)]
  if (is.null(sets) || !length(sets)) return(list(ok = FALSE, msg = "Unknown collection."))
  key <- names(sets)[match(gsub("_", " ", toupper(pathway)),
                           gsub("_", " ", toupper(names(sets))))]
  if (is.na(key)) return(list(ok = FALSE, msg = paste(pathway, "is not in", collection)))
  g <- intersect(sets[[key]], names(stats_vec))
  if (length(g) < 3) return(list(ok = FALSE, msg = sprintf(
    "%s has only %d of its genes in this matrix.", key, length(g))))
  s <- sort(stats_vec, decreasing = TRUE)
  hit <- which(names(s) %in% g)
  n <- length(s); nh <- length(hit)

  w <- abs(s[hit]); Nr <- sum(w)
  if (n == nh || Nr == 0) return(list(ok = FALSE, msg = "An enrichment curve needs non-hit features and nonzero hit statistics."))
  inc <- numeric(n); inc[hit] <- w / Nr
  dec <- rep(1 / (n - nh), n); dec[hit] <- 0
  run <- cumsum(inc - dec)
  list(ok = TRUE, pathway = key, collection = collection,
       curve = data.frame(rank = seq_len(n), running = run),
       hits = data.frame(rank = hit),
       es = run[which.max(abs(run))], n_genes = nh, n_ranked = n,
       ttl = sprintf("%s - %d of %d ranked features, ES = %.3f",
                     key, nh, n, run[which.max(abs(run))]))
}
byod_gsea_curve_gg <- function(o) {
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg %||% "not available"))
  cur <- o$curve
  h <- max(abs(cur$running)) * 0.12
  ggplot() +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = .3) +
    geom_line(data = cur, aes(rank, running), colour = PIG, linewidth = .8) +
    geom_segment(data = o$hits, aes(x = rank, xend = rank, y = -h, yend = 0),
                 colour = "grey25", linewidth = .25) +
    labs(x = "rank in the ordered feature list", y = "running enrichment score",
         title = gsub("·", "-", o$ttl)) +
    byod_theme_bw(base_size = 13) +
    theme(plot.title = element_text(size = 11), panel.grid.minor = element_blank())
}

source("compare.R", local = FALSE)
source("bundle.R", local = FALSE)

# Folder numbers are stable identifiers shared by session exports.
ANALYSIS_REGISTRY <- local({
  r <- function(n, module, folder) data.frame(number = n, module = module, folder = folder,
                                              stringsAsFactors = FALSE)
  do.call(rbind, list(
    r( 1L, "Normalization QC",             "Normalisation_QC"),
    r( 2L, "Landscape",                    "Landscape"),
    r( 3L, "DE (limma)",                   "Differential"),
    r( 4L, "Group x time interaction",     "Group_by_time"),
    r( 5L, "Trajectory",                   "Trajectory"),
    r( 6L, "Temporal programs",            "Temporal_programs"),
    r( 7L, "Association",                  "Association"),
    r( 8L, "Subject variables: features",  "Subject_variables_features"),
    r( 9L, "Subject variables: pathways",  "Subject_variables_pathways"),
    r(10L, "Subject variables: screen",    "Subject_variables_screen"),
    r(11L, "Clinical correlates",          "Clinical_correlates"),
    r(12L, "Enrichment",                   "Enrichment"),
    r(13L, "Pathways",                     "Pathway_activity"),
    r(14L, "Module score",                 "Module_score"),
    r(15L, "Heatmap",                      "Heatmap"),
    r(16L, "WGCNA",                        "Co_expression"),
    r(17L, "Multi-timepoint",              "Multi_timepoint"),
    r(18L, "Response composition",         "Response_composition"),
    r(19L, "Trajectory consistency",       "Trajectory_consistency"),
    r(20L, "Baseline sensitivity",         "Baseline_sensitivity"),
    r(21L, "Integration",                  "Multiomics_integration"),
    r(22L, "Concordance",                  "Multiomics_concordance"),
    r(23L, "DIABLO",                       "DIABLO"),
    r(24L, "DIABLO cross-validation",      "DIABLO_cross_validation"),
    r(25L, "Cross-study comparison",       "Cross_study_trajectory"),
    r(26L, "cross-dataset differential concordance", "Cross_dataset_differential"),
    r(27L, "cross-dataset pathway concordance",      "Cross_dataset_pathway")))
})

analysis_folder <- function(module) {
  i <- match(module, ANALYSIS_REGISTRY$module)
  if (is.na(i)) return(NA_character_)
  sprintf("%02d_%s", ANALYSIS_REGISTRY$number[i], ANALYSIS_REGISTRY$folder[i])
}

EXPORT_PANEL_MAP <- rbind(
  data.frame(panel = "byod_qc_appr",        module = "Normalization QC",            figure = "per_sample_distribution"),
  data.frame(panel = "byod_dr_appr",        module = "Landscape",                   figure = "projection"),
  data.frame(panel = "byod_dr_appr_heat",   module = "Landscape",                   figure = "projection_continuous"),
  data.frame(panel = "byod_de_appr",        module = "DE (limma)",                  figure = "volcano"),
  data.frame(panel = "byod_de_hm_appr",     module = "DE (limma)",                  figure = "top_features_heatmap"),
  data.frame(panel = "byod_de_ora_appr",    module = "DE (limma)",                  figure = "top_features_overrepresentation"),
  data.frame(panel = "byod_gti_ov_appr",    module = "Group x time interaction",    figure = "interaction_overview"),
  data.frame(panel = "byod_gti_tc_appr",    module = "Group x time interaction",    figure = "feature_trajectory"),
  data.frame(panel = "byod_gti_hm_appr",    module = "Group x time interaction",    figure = "top_features_heatmap"),
  data.frame(panel = "byod_gti_en_appr",    module = "Group x time interaction",    figure = "top_features_overrepresentation"),
  data.frame(panel = "byod_ge_appr",        module = "Trajectory",                  figure = "feature_trajectory"),
  data.frame(panel = "byod_tp_prof_appr",   module = "Temporal programs",           figure = "cluster_profiles"),
  data.frame(panel = "byod_tp_heat_appr",   module = "Temporal programs",           figure = "cluster_heatmap"),
  data.frame(panel = "byod_tp_en_appr",     module = "Temporal programs",           figure = "cluster_enrichment"),
  data.frame(panel = "byod_as_appr",        module = "Association",                 figure = "association"),
  data.frame(panel = "byod_cc_appr",        module = "Subject variables: features", figure = "volcano"),
  data.frame(panel = "byod_cc_sc_appr",     module = "Subject variables: features", figure = "single_feature"),
  data.frame(panel = "byod_cc_en_appr",     module = "Subject variables: features", figure = "enrichment"),
  data.frame(panel = "byod_pxc_appr",       module = "Subject variables: pathways", figure = "pathway_by_characteristic"),
  data.frame(panel = "byod_vs_appr",        module = "Subject variables: screen",   figure = "screen"),
  data.frame(panel = "byod_en_appr",        module = "Enrichment",                  figure = "enrichment"),
  data.frame(panel = "byod_en_curve_appr",  module = "Enrichment",                  figure = "running_enrichment"),
  data.frame(panel = "byod_pw_hm_appr",     module = "Pathways",                    figure = "activity_heatmap"),
  data.frame(panel = "byod_pw_tc_appr",     module = "Pathways",                    figure = "pathway_trajectory"),
  data.frame(panel = "byod_ms_appr",        module = "Module score",                figure = "module_score"),
  data.frame(panel = "byod_hm_appr",        module = "Heatmap",                     figure = "heatmap"),
  data.frame(panel = "byod_wg_appr",        module = "WGCNA",                       figure = "module_trait_heatmap"),
  data.frame(panel = "byod_mc_volc_appr",   module = "Multi-timepoint",             figure = "volcano_per_timepoint"),
  data.frame(panel = "byod_mc_upset_appr",  module = "Multi-timepoint",             figure = "shared_and_specific"),
  data.frame(panel = "byod_rc_appr",        module = "Response composition",        figure = "composition"),
  data.frame(panel = "byod_tcon_appr",      module = "Trajectory consistency",      figure = "consistency"),
  data.frame(panel = "byod_bs_appr",        module = "Baseline sensitivity",        figure = "effect_per_baseline"),
  data.frame(panel = "byod_bs_ov_appr",     module = "Baseline sensitivity",        figure = "overlap_between_baselines"),
  data.frame(panel = "byod_int_appr",       module = "Integration",                 figure = "layer_vs_layer"),
  data.frame(panel = "byod_conc_hist_appr", module = "Concordance",                 figure = "correlation_distribution"),
  data.frame(panel = "byod_pconc_appr",     module = "Concordance",                 figure = "pathway_level"),
  data.frame(panel = "byod_diablo_appr",    module = "DIABLO",                      figure = "sample_projection"),
  data.frame(panel = "byod_diablo_load_appr", module = "DIABLO",                    figure = "selected_features"),
  data.frame(panel = "byod_diablo_cv_appr", module = "DIABLO cross-validation",     figure = "confusion_matrix"),
  data.frame(panel = "cs_appr",             module = "Cross-study comparison",      figure = "one_feature_builtin"),
  data.frame(panel = "csr_appr",            module = "Cross-study comparison",      figure = "trajectory_scatter"),
  data.frame(panel = "csr_hist_appr",       module = "Cross-study comparison",      figure = "trajectory_agreement"),
  data.frame(panel = "csr_gene_appr",       module = "Cross-study comparison",      figure = "one_feature"),
  data.frame(panel = "xs_design_appr",      module = "cross-dataset differential concordance", figure = "design_comparison"),
  data.frame(panel = "xs_de_sc_appr",       module = "cross-dataset differential concordance", figure = "fold_change_both_studies"),
  data.frame(panel = "xs_de_quad_appr",     module = "cross-dataset differential concordance", figure = "agreement_of_significant"),
  data.frame(panel = "xs_de_null_appr",     module = "cross-dataset differential concordance", figure = "permutation_null"),
  data.frame(panel = "xs_pw_sc_appr",       module = "cross-dataset pathway concordance",      figure = "gene_set_scores"),
  data.frame(panel = "xs_pw_bar_appr",      module = "cross-dataset pathway concordance",      figure = "top_gene_sets"),
  stringsAsFactors = FALSE)

# Count analysis panels, excluding one occurrence of each top-level navigation page.
N_MODULES <- tryCatch({
  u <- paste(readLines("ui.R", warn = FALSE), collapse = "\n")

  titles <- regmatches(u, gregexpr('nav_panel\\(\\s*\n?\\s*"[^"]+"', u))[[1]]
  titles <- gsub('^nav_panel\\(\\s*\n?\\s*"|"$', "", titles)

  pages <- c("Overview", "Kidney Perfusion Reference", "Your data", "Analyse",
             "Cross-study comparison", "Reproducibility", "About")

  n <- length(titles)
  for (pg in pages) if (any(titles == pg)) n <- n - 1L
  n
}, error = function(e) NA_integer_)
