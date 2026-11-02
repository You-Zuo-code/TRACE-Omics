MIN_UNITS$cross_dataset <- 3L    # per side; each side is also gated by its own module

XS_MIN_SHARED <- 50L

XS_PERM_OFF <- paste("Not computed. Tick \u201cTest the agreement against a permutation null\u201d",
                     "in the sidebar and run the comparison again.")

xs_decl <- function(d, what = c("subject", "nest")) {
  v <- d[[match.arg(what)]]
  if (is.null(v) || !length(v) || !is.character(v) || !nzchar(v[1])) NULL else as.character(v[1])
}

xs_align_timepoints <- function(dA, dB, mode = c("label", "position", "manual"),
                                manual = NULL) {
  mode <- match.arg(mode)
  if (mode != "manual") {
    al <- cs_align_timepoints(dA, dB, mode)
    al$mode <- mode
    return(al)
  }
  a <- cs_timepoints(dA); b <- cs_timepoints(dB)
  if (is.null(manual) || !length(manual))
    return(list(align = "manual", mode = "manual", A = character(0), B = character(0),
                label = character(0), n = 0,
                note = "No timepoint correspondence has been declared yet."))
  manual <- manual[names(manual) %in% a]
  manual <- manual[!is.na(manual) & nzchar(manual) & manual %in% b]

  if (anyDuplicated(unname(manual)))
    return(list(align = "manual", mode = "manual", A = character(0), B = character(0),
                label = character(0), n = 0,
                note = sprintf("More than one timepoint of the first study is mapped onto %s of the second. Each timepoint may be used once.",
                               paste(unique(unname(manual)[duplicated(unname(manual))]), collapse = ", "))))

  keys <- a[a %in% names(manual)]
  list(align = "manual", mode = "manual", A = keys, B = unname(manual[keys]),
       label = sprintf("%s / %s", keys, unname(manual[keys])), n = length(keys),
       note = sprintf("Declared correspondence: %s. The engine takes this as given; it cannot check that the elapsed times are comparable.",
                      paste(sprintf("%s→%s", keys, unname(manual[keys])), collapse = ", ")))
}

xs_feature_bridge <- function(dA, dB, one_to_one_only = TRUE) {
  fp <- cs_feature_pairs(dA, dB, one_to_one_only)
  nA <- nrow(dA$expr); nB <- nrow(dB$expr); k <- length(fp$features)
  fp$n_A <- nA; fp$n_B <- nB; fp$n_shared <- k
  fp$cover_A <- if (nA) k / nA else NA_real_
  fp$cover_B <- if (nB) k / nB else NA_real_
  fp$enough <- k >= XS_MIN_SHARED
  fp$note <- if (!k)
      paste("No eligible feature pairs.", fp$rule, "Check the layer annotation reports. Pathway comparison is available when both layers have compatible gene-set annotations.")
    else sprintf("%s shared features (%.0f%% of %s, %.0f%% of %s), matched on %s.",
                 format(k, big.mark = ","), 100 * fp$cover_A, format(nA, big.mark = ","),
                 100 * fp$cover_B, format(nB, big.mark = ","), fp$rule)
  fp
}

xs_rows <- function(fp, side = c("A", "B")) {
  side <- match.arg(side)
  k <- if (side == "A") fp$keyA else fp$keyB
  if (is.null(k)) fp$features else unname(k[fp$features])
}

xs_spec <- function(groupvar, test, reference, at_tp = "All", paired = TRUE)
  list(groupvar = groupvar, test = test, reference = reference,
       at_tp = at_tp %||% "All", paired = isTRUE(paired))

xs_spec_text <- function(spec)
  sprintf("%s: %s vs %s%s", spec$groupvar, spec$test, spec$reference,
          if (!identical(spec$at_tp, "All")) sprintf(" (at %s)", spec$at_tp) else "")

# Translate test-minus-reference to the single-study engine's argument order.
xs_run_contrast <- function(d, spec, label = "study") {

  tt <- tryCatch(byod_run_limma(d, spec$groupvar, spec$reference, spec$test,
                                subject_col = xs_decl(d, "subject"), paired = spec$paired,
                                at_tp = spec$at_tp),
                 error = function(e) { z <- data.frame(); attr(z, "declined") <- conditionMessage(e); z })
  oc <- byod_outcome_of(tt)
  list(label = label, spec = spec, tt = tt,
       declared_subject = xs_decl(d, "subject"), declared_nest = xs_decl(d, "nest"),
       ok = identical(oc$status, "completed"),
       declined = if (identical(oc$status, "completed")) NULL else oc$reason,
       text = xs_spec_text(spec),
       n = attr(tt, "n") %||% NA_integer_,
       n_units = attr(tt, "n_units") %||% NA_integer_,
       unit = attr(tt, "unit") %||% NA_character_,
       paired = isTRUE(attr(tt, "paired")),
       excluded_profiles = attr(tt, "excluded_profiles") %||% 0L,
       excluded_subjects = attr(tt, "excluded_subjects") %||% character(),
       nested = isTRUE(attr(tt, "nested")))
}

xs_exclusion_text <- function(o) {
  paste(vapply(list(o$sideA, o$sideB), function(s) {
    if (is.null(s) || (s$excluded_profiles %||% 0L) == 0L) return("")
    sprintf("%s: excluded %d profiles from %d incomplete subjects under the declared pairing. ",
      s$label, s$excluded_profiles, length(s$excluded_subjects))
  }, character(1)), collapse = "")
}

xs_design_table <- function(sA, sB, dA, dB, labA = "A", labB = "B",
                            bridge = NULL, align = NULL) {
  fmt <- function(x, max_n = 6L) {
    if (is.null(x) || !length(x)) return("—")
    if (length(x) > max_n)
      return(sprintf("%s, and %d more", paste(utils::head(x, max_n), collapse = ", "),
                     length(x) - max_n))
    paste(x, collapse = ", ")
  }

  idlab <- function(d) d$ids$label %||%
    tryCatch(byod_detect_ids(rownames(d$expr), d$omics %||% "other")$label,
             error = function(e) "unknown")
  n_units <- function(d) {
    u <- resolved_unit_col(d)
    if (is.null(u)) return(list(n = ncol(d$expr), lab = "profile (no subject declared)"))
    list(n = length(unique(as.character(d$meta[[u]]))), lab = u)
  }
  uA <- n_units(dA); uB <- n_units(dB)
  scA <- cs_scale_info(dA); scB <- cs_scale_info(dB)
  br <- bridge %||% xs_feature_bridge(dA, dB, FALSE)
  al <- align %||% xs_align_timepoints(dA, dB, "label")

  rows <- list()
  add <- function(what, a, b, verdict, why = "")
    rows[[length(rows) + 1L]] <<- data.frame(
      Property = what, A = a, B = b, Verdict = verdict, Note = why,
      stringsAsFactors = FALSE)

  add("Profiles", format(ncol(dA$expr), big.mark = ","), format(ncol(dB$expr), big.mark = ","),
      "informational", "Profiles are not the unit of analysis; the row below is.")
  add("Independent units", sprintf("%d (%s)", uA$n, uA$lab), sprintf("%d (%s)", uB$n, uB$lab),
      if (min(uA$n, uB$n) >= MIN_UNITS$cross_dataset) "ok" else "blocking",
      sprintf("Each side needs at least %d for any contrast. Every test below is at this resolution, never at profile resolution.",
              MIN_UNITS$cross_dataset))
  add("Subject declared", fmt(xs_decl(dA, "subject")), fmt(xs_decl(dB, "subject")),
      if (!is.null(xs_decl(dA, "subject")) && !is.null(xs_decl(dB, "subject"))) "ok" else "informational",
      "Without a subject the two sides cannot be paired within study, so their contrasts are between-profile and less powerful.")
  add("Higher-level unit", fmt(xs_decl(dA, "nest")), fmt(xs_decl(dB, "nest")), "informational",
      "Declared nesting is honoured inside each study; it does not have to match between them.")

  add("Timepoints", fmt(cs_timepoints(dA)), fmt(cs_timepoints(dB)),
      if (al$n >= 2) "ok" else "time-contrast blocked",
      sprintf("%d matched under the current rule (%s). This is needed only for a contrast BETWEEN timepoints; a contrast between groups inside each study does not require the two time grids to line up.",
              al$n, al$mode %||% al$align))
  add("Features", format(nrow(dA$expr), big.mark = ","), format(nrow(dB$expr), big.mark = ","),
      "informational", "")
  add("Identifier namespace", idlab(dA), idlab(dB),
      if (isTRUE(br$enough)) "ok" else "feature-level blocked",
      br$note)
  add("Shared features", format(br$n_shared, big.mark = ","), format(br$n_shared, big.mark = ","),
      if (isTRUE(br$enough)) "ok" else "feature-level blocked",
      sprintf("%.0f%% of A and %.0f%% of B. Pathway-level comparison does not need this.",
              100 * (br$cover_A %||% 0), 100 * (br$cover_B %||% 0)))
  add("Measurement unit", dA$unit %||% "unstated", dB$unit %||% "unstated",
      if (identical(scA$family, scB$family)) "ok" else "check",
      sprintf("Transform family: %s vs %s. Log fold-changes are comparable across scales; raw values are not, which is why nothing below compares raw values.",
              scA$family %||% "unresolved", scB$family %||% "unresolved"))
  add("Normalisation", dA$norm %||% "unstated", dB$norm %||% "unstated", "informational",
      "Each study was normalised on its own; no cross-study normalisation is applied anywhere.")
  add("Grouping variables", fmt(setdiff(group_vars(dA), "Timepoint")),
      fmt(setdiff(group_vars(dB), "Timepoint")), "informational",
      "The two studies need not share a variable name; you declare which contrast in A answers the same question as which contrast in B.")
  out <- do.call(rbind, rows)
  names(out)[2:3] <- c(labA, labB)
  attr(out, "blocking") <- sum(out$Verdict == "blocking")
  out
}

xs_design_verdict <- function(tab) {
  blk <- tab[tab$Verdict == "blocking", , drop = FALSE]
  tb  <- tab[tab$Verdict == "time-contrast blocked", , drop = FALSE]
  fl  <- tab[tab$Verdict == "feature-level blocked", , drop = FALSE]
  chk <- tab[tab$Verdict == "check", , drop = FALSE]
  if (nrow(blk))
    return(list(type = "danger", txt = sprintf(
      "These two studies cannot be compared as declared: %s. Fix the declaration, or pick a different pair.",
      paste(tolower(blk$Property), collapse = "; "))))
  if (nrow(tb))
    return(list(type = "warning", txt = paste(
      "The two time grids do not line up, so a contrast between timepoints cannot be matched yet.",
      "Declare the correspondence by hand above - which timepoint of one study is which timepoint of",
      "the other - or compare a contrast between groups instead, which does not need the grids to agree.")))
  if (nrow(fl))
    return(list(type = "warning", txt = paste(
      "Feature-level comparison is unavailable because the two layers do not share enough identifiers.",
      "Pathway-level comparison still works: it compares each study's enrichment result, which lives in a",
      "namespace both studies can reach.")))
  if (nrow(chk))
    return(list(type = "warning", txt = sprintf(
      "Comparable, with one thing to check: %s. Everything below compares fold-changes and enrichment scores, which survive a scale difference; nothing compares raw values.",
      paste(tolower(chk$Property), collapse = "; "))))
  list(type = "success", txt = paste(
    "Comparable. The same contrast will be fitted separately in each study, under each study's own declared",
    "design and unit floors, and only the two answers are compared."))
}

xs_pi1 <- function(p) {
  p <- p[is.finite(p)]
  if (length(p) < 20) return(NA_real_)
  max(0, min(1, 1 - 2 * mean(p > 0.5)))
}

xs_exchange_plan <- function(d, spec) {
  # Permute at the independent unit; paired contrasts use unit-level sign flips.
  gf <- grp_factor(d, spec$groupvar)
  sel <- which(as.character(gf) %in% c(spec$test, spec$reference) & !is.na(gf))
  if (!identical(spec$at_tp, "All"))
    sel <- sel[as.character(d$meta$Timepoint[sel]) == spec$at_tp]
  if (length(sel) < 4) return(NULL)
  ucol <- xs_decl(d, "nest") %||% xs_decl(d, "subject")
  if (is.null(ucol) || !ucol %in% names(d$meta)) {

    return(list(scheme = "profile", sel = sel, ucol = NULL,
                units = as.character(seq_along(sel)),
                n_units = length(sel),
                note = "no subject declared, so profile labels were permuted"))
  }
  u <- as.character(d$meta[[ucol]])[sel]
  g <- as.character(gf[sel])
  paired_u <- names(which(tapply(g, u, function(z) length(unique(z)) > 1)))
  if (length(paired_u) >= 2)
    return(list(scheme = "flip", sel = sel, ucol = ucol, units = u,
                flip_units = paired_u, n_units = length(paired_u),
                note = sprintf("the two labels were swapped within each of the %d %ss that carry both (sign-flip null)",
                               length(paired_u), ucol)))
  list(scheme = "permute", sel = sel, ucol = ucol, units = u, flip_units = character(0),
       n_units = length(unique(u)),
       note = sprintf("the label was permuted across %ss (between-unit null)", ucol))
}

xs_permute_study <- function(d, spec, plan) {
  g <- as.character(grp_factor(d, spec$groupvar))
  gs <- g[plan$sel]
  if (identical(plan$scheme, "flip")) {

    for (u in plan$flip_units) {
      ii <- which(plan$units == u)
      if (stats::runif(1) < 0.5) {
        z <- gs[ii]
        gs[ii] <- ifelse(z == spec$test, spec$reference,
                         ifelse(z == spec$reference, spec$test, z))
      }
    }
  } else if (identical(plan$scheme, "permute")) {
    uu <- unique(plan$units)
    lab <- vapply(uu, function(u) gs[which(plan$units == u)[1]], character(1))
    lab <- sample(lab)
    for (i in seq_along(uu)) gs[plan$units == uu[i]] <- lab[i]
  } else {
    gs <- sample(gs)
  }
  g[plan$sel] <- gs
  dd <- d

  dd$meta[["..xs_perm.."]] <- factor(g, levels = c(spec$reference, spec$test))
  dd
}

xs_flip_engine <- function(d, spec, plan, ref, tol = 1e-6) {
  no <- function(why) list(ok = FALSE, why = why)
  if (!identical(plan$scheme, "flip") || length(plan$flip_units) < 2)
    return(no("the design does not admit a sign-flip draw"))
  g  <- as.character(grp_factor(d, spec$groupvar))[plan$sel]
  u  <- plan$units
  cn <- colnames(d$expr)[plan$sel]
  E  <- d$expr[, plan$sel, drop = FALSE]
  uu <- plan$flip_units
  D <- vapply(uu, function(z) {
    i1 <- which(u == z & g == spec$test); i2 <- which(u == z & g == spec$reference)
    rowMeans(E[, i1, drop = FALSE], na.rm = TRUE) - rowMeans(E[, i2, drop = FALSE], na.rm = TRUE)
  }, numeric(nrow(E)))
  dimnames(D) <- list(rownames(E), uu)
  fc <- rowMeans(D, na.rm = TRUE)
  common <- intersect(names(fc), names(ref))
  if (length(common) < 10) return(no("the observed fit and the shortcut share too few features to check"))
  gap <- max(abs(fc[common] - ref[common]), na.rm = TRUE)
  if (!is.finite(gap) || gap > tol)
    return(no(sprintf("the shortcut does not reproduce the fitted coefficients (largest difference %.3g)", gap)))

  list(ok = TRUE, D = D, k = length(uu), max_gap = gap,
       n_distinct = 2^(length(uu) - 1L),
       fc = function(sgn) stats::setNames(as.numeric(D %*% sgn) / ncol(D), rownames(D)),
       enumerate = function() {
         k <- ncol(D)

         lapply(seq_len(2^(k - 1L)) - 1L, function(i)
           c(1, ifelse(bitwAnd(i, 2^(seq_len(k - 1L) - 1L)) > 0, -1, 1)))
       },
       draw = function() {
         s <- sample(c(-1, 1), ncol(D), replace = TRUE)
         stats::setNames(as.numeric(D %*% s) / ncol(D), rownames(D))
       })
}

xs_de_compare <- function(dA, dB, specA, specB, labA = "A", labB = "B",
                          fdr = 0.05, one_to_one_only = TRUE,
                          n_perm = 0L, seed = trace_seed(), min_shared = XS_MIN_SHARED,
                          max_exact = 4096L) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  no <- function(msg, ...) c(list(ok = FALSE, msg = msg), list(...))
  requested <- function(d, spec, label) list(label = label, spec = spec,
    text = xs_spec_text(spec), ok = FALSE, executed = FALSE,
    declared_subject = xs_decl(d, "subject"), declared_nest = xs_decl(d, "nest"),
    n_units = NA_integer_, unit = "not fitted")
  br <- xs_feature_bridge(dA, dB, one_to_one_only)
  if (br$n_shared < min_shared)
    return(no(sprintf("Only %d feature name is shared by the two layers; %d are needed for a genome-wide concordance. %s",
                      br$n_shared, min_shared, br$note), bridge = br,
              sideA = requested(dA, specA, labA), sideB = requested(dB, specB, labB)))

  rA <- xs_run_contrast(dA, specA, labA)
  rB <- xs_run_contrast(dB, specB, labB)
  if (!rA$ok) return(no(sprintf("%s: %s", labA, rA$declined), bridge = br, sideA = rA, sideB = rB))
  if (!rB$ok) return(no(sprintf("%s: %s", labB, rB$declined), bridge = br, sideA = rA, sideB = rB))

  rowsA <- xs_rows(br, "A"); rowsB <- xs_rows(br, "B")
  ia <- match(rowsA, rA$tt$Feature); ib <- match(rowsB, rB$tt$Feature)
  keep <- !is.na(ia) & !is.na(ib)
  if (sum(keep) < min_shared)
    return(no(sprintf("Only %d shared feature survived both models (each study drops features with no variance).",
                      sum(keep)), bridge = br, sideA = rA, sideB = rB))
  tab <- data.frame(
    Feature = br$features[keep],

    Feature_A = rowsA[keep], Feature_B = rowsB[keep],
    logFC_A = rA$tt$logFC[ia[keep]], logFC_B = rB$tt$logFC[ib[keep]],
    t_A = rA$tt$t[ia[keep]],         t_B = rB$tt$t[ib[keep]],
    P_A = rA$tt$P.Value[ia[keep]],   P_B = rB$tt$P.Value[ib[keep]],
    FDR_A = rA$tt$adj.P.Val[ia[keep]], FDR_B = rB$tt$adj.P.Val[ib[keep]],
    stringsAsFactors = FALSE)

  tab$FDR_shared_A <- stats::p.adjust(tab$P_A, "BH")
  tab$FDR_shared_B <- stats::p.adjust(tab$P_B, "BH")
  tab$sig_A <- tab$FDR_shared_A < fdr
  tab$sig_B <- tab$FDR_shared_B < fdr
  tab$same_sign <- sign(tab$logFC_A) == sign(tab$logFC_B) & tab$logFC_A != 0 & tab$logFC_B != 0
  tab$Agreement <- ifelse(tab$sig_A & tab$sig_B & tab$same_sign, "both, same direction",
                   ifelse(tab$sig_A & tab$sig_B, "both, opposite direction",
                   ifelse(tab$sig_A, paste(labA, "only"),
                   ifelse(tab$sig_B, paste(labB, "only"), "neither"))))

  r   <- suppressWarnings(stats::cor(tab$logFC_A, tab$logFC_B, use = "complete.obs"))
  rho <- suppressWarnings(stats::cor(tab$logFC_A, tab$logFC_B, use = "complete.obs", method = "spearman"))
  nA_sig <- sum(tab$sig_A); nB_sig <- sum(tab$sig_B)
  r_sig <- if (nA_sig >= 10)
    suppressWarnings(stats::cor(tab$logFC_A[tab$sig_A], tab$logFC_B[tab$sig_A], use = "complete.obs"))
    else NA_real_
  both <- sum(tab$sig_A & tab$sig_B & tab$same_sign)
  rep_rate  <- if (nA_sig) both / nA_sig else NA_real_
  sign_rate <- if (nA_sig) sum(tab$sig_A & tab$P_B < 0.05 & tab$same_sign) / nA_sig else NA_real_

  dir_n <- sum(tab$sig_A); dir_k <- sum(tab$sig_A & tab$same_sign)
  q2 <- table(factor(tab$sig_A, c(FALSE, TRUE)), factor(tab$sig_B, c(FALSE, TRUE)))
  fish <- tryCatch(stats::fisher.test(q2), error = function(e) NULL)

  perm <- NULL
  if (n_perm > 0L) {
    plan <- xs_exchange_plan(dB, specB)
    if (is.null(plan)) {
      perm <- list(ok = FALSE, msg = "The second study's contrast has too few profiles to permute.")
    } else {

      poss <- if (identical(plan$scheme, "flip")) 2^(plan$n_units - 1L)
              else factorial(min(plan$n_units, 10))

      dBs <- dB
      dBs$expr <- dB$expr[intersect(rowsB, rownames(dB$expr)), , drop = FALSE]
      ref <- stats::setNames(rB$tt$logFC, rB$tt$Feature)
      fast <- xs_flip_engine(dBs, specB, plan, ref)

      exact <- isTRUE(fast$ok) && fast$n_distinct <= max(n_perm, max_exact)
      null_r <- if (isTRUE(fast$ok) && exact) {
        vapply(fast$enumerate(), function(sgn) {
          fc <- fast$fc(sgn)
          j <- match(rowsB[keep], names(fc)); ok <- !is.na(j)
          if (sum(ok) < min_shared) return(NA_real_)
          suppressWarnings(stats::cor(tab$logFC_A[ok], fc[j[ok]], use = "complete.obs"))
        }, numeric(1))
      } else if (isTRUE(fast$ok)) {
        vapply(seq_len(n_perm), function(i) {
          fc <- fast$draw()
          j <- match(rowsB[keep], names(fc)); ok <- !is.na(j)
          if (sum(ok) < min_shared) return(NA_real_)
          suppressWarnings(stats::cor(tab$logFC_A[ok], fc[j[ok]], use = "complete.obs"))
        }, numeric(1))
      } else {
        vapply(seq_len(n_perm), function(i) {
          dp <- xs_permute_study(dBs, specB, plan)
          sp <- xs_spec("..xs_perm..", specB$test, specB$reference, specB$at_tp, specB$paired)
          z <- tryCatch(byod_run_limma(dp, sp$groupvar, sp$reference, sp$test, subject_col = xs_decl(dp, "subject"),
                                       paired = sp$paired, at_tp = sp$at_tp),
                        error = function(e) data.frame())
          if (!is.data.frame(z) || !nrow(z)) return(NA_real_)
          j <- match(rowsB[keep], z$Feature)
          ok <- !is.na(j)
          if (sum(ok) < min_shared) return(NA_real_)
          suppressWarnings(stats::cor(tab$logFC_A[ok], z$logFC[j[ok]], use = "complete.obs"))
        }, numeric(1))
      }
      null_r <- null_r[is.finite(null_r)]
      perm <- if (!length(null_r))
        list(ok = FALSE, msg = "No permutation could be fitted; the null is not available for this contrast.")
      else {

        pv <- if (exact) mean(abs(null_r) >= abs(r) - 1e-12)
              else (1 + sum(abs(null_r) >= abs(r))) / (length(null_r) + 1)
        list(ok = TRUE, null = null_r, n = length(null_r), scheme = plan$scheme,
             note = plan$note, n_units = plan$n_units,
             possible = poss, exact = exact,
             engine = if (isTRUE(fast$ok)) "closed form" else "refit",
             engine_check = if (isTRUE(fast$ok)) fast$max_gap else NA_real_,
             engine_why = if (isTRUE(fast$ok)) NULL else fast$why,
             p = pv,
             min_p = if (exact) 1 / length(null_r) else 1 / (min(n_perm, poss) + 1))
      }
    }
  }

  if (all(tab$Feature_A == tab$Feature) && all(tab$Feature_B == tab$Feature))
    tab$Feature_A <- tab$Feature_B <- NULL
  list(ok = TRUE, tab = tab[order(tab$P_A), ], bridge = br,
       sideA = rA, sideB = rB, labA = labA, labB = labB, fdr = fdr,
       n_shared = nrow(tab), r = r, rho = rho, r_sig = r_sig,
       n_sig_A = nA_sig, n_sig_B = nB_sig, n_both = both,
       rep_rate = rep_rate, sign_rate = sign_rate,
       pi1 = xs_pi1(tab$P_B[tab$sig_A]),
       dir_k = dir_k, dir_n = dir_n,
       quad = q2, or = if (is.null(fish)) NA_real_ else unname(fish$estimate),
       perm = perm)
}

xs_de_text <- function(o) {
  if (!isTRUE(o$ok)) return(o$msg)
  p <- function(x) if (!is.finite(x)) "n/a" else sprintf("%.0f%%", 100 * x)
  paste0(
    xs_exclusion_text(o),
    sprintf("%s (%s; %s, %d %s) and %s (%s; %s, %d %s) were fitted separately. ",
            o$labA, o$sideA$text, if (o$sideA$paired) "paired" else "unpaired",
            o$sideA$n_units, o$sideA$unit,
            o$labB, o$sideB$text, if (o$sideB$paired) "paired" else "unpaired",
            o$sideB$n_units, o$sideB$unit),
    sprintf("Across %s shared features (%s), the two effect-estimate vectors correlate r = %.2f (Spearman rho = %.2f)%s. ",
            format(o$n_shared, big.mark = ","), o$bridge$rule, o$r, o$rho,
            if (is.finite(o$r_sig))
              sprintf("; over the %s features %s called significant, r = %.2f",
                      format(o$n_sig_A, big.mark = ","), o$labA, o$r_sig) else ""),
    sprintf("%s called %s features at FDR < %.2g and %s called %s; %s reached that threshold in both with the same direction (%s of %s's). ",
            o$labA, format(o$n_sig_A, big.mark = ","), o$fdr, o$labB, format(o$n_sig_B, big.mark = ","),
            format(o$n_both, big.mark = ","), p(o$rep_rate), o$labA),
    if (o$dir_n > 0)
      sprintf("Of %s's significant features, %s have the same sign in %s. ",
              o$labA, p(o$dir_k / o$dir_n), o$labB) else "",
    if (is.finite(o$pi1))
      sprintf("Storey pi1 of %s's p-values on that set is %.2f. ", o$labB, o$pi1) else "",
    if (!is.null(o$perm) && isTRUE(o$perm$ok))
      sprintf("Under a null in which %s (%s), the observed correlation gives p = %.3g%s.",
              o$perm$note,
              if (isTRUE(o$perm$exact)) sprintf("all %s distinct rearrangements enumerated",
                                                format(o$perm$n, big.mark = ","))
              else sprintf("%d draws", o$perm$n),
              o$perm$p,
              if (o$perm$min_p > 0.01)
                sprintf(" (the design admits %s distinct rearrangements, so p cannot fall below %.3g)",
                        format(o$perm$possible, big.mark = ","), o$perm$min_p) else "")
    else "No permutation null was computed, so no p-value is attached to the correlation.")
}

XS_DE_HELP <- list(
  list("r over all shared features",
       "the correlation of the two effect-estimate vectors over every feature the two layers share. Features that changed in neither study are included, and in a genome-wide matrix they are the majority."),
  list("r over the discovery set",
       "the same correlation restricted to the features the first study called significant. Conditioned on the first study, so it is not symmetric, and it is not a test."),
  list("replication rate",
       "the proportion of the first study's significant features that also reach the FDR threshold in the second with the same sign. Depends on the second study's power as well as on agreement."),
  list("direction agreement",
       "the proportion of the first study's significant features that have the same sign in the second. A proportion, not a test: co-regulated features move together, so the signs are not independent trials and a binomial p against 0.5 would be far too small. Inference comes from the permutation p below, which preserves the correlation between features."),
  list("Storey pi1",
       "an estimate of the proportion of the first study's significant features that are non-null in the second, from the distribution of the second study's p-values on that set. Does not require them to cross a threshold."),
  list("permutation p",
       "the proportion of rearrangements of the second study's own labels that produce a cross-study correlation at least as large in absolute value as the observed one. The rearrangement respects the design: labels are swapped within unit for a within-unit contrast and permuted across units otherwise. Rearranging a paired contrast rescales the effect rather than removing it, and correlation is invariant to scale, so the null distribution is wide and the test has low power against a study with a large effect.")
)

xs_pathway_compare <- function(dA, dB, specA, specB, collection = "Hallmark",
                               mode = c("gsea", "ora"), labA = "A", labB = "B",
                               fdr = 0.05, top_n = 250L, seed = trace_seed()) {
  mode <- match.arg(mode)
  no <- function(msg, ...) c(list(ok = FALSE, msg = msg), list(...))
  if (!annotation_compatible(dA, dB, collection))
    return(no("The selected pathway collection is not compatible across these layers. Custom GMT comparison requires the same organism and identical set definitions."))
  rA <- xs_run_contrast(dA, specA, labA)
  rB <- xs_run_contrast(dB, specB, labB)
  if (!rA$ok) return(no(sprintf("%s: %s", labA, rA$declined), sideA = rA, sideB = rB))
  if (!rB$ok) return(no(sprintf("%s: %s", labB, rB$declined), sideA = rA, sideB = rB))

  one <- function(r, d) {
    if (identical(mode, "gsea")) {
      st <- stats::setNames(r$tt$t, r$tt$Feature)
      z <- run_gsea(st, collection, seed = seed, d = d)
      if (!nrow(z)) return(list(ok = FALSE, msg = attr(z, "declined") %||% "no gene set could be scored"))
      list(ok = TRUE, res = data.frame(pathway = z$pathway, score = z$NES,
                                       p = z$pval, padj = z$padj, size = z$size,
                                       stringsAsFactors = FALSE),
           score_name = "NES")
    } else {
      sig <- r$tt$Feature[r$tt$adj.P.Val < fdr]
      if (length(sig) < 3) return(list(ok = FALSE,
        msg = sprintf("ORA requires at least 3 features below FDR %.3g; this contrast has %d. Use GSEA for a threshold-free ranking.", fdr, length(sig))))
      z <- run_ora(sig, collection, universe = rownames(d$expr), d = d)
      if (!nrow(z)) return(list(ok = FALSE, msg = "no gene set was over-represented (fewer than 3 matched genes, or no overlap)"))

      lfc <- stats::setNames(r$tt$logFC, r$tt$Feature)
      sets <- annotation_sets(d, collection)
      key <- stats::setNames(names(sets), gsub("_", " ", names(sets)))
      dir <- vapply(z$pathway, function(pw) {
        g <- sets[[key[[pw]] %||% pw]]
        if (is.null(g)) return(1)
        v <- lfc[intersect(intersect(g, sig), names(lfc))]
        d <- sign(stats::median(v, na.rm = TRUE))
        if (!is.finite(d) || d == 0) 1 else d
      }, numeric(1))
      list(ok = TRUE, res = data.frame(pathway = z$pathway,
                                       score = unname(dir) * -log10(pmax(z$padj, 1e-300)),
                                       p = z$p, padj = z$padj, size = z$set_size,
                                       stringsAsFactors = FALSE),
           score_name = "signed -log10 FDR")
    }
  }
  eA <- one(rA, dA); eB <- one(rB, dB)
  if (!isTRUE(eA$ok)) return(no(sprintf("%s: %s", labA, eA$msg), sideA = rA, sideB = rB))
  if (!isTRUE(eB$ok)) return(no(sprintf("%s: %s", labB, eB$msg), sideA = rA, sideB = rB))

  m <- merge(eA$res, eB$res, by = "pathway", suffixes = c("_A", "_B"))
  if (nrow(m) < 5)
    return(no(sprintf("Only %d gene set was scored in both studies; too few to compare.", nrow(m)),
              sideA = rA, sideB = rB))
  names(m) <- sub("^score_", paste0(eA$score_name, "_"), names(m))
  sa <- paste0(eA$score_name, "_A"); sb <- paste0(eA$score_name, "_B")
  m$same_sign <- sign(m[[sa]]) == sign(m[[sb]])
  m$sig_A <- m$padj_A < fdr; m$sig_B <- m$padj_B < fdr
  m$Agreement <- ifelse(m$sig_A & m$sig_B & m$same_sign, "both, same direction",
                 ifelse(m$sig_A & m$sig_B, "both, opposite direction",
                 ifelse(m$sig_A, paste(labA, "only"),
                 ifelse(m$sig_B, paste(labB, "only"), "neither"))))
  r   <- suppressWarnings(stats::cor(m[[sa]], m[[sb]], use = "complete.obs"))
  rho <- suppressWarnings(stats::cor(m[[sa]], m[[sb]], use = "complete.obs", method = "spearman"))
  eith <- m$sig_A | m$sig_B
  agree_either <- if (any(eith)) mean(m$same_sign[eith]) else NA_real_
  m <- m[order(pmin(m$padj_A, m$padj_B)), ]
  list(ok = TRUE, tab = m, score_name = eA$score_name, mode = mode,
       collection = collection, labA = labA, labB = labB, fdr = fdr,
       sideA = rA, sideB = rB, sa = sa, sb = sb,
       n = nrow(m), r = r, rho = rho,
       n_sig_A = sum(m$sig_A), n_sig_B = sum(m$sig_B),
       n_both = sum(m$sig_A & m$sig_B & m$same_sign),
       n_either = sum(eith),
       agree_either = agree_either)
}

xs_pw_text <- function(o) {
  if (!isTRUE(o$ok)) return(o$msg)
  paste0(
    xs_exclusion_text(o),
    sprintf("Each study was scored against %s by %s using its own ranking, and only the resulting %s values were compared. ",
            o$collection, if (identical(o$mode, "gsea")) "GSEA" else "over-representation", o$score_name),
    sprintf("%d gene sets were scored in both. Their scores correlate r = %.2f (Spearman rho = %.2f). ",
            o$n, o$r, o$rho),
    sprintf("%d sets reached FDR < %.2g in %s and %d in %s; %d in both with the same direction",
            o$n_sig_A, o$fdr, o$labA, o$n_sig_B, o$labB, o$n_both),
    if (is.finite(o$agree_either))
      sprintf(", and among the %d significant in either study the direction agreed for %.0f%%.",
              o$n_either, 100 * o$agree_either) else ".",
    " Pathway scores do not require the two studies to share feature identifiers, so this comparison is available even when the feature-level one is not.")
}

xs_wrap <- function(x, width = 78L) {

  if (is.null(x) || !length(x)) return(x)
  x <- paste(as.character(x), collapse = " ")
  if (!nzchar(x)) return(x)
  paste(strwrap(x, width = width), collapse = "\n")
}
xs_legend_rows <- function(labs, per_row = 62L) {
  w <- sum(nchar(labs)) + 4L * length(labs)
  max(1L, min(length(labs), ceiling(w / per_row)))
}

xs_two_labels <- function(labA, labB) {
  if (!identical(labA, labB)) return(c(labA, labB))
  c(paste(labA, "(A)"), paste(labB, "(B)"))
}

xs_agree_levels <- function(labA, labB) {
  lb <- xs_two_labels(labA, labB)
  c("both, same direction", "both, opposite direction",
    paste(lb[1], "only"), paste(lb[2], "only"), "neither")
}

xs_de_scatter_gg <- function(o, cols = NULL, max_points = 8000L, label_n = 0L, seed = trace_seed()) {
  restore_rng <- trace_rng_scope(seed); on.exit(restore_rng(), add = TRUE)
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg))
  d <- o$tab
  lv <- xs_agree_levels(o$labA, o$labB)
  d$Agreement <- factor(d$Agreement, levels = lv)
  cc <- fb_cols(lv, cols, PAL_NAMES[1]); cc[is.na(cc)] <- "#C9D2DA"

  n_all <- nrow(d); thinned <- 0L
  if (nrow(d) > max_points) {
    dull <- which(d$Agreement == "neither")
    if (length(dull) > max_points) {
      d <- d[c(setdiff(seq_len(nrow(d)), dull), sample(dull, max_points)), , drop = FALSE]
      thinned <- n_all - nrow(d)
    }
  }
  lim <- max(abs(c(d$logFC_A, d$logFC_B)), na.rm = TRUE)
  g <- ggplot2::ggplot(d, ggplot2::aes(logFC_A, logFC_B)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#B9C4CE", linewidth = .3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#B9C4CE", linewidth = .3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#8FA0AE", linewidth = .4) +
    ggplot2::geom_point(ggplot2::aes(colour = Agreement), size = 1, alpha = .75) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                         colour = INK, linewidth = .5) +
    ggplot2::scale_colour_manual(values = cc, drop = FALSE, name = NULL,
                                 guide = ggplot2::guide_legend(nrow = xs_legend_rows(lv),
                                                               override.aes = list(size = 2.4))) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = xs_wrap(sprintf("%s  log2 FC", o$labA), 38),
                  y = xs_wrap(sprintf("%s  log2 FC", o$labB), 38),
                  title = sprintf("%s shared features \u00b7 r = %.2f, rho = %.2f", 
                                  format(o$n_shared, big.mark = ","), o$r, o$rho),
                  subtitle = paste(xs_wrap(paste("A:", o$sideA$text), 80),
                                   xs_wrap(paste("B:", o$sideB$text), 80), sep = "\n")) +
    ggplot2::theme_minimal(base_family = BYOD_UI_FONT, base_size = 11) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8),
                   plot.title.position = "plot",
                   plot.subtitle = ggplot2::element_text(size = 9, colour = MUTED,
                                                         margin = ggplot2::margin(b = 10)),
                   axis.title = ggplot2::element_text(size = 10),
                   plot.caption = ggplot2::element_text(size = 7.5, colour = MUTED, hjust = 0),
                   panel.grid.minor = ggplot2::element_blank())
  if (thinned > 0)
    g <- g + ggplot2::labs(caption = xs_wrap(sprintf(
      "%s of the %s features that were significant in neither study are drawn; every significant feature is. All %s were used for every number above.",
      format(sum(d$Agreement == "neither"), big.mark = ","),
      format(n_all - sum(d$Agreement != "neither"), big.mark = ","),
      format(n_all, big.mark = ",")), 104L))
  if (label_n > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    top <- utils::head(d[d$sig_A & d$sig_B, , drop = FALSE], label_n)
    if (nrow(top))
      g <- g + ggrepel::geom_text_repel(data = top, ggplot2::aes(label = Feature),
                                        size = 2.7, max.overlaps = 20, seed = seed)
  }
  g
}

xs_de_scatter <- function(o, cols = NULL) {
  if (!isTRUE(o$ok)) return(emptyplot(o$msg))
  d <- o$tab
  lv <- xs_agree_levels(o$labA, o$labB)
  cc <- fb_cols(lv, cols, PAL_NAMES[1]); cc[is.na(cc)] <- "#C9D2DA"
  d$col <- unname(cc[as.character(d$Agreement)])
  plot_ly(d, x = ~logFC_A, y = ~logFC_B, type = "scattergl", mode = "markers",
          marker = list(size = 4, color = ~col, opacity = .75),
          text = ~sprintf("%s<br>%s: %.2f (FDR %.3g)<br>%s: %.2f (FDR %.3g)",
                          Feature, o$labA, logFC_A, FDR_shared_A, o$labB, logFC_B, FDR_shared_B),
          hoverinfo = "text") %>%
    ply(xaxis = list(title = sprintf("%s log2 FC", o$labA), zeroline = TRUE),
        yaxis = list(title = sprintf("%s log2 FC", o$labB), zeroline = TRUE),
        title = list(text = sprintf("r = %.2f over all %s shared features%s",
                                    o$r, format(o$n_shared, big.mark = ","),
                                    if (is.finite(o$r_sig))
                                      sprintf("; r = %.2f over the %s significant in %s",
                                              o$r_sig, format(o$n_sig_A, big.mark = ","), o$labA) else ""),
                     x = .02, font = list(size = 13)))
}

xs_de_quadrant_gg <- function(o, cols = NULL) {
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg))
  lv <- xs_agree_levels(o$labA, o$labB)
  d <- data.frame(Agreement = factor(lv, levels = rev(lv)),
                  n = as.integer(table(factor(o$tab$Agreement, levels = lv))))
  d <- d[d$Agreement != "neither", , drop = FALSE]
  cc <- fb_cols(lv, cols, PAL_NAMES[1]); cc[is.na(cc)] <- "#C9D2DA"
  ggplot2::ggplot(d, ggplot2::aes(n, Agreement, fill = Agreement)) +
    ggplot2::geom_col(width = .68) +
    ggplot2::geom_text(ggplot2::aes(label = format(n, big.mark = ",")),
                       hjust = -0.15, size = 3.1, colour = INK) +
    ggplot2::scale_fill_manual(values = cc, guide = "none") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, .18))) +
    ggplot2::labs(x = sprintf("features (FDR < %.2g within the shared set)", o$fdr), y = NULL,
                  title = xs_wrap(sprintf("Features called in %s, by agreement with %s",
                                          o$labA, o$labB), 74),
                  subtitle = xs_wrap(sprintf("%s of %s reach FDR < %.2g in both with the same direction%s",
                                             format(o$n_both, big.mark = ","),
                                             format(o$n_sig_A, big.mark = ","), o$fdr,
                                             if (is.finite(o$rep_rate)) sprintf(" (%.0f%%)", 100 * o$rep_rate) else ""), 84)) +
    ggplot2::theme_minimal(base_family = BYOD_UI_FONT, base_size = 11) +
    ggplot2::theme(plot.title.position = "plot",
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank())
}

xs_de_null_gg <- function(o) {
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg))
  if (is.null(o$perm) || !isTRUE(o$perm$ok))
    return(emptyplot_gg(o$perm$msg %||% XS_PERM_OFF))
  d <- data.frame(r = o$perm$null)
  ggplot2::ggplot(d, ggplot2::aes(r)) +
    ggplot2::geom_histogram(bins = 30, fill = "#C9D2DA", colour = "white", linewidth = .2) +
    ggplot2::geom_vline(xintercept = o$r, colour = PIG, linewidth = .8) +
    ggplot2::annotate("text", x = o$r, y = Inf, label = sprintf("  observed r = %.2f", o$r),
                      hjust = if (o$r >= median(d$r)) 1 else 0, vjust = 1.6, size = 3.2, colour = PIG) +
    ggplot2::labs(x = xs_wrap(sprintf("cross-study r if %s had no effect (%s)", o$labB, o$perm$note), 80),
                  y = "permutations",
                  title = sprintf("%s \u00b7 p = %.3g",
                                  if (isTRUE(o$perm$exact))
                                    sprintf("all %s rearrangements", format(o$perm$n, big.mark = ","))
                                  else sprintf("%d draws", o$perm$n), o$perm$p),
                  caption = xs_wrap(sprintf("Null distribution of the cross-study correlation over %s of %s's own labels. Observed value marked.",
                                             if (isTRUE(o$perm$exact))
                                               sprintf("all %s distinct rearrangements",
                                                       format(o$perm$n, big.mark = ","))
                                             else sprintf("%d sampled rearrangements", o$perm$n),
                                             o$labB), 100)) +
    ggplot2::theme_minimal(base_family = BYOD_UI_FONT, base_size = 11) +
    ggplot2::theme(plot.title.position = "plot",
                   plot.caption = ggplot2::element_text(size = 7.5, colour = MUTED, hjust = 0),
                   panel.grid.minor = ggplot2::element_blank())
}

xs_de_null_plot <- function(o) {
  if (!isTRUE(o$ok)) return(emptyplot(o$msg))
  if (is.null(o$perm) || !isTRUE(o$perm$ok))
    return(emptyplot(o$perm$msg %||% XS_PERM_OFF))
  plot_ly(x = o$perm$null, type = "histogram", nbinsx = 30,
          marker = list(color = "#C9D2DA")) %>%
    ply(xaxis = list(title = "cross-study r under the null"), yaxis = list(title = "permutations"),
        shapes = list(list(type = "line", x0 = o$r, x1 = o$r, y0 = 0, y1 = 1,
                           yref = "paper", line = list(color = PIG, width = 2))),
        title = list(text = sprintf("observed r = %.2f · p = %.3g", o$r, o$perm$p),
                     x = .02, font = list(size = 13)))
}

xs_pw_scatter_gg <- function(o, cols = NULL, label_n = 12L) {
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg))
  d <- o$tab
  lv <- xs_agree_levels(o$labA, o$labB)
  d$Agreement <- factor(d$Agreement, levels = lv)
  cc <- fb_cols(lv, cols, PAL_NAMES[1]); cc[is.na(cc)] <- "#C9D2DA"
  d$x <- d[[o$sa]]; d$y <- d[[o$sb]]
  lim <- max(abs(c(d$x, d$y)), na.rm = TRUE)
  g <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_hline(yintercept = 0, colour = "#B9C4CE", linewidth = .3) +
    ggplot2::geom_vline(xintercept = 0, colour = "#B9C4CE", linewidth = .3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "#8FA0AE", linewidth = .4) +
    ggplot2::geom_point(ggplot2::aes(colour = Agreement, size = pmin(size_A, 500)), alpha = .85) +
    ggplot2::scale_colour_manual(values = cc, drop = FALSE, name = NULL,
                                 guide = ggplot2::guide_legend(nrow = xs_legend_rows(lv),
                                                               override.aes = list(size = 2.4))) +
    ggplot2::scale_size_continuous(range = c(1.2, 4.5), guide = "none") +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    ggplot2::labs(x = xs_wrap(sprintf("%s  %s", o$labA, o$score_name), 38),
                  y = xs_wrap(sprintf("%s  %s", o$labB, o$score_name), 38),
                  title = xs_wrap(sprintf("%s · %d sets scored in both · r = %.2f, rho = %.2f",
                                          o$collection, o$n, o$r, o$rho), 82)) +
    ggplot2::theme_minimal(base_family = BYOD_UI_FONT, base_size = 11) +
    ggplot2::theme(legend.position = "bottom", legend.text = ggplot2::element_text(size = 8),
                   plot.title.position = "plot",
                   panel.grid.minor = ggplot2::element_blank())
  if (label_n > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
    top <- utils::head(d[d$sig_A | d$sig_B, , drop = FALSE], label_n)
    if (nrow(top))
      g <- g + ggrepel::geom_text_repel(data = top, ggplot2::aes(label = pathway),
                                        size = 2.5, max.overlaps = 24, seed = trace_seed())
  }
  g
}

xs_pw_scatter <- function(o, cols = NULL) {
  if (!isTRUE(o$ok)) return(emptyplot(o$msg))
  d <- o$tab
  lv <- xs_agree_levels(o$labA, o$labB)
  cc <- fb_cols(lv, cols, PAL_NAMES[1]); cc[is.na(cc)] <- "#C9D2DA"
  d$col <- unname(cc[as.character(d$Agreement)])
  plot_ly(d, x = d[[o$sa]], y = d[[o$sb]], type = "scatter", mode = "markers",
          marker = list(size = 8, color = d$col, opacity = .85),
          text = ~sprintf("%s<br>%s: %.2f (FDR %.3g)<br>%s: %.2f (FDR %.3g)",
                          pathway, o$labA, d[[o$sa]], padj_A, o$labB, d[[o$sb]], padj_B),
          hoverinfo = "text") %>%
    ply(xaxis = list(title = sprintf("%s %s", o$labA, o$score_name), zeroline = TRUE),
        yaxis = list(title = sprintf("%s %s", o$labB, o$score_name), zeroline = TRUE),
        title = list(text = sprintf("%s · r = %.2f · rho = %.2f", o$collection, o$r, o$rho),
                     x = .02, font = list(size = 13)))
}

xs_pw_bar_gg <- function(o, n = 20L, cols = NULL) {
  if (!isTRUE(o$ok)) return(emptyplot_gg(o$msg))
  d <- utils::head(o$tab, n)
  lb <- xs_two_labels(o$labA, o$labB)
  long <- rbind(
    data.frame(pathway = d$pathway, study = lb[1], score = d[[o$sa]], stringsAsFactors = FALSE),
    data.frame(pathway = d$pathway, study = lb[2], score = d[[o$sb]], stringsAsFactors = FALSE))
  ord <- d$pathway[order(d[[o$sa]] + d[[o$sb]])]
  long$pathway <- factor(long$pathway, levels = ord)
  long$study <- factor(long$study, levels = lb)
  cc <- fb_cols(lb, cols, PAL_NAMES[1])
  ggplot2::ggplot(long, ggplot2::aes(score, pathway, fill = study)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = .72), width = .68) +
    ggplot2::geom_vline(xintercept = 0, colour = "#8FA0AE", linewidth = .3) +
    ggplot2::scale_fill_manual(values = cc, name = NULL) +
    ggplot2::labs(x = o$score_name, y = NULL,
                  title = xs_wrap(sprintf("%s · top %d sets by best FDR", o$collection, nrow(d)), 80)) +
    ggplot2::theme_minimal(base_family = BYOD_UI_FONT, base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_text(size = 7.5))
}

xs_design_gg <- function(tab, labA = "A", labB = "B") {
  if (is.null(tab) || !nrow(tab)) return(emptyplot_gg("no design comparison"))
  d <- data.frame(Property = factor(tab$Property, levels = rev(tab$Property)),
                  A = tab[[2]], B = tab[[3]], Verdict = tab$Verdict,
                  stringsAsFactors = FALSE)
  vc <- c(ok = "#2C8C99", informational = "#6E7E8D", check = "#E0A03A",
          blocking = "#B0413E", "feature-level blocked" = "#B0413E",
          "time-contrast blocked" = "#E0A03A")
  lb <- xs_two_labels(labA, labB)
  long <- rbind(
    data.frame(Property = d$Property, side = lb[1], value = d$A, Verdict = d$Verdict),
    data.frame(Property = d$Property, side = lb[2], value = d$B, Verdict = d$Verdict))
  long$side <- factor(long$side, levels = lb)
  long$value <- gsub("\u2014", "not declared", long$value, fixed = TRUE)
  ggplot2::ggplot(long, ggplot2::aes(side, Property)) +
    ggplot2::geom_tile(ggplot2::aes(fill = Verdict), colour = "white", linewidth = 1, alpha = .16) +
    ggplot2::geom_text(ggplot2::aes(label = value, colour = Verdict), size = 2.9) +
    ggplot2::scale_fill_manual(values = vc, guide = "none") +
    ggplot2::scale_colour_manual(values = vc, name = NULL) +
    ggplot2::labs(x = NULL, y = NULL, title = "Design comparison") +
    ggplot2::theme_minimal(base_family = BYOD_UI_FONT, base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(face = "bold"))
}
