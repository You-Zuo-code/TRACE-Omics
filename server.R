function(input, output, session) {
  session_started <- Sys.time()

  analysis_seed <- reactiveVal(TRACE_DEFAULT_SEED)
  session$userData$analysis_seed <- analysis_seed
  observeEvent(input$apply_seed, {
    value <- tryCatch(trace_validate_seed(input$analysis_seed), error = function(e) e)
    if (inherits(value, "error")) {
      showNotification(conditionMessage(value), type = "error"); return()
    }
    if (!identical(value, analysis_seed())) {
      analysis_seed(value)
    }
  })

  upload_paths <- character()
  observe({
    paths <- unlist(lapply(c("byod_expr1", "byod_expr2", "byod_expr3", "byod_meta",
      paste0("byod_mapping", 1:3), paste0("byod_gmt", 1:3)),
      function(id) input[[id]]$datapath), use.names = FALSE)
    upload_paths <<- unique(c(upload_paths, paths))
  })
  session$onSessionEnded(function() {
    tmp <- paste0(normalizePath(tempdir()), .Platform$file.sep)
    for (f in upload_paths) {
      if (!is.na(f) && file.exists(f) && startsWith(normalizePath(f), tmp))
        unlink(f, force = TRUE)
    }
    br <- tryCatch(session$userData$bundle_root, error = function(e) NULL)
    if (!is.null(br) && nzchar(br) && dir.exists(br) &&
        startsWith(normalizePath(br), tmp)) unlink(br, recursive = TRUE, force = TRUE)
  })

  analysis_log <- reactiveVal(list())   # list of records: list(module, params=list(...), t)
  MODULE_ENGINE <- c(
    "DE (limma)"               = "limma",
    "Group x time interaction" = "limma",
    "Subject variables: features" = "stats",
    "Enrichment"               = "fgsea",
    "Pathways"                 = "GSVA",
    "WGCNA"                    = "WGCNA",
    "DIABLO"                   = "mixOmics",
    "Integration"              = "mixOmics",
    "Concordance"              = "stats",
    "Temporal programs"        = "limma",
    "Clinical correlates"      = "limma",
    "DIABLO cross-validation"  = "mixOmics",
    "Subject variables: screen" = "limma",
    "Subject variables: pathways" = "stats",
    "Multi-timepoint"          = "limma",
    "Trajectory consistency"   = "stats",
    "Baseline sensitivity"     = "limma",
    "Response composition"     = "limma",
    "Landscape"                = "stats",
    "Heatmap"                  = "stats",
    "Module score"             = "stats",
    "Trajectory"               = "stats",
    "Association"              = "stats",
    "Cross-study comparison"   = "stats",

    "cross-dataset differential concordance" = "limma",
    "cross-dataset pathway concordance"      = "fgsea")
  .pkg_ver <- function(pkg) tryCatch(as.character(utils::packageVersion(pkg)),
                                     error = function(e) NA_character_)
  .design_context <- function() {
    s <- tryCatch(byod_study(), error = function(e) NULL)
    if (is.null(s)) return(list())
    ln <- isolate(input$byod_layer) %||% s$layer_names[1]
    if (!ln %in% s$layer_names) ln <- s$layer_names[1]
    d  <- s$layers[[ln]]
    u  <- d$nest %||% d$subject

    slot <- tryCatch(byod_aslot(), error = function(e) NULL)
    out <- list(study = s$name, study_id = bundle_study_id(s),
                dataset = if (is.null(slot)) "built-in cohort" else xs_slot_label(slot),
                declared_subject = d$subject %||% "(none declared)",
                declared_nest    = d$nest %||% "(none declared)",
                annotation = annotation_summary(d))
    if (!is.null(u) && u %in% names(d$meta))
      out$study_units <- sprintf("%d %s", length(unique(d$meta[[u]])), u)
    out$study_profiles <- ncol(d$expr)
    out
  }
  log_event <- function(module, ..., outcome = "completed", engine = NULL) {
    eng <- if (!is.null(engine)) engine else unname(MODULE_ENGINE[module])
    if (!length(eng) || anyNA(eng))
      warning(sprintf("no engine registered for module '%s'; its record will be incomplete",
                      module), call. = FALSE)
    prov <- .design_context()
    if (length(eng) && !anyNA(eng))
      prov <- c(list(engine = paste(vapply(eng, function(p)
        sprintf("%s %s", p, .pkg_ver(p)), character(1)), collapse = "; ")), prov)
    if (identical(outcome, "declined")) {
      prov$requested_engine <- prov$engine
      prov$engine <- "not completed; see outcome and per-side status"
    }
    prov <- c(prov, list(outcome = outcome, random_seed = trace_seed(),
                         R = as.character(getRversion())))
    rec <- list(module = module, params = c(list(...), prov),
                t = format(Sys.time(), "%H:%M:%S"))
    isolate(analysis_log(c(analysis_log(), list(rec))))
  }

  log_declined <- function(module, reason, ...)
    log_event(module, ..., reason = reason, outcome = "declined")

  .mode_engine <- function(mode) if (identical(mode, "ora")) "stats" else "fgsea"
  .log_line <- function(r) {
    kv <- if (length(r$params)) paste(sprintf("%s=%s", names(r$params),
              vapply(r$params, function(x) paste(as.character(x), collapse = ","), character(1))), collapse = "; ") else ""
    sprintf("[%s] %-16s %s", r$t, r$module, kv)
  }
  output$repro_log <- renderText({
    lg <- analysis_log()
    if (!length(lg)) "No analyses run yet. Every parameterised run (DE, enrichment, integration, pathways, WGCNA, DIABLO, subject-variable association, group×time, temporal programs) is logged here for reproducibility." else paste(vapply(lg, .log_line, character(1)), collapse = "\n")
  })

  .log_studies <- reactive({
    z <- unlist(lapply(analysis_log(), function(r) r$params$study))
    unique(z[!is.na(z) & nzchar(z)])
  })
  output$repro_methods <- renderText(byod_methods_text(analysis_log(), .log_studies()))
  output$repro_methods_dl <- downloadHandler(
    filename = function() paste0("TRACE-Omics_methods_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) writeLines(byod_methods_text(isolate(analysis_log()), isolate(.log_studies())), file))
  output$repro_env <- renderText({
    paste0(
      "R ", getRversion(), "  |  platform ", R.version$platform, "\n",
      "shiny ", as.character(utils::packageVersion("shiny")),
      " · limma ", as.character(utils::packageVersion("limma")),
      " · fgsea ", as.character(utils::packageVersion("fgsea")), "\n",
      "GSVA ", if (requireNamespace("GSVA", quietly=TRUE)) as.character(utils::packageVersion("GSVA")) else "not installed (precomputed used)",
      " · WGCNA ", if (requireNamespace("WGCNA", quietly=TRUE)) as.character(utils::packageVersion("WGCNA")) else "not installed",
      " · mixOmics ", if (requireNamespace("mixOmics", quietly=TRUE)) as.character(utils::packageVersion("mixOmics")) else "not installed",
      "\nSession start: ", format(session_started, "%Y-%m-%d %H:%M"))
  })
  output$repro_dl <- downloadHandler(
    filename = function() paste0("TRACE-Omics_session_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"),
    content = function(file) {
      writeLines(c("TRACE-Omics — analysis session log",
        paste("Generated:", format(Sys.time())),
        paste("R version:", getRversion(), "| platform:", R.version$platform),
        strrep("-", 60),
        vapply(analysis_log(), .log_line, character(1)),
        strrep("-", 60),
        "Note: live analyses (limma DE, GSEA/ORA, GSVA, WGCNA, DIABLO) are",
        "convenient re-analyses and may differ slightly from the companion",
        "paper's full pipeline. See Methods."), file)
    })
  observeEvent(input$repro_reset, {
    analysis_log(list()); byod_msg(list(type = "info", txt = "Session log cleared."))
  })
  output$repro_runtime <- renderDT({
    f <- file.path("data", "runtime_table.csv")
    if (!file.exists(f)) return(NULL)
    datatable(utils::read.csv(f, check.names = FALSE), rownames = FALSE,
      options = list(pageLength = 12, dom = "t"), class = "compact stripe")
  })

  feat_union <- function(sp) sort(union(features_of(sp, "mRNA"), features_of(sp, "Protein")))
  parse_genes <- function(txt) { t <- unlist(strsplit(txt %||% "", "[,;[:space:]]+")); t[nzchar(t)] }
  match_feats <- function(d, toks) {
    idx <- match(toupper(toks), toupper(rownames(d$expr))); rownames(d$expr)[idx[!is.na(idx)]] }
  `%||%` <- function(a, b) if (is.null(a)) b else a
  def_gene <- function(ch, pref = "C3") if (pref %in% ch) pref else ch[1]

  cs_ch <- sort(intersect(feat_union("Pig"), feat_union("Human")))
  updateSelectizeInput(session, "cs_gene", choices = cs_ch,
                       selected = def_gene(cs_ch, "C3"), server = TRUE)

  ov_card <- function(s) {
    accent <- if (s$species == "Pig") PIG else HUMAN
    div(class = "card", style = "padding:16px 18px;height:100%;",
      div(style = paste0("font-size:.72rem;font-weight:700;letter-spacing:1px;color:", accent, ";"),
          toupper(s$species)),
      div(style = "font-weight:600;color:#1E3A5F;margin-bottom:6px;", paste0(s$omic, " expression")),
      div(class = "statnum", format(s$n_features, big.mark = ",")),
      div(class = "statlab", "features"),
      div(style = "margin-top:6px;", span(class = "statnum", style = "font-size:1.3rem;", s$n_kidneys),
          span(class = "statlab", " kidneys · ", length(s$tp), " timepoints")))
  }
  output$ov_card_1 <- renderUI(ov_card(atlas$summary$pig_mrna))
  output$ov_card_2 <- renderUI(ov_card(atlas$summary$pig_protein))
  output$ov_card_3 <- renderUI(ov_card(atlas$summary$hum_mrna))
  output$ov_card_4 <- renderUI(ov_card(atlas$summary$hum_protein))

  output$ov_design <- renderPlotly({
    df <- do.call(rbind, lapply(names(atlas$summary), function(k) {
      s <- atlas$summary[[k]]
      data.frame(ds = paste(s$species, s$omic), tp = names(s$tp), n = as.integer(s$tp)) }))
    df$tp <- factor(df$tp, levels = atlas$tp_pig)
    plot_ly(df, x = ~tp, y = ~n, color = ~ds,
      colors = c("Pig mRNA" = PIG, "Pig Protein" = "#D98A88",
                 "Human mRNA" = HUMAN, "Human Protein" = "#8FB4DA"),
      type = "bar") %>%
      ply(barmode = "group", xaxis = list(title = ""), yaxis = list(title = "samples"),
          legend = list(orientation = "h", y = 1.15))
  })

  output$about_txt <- renderUI(HTML(paste0(
    "<p><b>TRACE-Omics</b> is a platform for longitudinal and paired multi-omics studies &mdash; ",
    "designs in which the same organ, patient or culture is measured more than once, so the ",
    "observations are <b>not independent</b>. You declare the design once: the layers and their ",
    "measurement units, the ordered time, the repeatedly sampled subject, any higher-level unit ",
    "that subject nests inside, and the grouping variable. From that declaration the engine ",
    "resolves the <b>independent unit for every contrast</b>, the level at which pairing exists, ",
    "and whether two layers join per sample or only per condition.</p>",
    "<p><b>Analysis modules</b> are organised in three blocks &mdash; <i>single-omics analysis</i> ",
    "within one layer, <i>multi-omics integration</i> across two layers, and <i>cross-study ",
    "comparison</i> between two studies. A module is offered only where the resolved design ",
    "supports it, and declines with a stated reason otherwise. Eight assay types are accepted: ",
    "transcriptomics, proteomics, miRNA, metabolomics, lipidomics, methylation, chromatin ",
    "accessibility, and any other numeric feature &times; sample matrix.</p>",
    "<p><b>Built-in demonstration.</b> The <i>Kidney Perfusion Reference</i> ships with the app: ",
    "cortical biopsies from porcine kidney pairs ",
    "(<a href=\"https://pre-image.eu/prepair-1-2/\" target=\"_blank\" rel=\"noopener\">PREPAIR-1</a>) and ",
    "discarded human kidneys (<a href=\"https://pre-image.eu/prepair-2/\" target=\"_blank\" rel=\"noopener\">PREPAIR-2</a>), ",
    "sampled at several timepoints during normothermic machine perfusion, with paired mRNA and ",
    "protein on overlapping samples. Porcine kidney pairs contain a warm-ischaemia kidney and ",
    "its control, so treatment varies <i>within</i> animal ",
    "and the animal is the upper-level unit. Uploads and built-in cohorts are the same object ",
    "type driven by the same engine, so every module behaves identically on either.</p>",
    "<p><b>Processing.</b> The built-in mRNA layers use log2(TPM+1); protein uses log2 intensity. ",
    "Uploaded matrices retain their declared scale after the selected preprocessing.</p>",
    "<p><b>Statistics.</b> Contrasts are computed <i>live</i> with limma at the resolved unit ",
    "(moderated t-test; within-subject differences averaged inside any declared nesting level ",
    "before testing across units). mRNA&ndash;protein and cross-species correlations use ",
    "descriptive correlation summaries; module score = mean per-gene z-score. Analysis records ",
    "include selected settings, available model and unit information, and execution outcomes. ",
    "Tables, figures and methods text can be exported.</p>",
    "<p><b>Developed by</b> You Zuo and Renhui Guan in the group of Cyril Moers, ",
    "University Medical Center Groningen &mdash; part of the ",
    "<a href=\"https://pre-image.eu\" target=\"_blank\" rel=\"noopener\">PRE-IMAGE</a> project (ERC grant 851368).</p>",
    "<p><b>Built with</b> ",
    "<a href=\"https://shiny.posit.co\" target=\"_blank\" rel=\"noopener\">R / Shiny</a>, ",
    "bslib, plotly, DT, limma, fgsea, GSVA, WGCNA and mixOmics. Released under the MIT ",
    "licence.</p>")))

  cs_align <- reactive({
    dP <- get_ds("Pig", input$cs_om %||% "mRNA"); dH <- get_ds("Human", input$cs_om %||% "mRNA")
    al <- cs_align_timepoints(dP, dH, "label")
    if (al$n < 2) al <- cs_align_timepoints(dP, dH, "position")
    al
  })
  observeEvent(cs_align(), {
    al <- cs_align(); if (al$n < 2) return()
    ch <- al$label[seq_len(al$n - 1L)]
    cur <- isolate(input$cs_base)
    updateSelectInput(session, "cs_base", choices = ch,
                      selected = if (!is.null(cur) && cur %in% ch) cur else ch[1])
  }, ignoreNULL = FALSE)

  cs_base <- reactive({
    req(input$cs_gene)
    disp <- input$cs_display %||% "raw"
    if (identical(disp, "relative")) req(input$cs_base)
    plot_crossspecies(input$cs_om, input$cs_gene, display = disp,
                      baseline = input$cs_base, measure = input$cs_measure %||% "difference")
  })
  cs_levels <- reactive({ s <- cs_base()$summary; if (is.null(s) || !nrow(s)) "all" else unique(as.character(s$Group)) })
  cs_export <- reactive({ o <- cs_base()
    list(summary = o$summary, ylab = o$ylab, ttl = o$ttl, splitby = NULL,
         prefix = paste0(input$cs_gene, "_crossspecies",
                         if (identical(o$display, "relative")) paste0("_rel_", input$cs_measure) else "")) })
  cs_cols <- apprServer("cs_appr", cs_levels, cs_export)
  output$cs_plot <- renderPlotly({ o <- cs_base(); if (is.null(o$summary) || !nrow(o$summary)) return(o$plot)
    p <- tc_plotly(list(summary = o$summary, splits = cs_levels()), o$ylab, o$ttl, cs_cols())
    if (identical(o$display, "relative"))
      p <- plotly::layout(p, shapes = list(list(type = "line", xref = "paper", x0 = 0, x1 = 1,
        yref = "y", y0 = 0, y1 = 0, line = list(color = MUTED, width = 1, dash = "dot"))))
    p })

  output$cs_display_help <- renderUI({
    o <- cs_base()
    if (identical(input$cs_display %||% "raw", "raw"))
      return(helpText("Each cohort's own expression level. ",
        tags$b("Compare the shape of the change, not the absolute levels:"),
        " the two cohorts differ in platform, batch and baseline state, so a difference ",
        "in height between the two curves is not interpretable. Switch to ",
        tags$b("Relative to baseline"), " to remove that offset and put both cohorts on one scale."))
    tagList(
      helpText("Each cohort against ", tags$b("its own"), " baseline, computed per animal or donor ",
               "and then averaged, so the error bars still describe variation between independent ",
               "units. The baseline sits at zero by definition."),
      if (!is.null(o$note)) div(class = "alert alert-warning py-2 px-3 small", o$note))
  })

  output$cs_measure_help <- renderUI({
    si <- cs_scale_info(get_ds("Pig", input$cs_om %||% "mRNA"))
    ad <- cs_measure_advice(si)[[input$cs_measure %||% "difference"]]
    cls <- c(recommended = "text-success", usable = "text-muted",
             discouraged = "text-warning", "not meaningful" = "text-danger")
    div(class = paste("small", cls[[ad$status]] %||% "text-muted"),
        tags$b(paste0(ad$status, ". ")), ad$why)
  })

  csr_registry <- reactive({
    reg <- lapply(names(STUDIES), function(k) list(key = k, name = STUDIES[[k]]$name,
                                                   study = STUDIES[[k]], source = "built-in"))
    ws <- byod_ws()
    for (k in xs_slot_keys()) {
      r <- ws[[k]]; if (is.null(r)) next
      reg <- c(reg, list(list(key = k, name = sprintf("%s (%s)", r$study$name, xs_slot_label(k)),
                              study = r$study,
                              source = if (is.null(r$study$reference_key)) "upload" else "reference-example")))
    }
    setNames(reg, vapply(reg, `[[`, character(1), "key"))
  })

  observeEvent(csr_registry(), {
    reg <- csr_registry()
    ch <- setNames(names(reg), vapply(reg, `[[`, character(1), "name"))
    keep <- function(cur, dflt) if (!is.null(cur) && cur %in% ch) cur else dflt
    updateSelectInput(session, "csr_A", choices = ch,
                      selected = keep(isolate(input$csr_A), ch[[1]]))
    updateSelectInput(session, "csr_B", choices = ch,
                      selected = keep(isolate(input$csr_B), ch[[min(2, length(ch))]]))
  }, ignoreNULL = TRUE)

  csr_layers_of <- function(key) { reg <- csr_registry(); if (is.null(reg[[key]])) NULL else reg[[key]]$study$layer_names }
  observeEvent(input$csr_A, {
    ln <- csr_layers_of(input$csr_A); req(length(ln))
    cur <- isolate(input$csr_layA)
    updateSelectInput(session, "csr_layA", choices = ln,
                      selected = if (!is.null(cur) && cur %in% ln) cur else ln[1])
  })
  observeEvent(input$csr_B, {
    ln <- csr_layers_of(input$csr_B); req(length(ln))
    cur <- isolate(input$csr_layB)
    updateSelectInput(session, "csr_layB", choices = ln,
                      selected = if (!is.null(cur) && cur %in% ln) cur else ln[1])
  })

  csr_dA <- reactive({ req(input$csr_A, input$csr_layA); reg <- csr_registry()
    req(reg[[input$csr_A]]); reg[[input$csr_A]]$study$layers[[input$csr_layA]] })
  csr_dB <- reactive({ req(input$csr_B, input$csr_layB); reg <- csr_registry()
    req(reg[[input$csr_B]]); reg[[input$csr_B]]$study$layers[[input$csr_layB]] })
  csr_labA <- reactive({ reg <- csr_registry(); req(reg[[input$csr_A]]); reg[[input$csr_A]]$study$name })
  csr_labB <- reactive({ reg <- csr_registry(); req(reg[[input$csr_B]]); reg[[input$csr_B]]$study$name })
  csr_idA <- reactive(bundle_study_id(csr_sA()))
  csr_idB <- reactive(bundle_study_id(csr_sB()))

  xs_ortho_ok <- reactive({ dA <- csr_dA(); dB <- csr_dB(); req(dA, dB)
    xs_ortholog_applicable(dA, dB) })
  output$csr_ortho_ui <- renderUI({
    a <- tryCatch(xs_ortho_ok(), error = function(e) NULL); req(a)
    if (isTRUE(a$ok))
      return(tagList(
        if (isTRUE(csr_dA()$atlas_annotation) && isTRUE(csr_dB()$atlas_annotation))
          checkboxInput("csr_1to1", "Restrict the atlas comparison to curated 1:1 orthologs", TRUE),
        div(class = "small text-muted", sprintf("%s. Uploaded mappings always exclude ambiguous and duplicate-target features.", a$why))))
    div(class = "small text-muted",
        tags$b("Feature matching. "), a$why, ". See the correspondence summary for the applied rule.")
  })

  output$csr_upload_hint <- renderUI({
    reg <- csr_registry()
    up <- Filter(function(z) identical(z$source, "upload"), reg)
    if (length(up) >= 2)
      return(div(class = "small text-success mt-1",
                 sprintf("%d uploaded studies are available, so both sides of this comparison can be your own: %s.",
                         length(up), paste(vapply(up, `[[`, character(1), "name"), collapse = ", "))))
    if (length(up) == 1)
      return(div(class = "small text-warning mt-1",
        "One uploaded study is loaded. To compare two studies of your own \u2014 yours against a published ",
        "series, say \u2014 build the second into another dataset slot under ", tags$b("Your data"),
        ". Until then either side can still be a built-in cohort.", br(),
        actionLink("csr_goto_upload", "Take me there", class = "small")))
    if (any(vapply(reg, function(z) identical(z$source, "reference-example"), logical(1))))
      return(div(class = "small text-muted mt-1",
        "Reference examples are workspace views of the built-in human cohort, not independent studies. ",
        "Upload separate studies under Your data to compare independent datasets."))
    div(class = "small text-muted mt-1",
        "Only the built-in cohorts are loaded. Under ", tags$b("Your data"),
        " you can build several studies at once \u2014 each with its own design declaration \u2014 and every one ",
        "of them appears in both pickers here.", br(),
        actionLink("csr_goto_upload", "Take me there", class = "small"))
  })
  observeEvent(input$csr_goto_upload, {
    nav_select("main_nav", "Your data", session = session)
  })

  csr_align_obj <- reactive({
    dA <- csr_dA(); dB <- csr_dB(); req(dA, dB)

    xs_align_timepoints(dA, dB, input$csr_align %||% "label",
                        manual = if (identical(input$csr_align, "manual")) xs_manual_map() else NULL)
  })

  observe({
    dA <- csr_dA(); dB <- csr_dB(); req(dA, dB)
    if (!length(cs_align_timepoints(dA, dB, "label")$A) &&
        identical(isolate(input$csr_align), "label"))
      updateRadioButtons(session, "csr_align", selected = "position")
  })

  output$csr_align_help <- renderUI({
    al <- csr_align_obj(); req(al)
    if (identical(al$align, "position") && al$n >= 2)
      return(div(class = "alert alert-warning py-2 px-2 small mb-2",
                 tags$b("Positions are being treated as comparable. "), al$note,
                 tags$br(), "Pairs: ", paste(al$label, collapse = "; ")))
    if (!al$n) return(div(class = "alert alert-danger py-2 px-2 small mb-2", al$note))
    div(class = "small text-muted mb-2",
        sprintf("%d timepoints share a label: %s", al$n, paste(al$A, collapse = ", ")))
  })

  observe({
    al <- csr_align_obj(); req(al)
    if (al$n < 2) { updateSelectInput(session, "csr_base", choices = character(0)); return() }
    ch <- head(al$label, -1)
    cur <- isolate(input$csr_base)
    updateSelectInput(session, "csr_base", choices = ch,
                      selected = if (!is.null(cur) && cur %in% ch) cur else ch[1])
  })

  csr_scales <- reactive({ dA <- csr_dA(); dB <- csr_dB(); req(dA, dB)
    list(A = cs_scale_info(dA), B = cs_scale_info(dB)) })

  observe({
    sc <- csr_scales(); req(sc)
    adA <- cs_measure_advice(sc$A); adB <- cs_measure_advice(sc$B)
    worst <- function(m) {
      rank <- c("not meaningful" = 1, "discouraged" = 2, "usable" = 3, "recommended" = 4)
      a <- adA[[m]]$status; b <- adB[[m]]$status
      if (rank[[a]] <= rank[[b]]) a else b
    }
    base <- c(difference = "Difference on the log scale",
              ratio = "Ratio to baseline (log2)",
              percent = "Percent change from baseline")

    mark <- c("recommended" = "  [recommended]", "usable" = "",
              "discouraged" = "  [discouraged for this data]",
              "not meaningful" = "  [not meaningful for this data]")
    lab <- vapply(CS_MEASURES, function(m) paste0(base[[m]], mark[[worst(m)]]), character(1))
    ch <- setNames(CS_MEASURES, unname(lab))
    cur <- isolate(input$csr_measure)
    best <- CS_MEASURES[which.max(vapply(CS_MEASURES, function(m)
      c("not meaningful" = 1, "discouraged" = 2, "usable" = 3, "recommended" = 4)[[worst(m)]], numeric(1)))]
    updateSelectInput(session, "csr_measure", choices = ch,
                      selected = if (!is.null(cur) && cur %in% CS_MEASURES) cur else best)
  })

  output$csr_scale_banner <- renderUI({
    sc <- csr_scales(); req(sc)
    row <- function(lab, si) tags$li(tags$b(lab), ": ", si$label,

      tags$span(class = "text-muted",
                sprintf(" (declared as '%s', %s; range %.2f to %.2f)",
                        si$declared, si$source, si$min, si$max)))
    same <- sc$A$family == sc$B$family
    div(class = paste("alert", if (same) "alert-light" else "alert-warning", "py-2 px-3 mb-3"),
        tags$b("Measurement scale detected"), tags$ul(class = "mb-1 mt-1",
          row(csr_labA(), sc$A), row(csr_labB(), sc$B)),
        if (!same) tags$div(class = "small",
          tags$b("The two studies are on different scales."),
          " A relative measure is computed within each study before they are compared, so the",
          " comparison is still defined — but check that the measure below is appropriate for both."))
  })

  output$csr_measure_help <- renderUI({
    sc <- csr_scales(); req(sc, input$csr_measure)
    m <- input$csr_measure
    aA <- cs_measure_advice(sc$A)[[m]]; aB <- cs_measure_advice(sc$B)[[m]]
    cls <- function(st) switch(st, recommended = "text-success", usable = "text-muted",
                               discouraged = "text-warning", "text-danger")
    div(class = "small mb-2",
        div(class = cls(aA$status), tags$b(paste0(csr_labA(), ": ", aA$status, ". ")), aA$why),
        if (!identical(aA$why, aB$why))
          div(class = cls(aB$status), tags$b(paste0(csr_labB(), ": ", aB$status, ". ")), aB$why))
  })

  output$csr_method_help <- renderUI({
    ad <- cs_method_advice()[[input$csr_method %||% "pearson"]]
    div(class = "small text-muted mb-2", tags$b(paste0(ad$status, ". ")), ad$why)
  })

  csr_res <- eventReactive(input$csr_go, {
    req(input$csr_base)
    withProgress(message = "Comparing trajectories", value = .3, {
      cs_trajectory(csr_dA(), csr_dB(), input$csr_base, input$csr_measure,
        manual = if (identical(input$csr_align, "manual")) xs_manual_map() else NULL,
                    isTRUE(input$csr_1to1), input$csr_method, csr_labA(), csr_labB(),
                    align = input$csr_align %||% "label")
    })
  }, ignoreNULL = FALSE)

  observeEvent(csr_res(), {
    o <- csr_res(); req(isTRUE(o$ok))
    updateSelectizeInput(session, "csr_gene", choices = o$per_gene$feature,
                         selected = o$per_gene$feature[1], server = TRUE)
  })

  csr_ceiling <- eventReactive(input$csr_go, {
    if (!isTRUE(input$csr_ceiling)) return(NULL)
    withProgress(message = "Split-half noise ceiling", value = .3, {
      cs_noise_ceiling(csr_dA(), csr_dB(), input$csr_base, input$csr_measure,
        manual = if (identical(input$csr_align, "manual")) xs_manual_map() else NULL,
                       isTRUE(input$csr_1to1), csr_labA(), csr_labB(),
                       align = input$csr_align %||% "label")
    })
  }, ignoreNULL = FALSE)

  csr_levels <- reactive({ c("same direction", "opposite") })

  csr_export <- reactive({ o <- csr_res(); if (!isTRUE(o$ok)) return(NULL)
    list(summary = o$long, ylab = CS_MEASURE_LAB[[o$measure]], ttl = o$ttl,
         splitby = NULL, csv = o$long,
         gg = function(s, yl, tt, co) cs_traj_scatter_gg(csr_res(), co),
         prefix = sprintf("trajectory_%s_vs_%s_%s", o$labA, o$labB, o$measure)) })
  csr_cols <- apprServer("csr_appr", csr_levels, csr_export)

  csr_hist_export <- reactive({ o <- csr_res(); req(isTRUE(o$ok))
    list(summary = o$per_gene, ylab = "features", ttl = o$ttl, splitby = NULL,
         csv = o$per_gene, prefix = "trajectory_correlation_distribution",
         gg = function(s, yl, tt, co) cs_traj_hist_gg(csr_res())) })
  csr_hist_cols <- apprServer("csr_hist_appr", reactive(character(0)), csr_hist_export)

  csr_gene_export <- reactive({ o <- csr_res(); req(isTRUE(o$ok), input$csr_gene)
    list(summary = o$per_gene[o$per_gene$feature == input$csr_gene, , drop = FALSE],
         ylab = CS_MEASURE_LAB[[o$measure]], ttl = input$csr_gene, splitby = NULL,
         prefix = paste0("trajectory_", gsub("[^A-Za-z0-9]+", "_", input$csr_gene)),
         gg = function(s, yl, tt, co) cs_traj_gene_gg(csr_res(), isolate(input$csr_gene), co)) })
  csr_gene_cols <- apprServer("csr_gene_appr",
    reactive({ o <- csr_res(); if (!isTRUE(o$ok)) character(0) else c(o$labA, o$labB) }),
    csr_gene_export)

  output$csr_scatter   <- renderPlotly({ o <- csr_res(); req(o); cs_traj_scatter(o, csr_cols()) })
  output$csr_hist      <- renderPlotly({ o <- csr_res(); req(o); cs_traj_hist(o) })
  output$csr_gene_plot <- renderPlotly({ o <- csr_res(); req(o, isTRUE(o$ok), input$csr_gene)
                                         cs_traj_gene(o, input$csr_gene) })
  output$csr_series    <- renderPlotly({ o <- csr_res(); req(o, isTRUE(o$ok))
    plot_ly(o$by_tp, x = ~timepoint, y = ~r, type = "scatter", mode = "lines+markers",
            line = list(color = PIG, width = 2.5), marker = list(color = PIG, size = 9),
            text = ~sprintf("r = %.3f · R2 = %.3f<br>%s%% same direction", r, r2, same_direction_pct),
            hoverinfo = "text") %>%
      ply(xaxis = list(title = paste0("timepoint (relative to ", o$baseline, ")")),
          yaxis = list(title = "cross-study r", range = c(-1, 1), zeroline = TRUE),
          title = list(text = "Agreement at each timepoint", x = .02, font = list(size = 12.5)))
  })

  output$csr_ceiling_note <- renderUI({
    ce <- csr_ceiling(); o <- csr_res(); req(isTRUE(o$ok))
    base <- sprintf("Pooled r over all %s feature x timepoint values: %.3f. Median per-feature trajectory r: %s. %s%% of features peak at the same timepoint in both studies.%s",
                    format(nrow(o$long), big.mark = ","), o$r_pooled,
                    if (is.na(o$median_trajectory_r)) "not available (needs at least three timepoints)" else sprintf("%.3f", o$median_trajectory_r),
                    o$pct_same_peak,

                    if ((o$n_dropped %||% 0) > 0)
                      sprintf(" %s of %s shared features were dropped for having a zero or missing value at the baseline or at some later timepoint.",
                              format(o$n_dropped, big.mark = ","),
                              format(o$n_dropped + o$n_features, big.mark = ",")) else "")

    pn <- if (!is.null(o$pairing_note))
      div(class = "alert alert-warning py-2 px-3 small mb-2",
          tags$b("Unpaired baseline. "), o$pairing_note) else NULL
    if (is.null(ce)) return(tagList(pn, div(class = "text-muted small mb-2", base)))
    tagList(pn,
      div(class = "text-muted small mb-2", base, tags$br(),
          tags$b("Noise ceiling (split-half within each study): "),
          paste(sprintf("%s median r = %.3f", ce$study, ce$median_r), collapse = "; "),
          " — a between-study correlation cannot exceed the reproducibility of the studies it is computed from."))
  })

  output$csr_tab <- renderDT({
    o <- csr_res(); req(isTRUE(o$ok))
    datatable(o$per_gene, rownames = FALSE, selection = "none",
              options = list(pageLength = 10, scrollX = TRUE, dom = "ftip"),
              caption = htmltools::tags$caption(
                style = "caption-side: top; text-align: left; padding: 0 0 8px 0; color:#555;",
                sprintf("%s · matched by %s · %s", o$ttl, o$rule, CS_MEASURE_LAB[[o$measure]])))
  })
  output$csr_dl <- downloadHandler(
    filename = function() { o <- csr_res(); sprintf("trajectory_perfeature_%s_vs_%s.csv", o$labA, o$labB) },
    content = function(f) { o <- csr_res(); req(isTRUE(o$ok)); write.csv(o$per_gene, f, row.names = FALSE) })
  output$csr_dl_long <- downloadHandler(
    filename = function() { o <- csr_res(); sprintf("trajectory_long_%s_vs_%s.csv", o$labA, o$labB) },
    content = function(f) { o <- csr_res(); req(isTRUE(o$ok)); write.csv(o$long, f, row.names = FALSE) })

  dl_mat <- function(key) downloadHandler(
    filename = function() paste0(key, "_expression.csv"),
    content = function(f) write.csv(DS[[key]]$expr, f))
  dl_meta <- function(key) downloadHandler(
    filename = function() paste0(key, "_metadata.csv"),
    content = function(f) write.csv(DS[[key]]$meta, f))
  output$dl_pig_mrna <- dl_mat("pig_mrna"); output$dl_pig_pro <- dl_mat("pig_protein")
  output$dl_hum_mrna <- dl_mat("hum_mrna"); output$dl_hum_pro <- dl_mat("hum_protein")
  output$dl_meta_pig <- dl_meta("pig_mrna"); output$dl_meta_hum <- dl_meta("hum_mrna")

  .per_sample <- function(d, values, value_name = "value") {
    samp <- intersect(names(values), rownames(d$meta))
    if (!length(samp)) return(NULL)
    out <- data.frame(sample = samp, stringsAsFactors = FALSE)
    out[[value_name]] <- as.numeric(values[samp])
    cbind(out, d$meta[samp, , drop = FALSE])
  }

  byod_meta_raw <- reactiveVal(NULL)   # raw uploaded metadata (for the config UI)
  byod_study    <- reactiveVal(NULL)   # the ACTIVE study (layers + shared metadata)
  byod_msg      <- reactiveVal(NULL)   # status alert

  session$userData$bundle_ctx <- function(panel) {
    i <- match(panel, EXPORT_PANEL_MAP$panel)
    if (is.na(i)) return(NULL)
    if (bundle_is_pair(EXPORT_PANEL_MAP$module[i])) {

      pr <- tryCatch(bundle_pair_folder(csr_idA(), csr_idB()), error = function(e) NULL)
      if (is.null(pr)) return(NULL)
      return(list(dataset = BUNDLE_PAIR_DATASET, layer = pr))
    }
    s <- tryCatch(byod_study(), error = function(e) NULL)
    if (is.null(s)) return(NULL)
    ln <- isolate(input$byod_layer) %||% s$layer_names[1]
    if (!ln %in% s$layer_names) ln <- s$layer_names[1]
    list(dataset = bundle_study_id(s), layer = ln)
  }

  BYOD_SLOT_KEYS <- paste0("ds", seq_len(6L))
  byod_ws      <- reactiveVal(setNames(vector("list", length(BYOD_SLOT_KEYS)), BYOD_SLOT_KEYS))
  byod_nslots  <- reactiveVal(2L)      # how many slots the interface offers
  byod_aslot   <- reactiveVal(NULL)    # active slot key; NULL when a built-in is active

  xs_slot_label <- function(k) {
    i <- match(k, BYOD_SLOT_KEYS)
    if (is.na(i)) k else paste("Dataset", LETTERS[i])
  }
  xs_slot_keys <- reactive(BYOD_SLOT_KEYS[seq_len(byod_nslots())])
  xs_ws_filled <- reactive({ ws <- byod_ws(); Filter(Negate(is.null), ws[xs_slot_keys()]) })
  output$xs_ws_n <- renderText(as.character(length(xs_ws_filled())))
  outputOptions(output, "xs_ws_n", suspendWhenHidden = FALSE)

  xs_ws_bytes <- reactive({
    f <- xs_ws_filled()
    if (!length(f)) return(0)
    sum(vapply(f, function(r) as.numeric(utils::object.size(r$study)), numeric(1)))
  })

  observe({
    sl <- xs_slot_keys(); ws <- byod_ws()
    ch <- setNames(sl, vapply(sl, function(k)
      sprintf("%s %s", xs_slot_label(k),
              if (is.null(ws[[k]])) "\u00b7 empty" else paste0("\u00b7 ", ws[[k]]$study$name)),
      character(1)))
    cur <- isolate(input$byod_target)
    updateRadioButtons(session, "byod_target", choices = ch,
                       selected = if (!is.null(cur) && cur %in% sl) cur else sl[1])
  })

  observeEvent(input$byod_addslot, {
    n <- byod_nslots()
    if (n >= length(BYOD_SLOT_KEYS))
      return(byod_msg(list(type = "warning",
        txt = sprintf("%d dataset slots is the maximum; every one is held in memory for the session.",
                      length(BYOD_SLOT_KEYS)))))
    byod_nslots(n + 1L)
  })

  xs_activate <- function(k) {
    ws <- byod_ws(); r <- ws[[k]]
    if (is.null(r)) return(invisible(FALSE))
    byod_aslot(k); byod_study(r$study)
    byod_msg(list(type = "success",
      txt = sprintf("%s (%s) is now the active study; every module on Analyse runs on it.",
                    r$study$name, xs_slot_label(k))))
    invisible(TRUE)
  }
  xs_drop <- function(k) {
    ws <- byod_ws(); if (is.null(ws[[k]])) return(invisible(FALSE))
    nm <- ws[[k]]$study$name
    ws[k] <- list(NULL); byod_ws(ws)
    if (identical(isolate(byod_aslot()), k)) {
      left <- Filter(Negate(is.null), ws)
      byod_aslot(NULL)
      if (length(left)) xs_activate(names(left)[1]) else {
        byod_study(NULL)
        byod_msg(list(type = "warning",
          txt = sprintf("Removed %s. No study is active; load a built-in cohort or build another.", nm)))
      }
    } else byod_msg(list(type = "success", txt = sprintf("Removed %s.", nm)))
    invisible(TRUE)
  }

  for (k_ in BYOD_SLOT_KEYS) local({
    k <- k_
    observeEvent(input[[paste0("xs_use_", k)]], xs_activate(k), ignoreInit = TRUE)
    observeEvent(input[[paste0("xs_drop_", k)]], xs_drop(k), ignoreInit = TRUE)
  })

  output$xs_workspace <- renderUI({
    sl <- xs_slot_keys(); ws <- byod_ws(); act <- byod_aslot()
    body <- lapply(sl, function(k) {
      r <- ws[[k]]
      if (is.null(r))
        return(div(class = "border rounded p-2 mb-2",
          style = "border-style:dashed !important;",
          tags$b(class = "text-muted", xs_slot_label(k)),
          div(class = "small text-muted",
              "Empty. Pick this slot above, then upload matrices and declare its design.")))
      st <- r$study
      u  <- resolved_unit_col(st$layers[[1]])
      nu <- if (is.null(u)) NA_integer_ else length(unique(as.character(st$layers[[1]]$meta[[u]])))
      div(class = paste("border rounded p-2 mb-2", if (identical(act, k)) "border-primary" else ""),
        div(class = "d-flex align-items-center gap-2",
            tags$b(xs_slot_label(k)), tags$span(class = "pill", st$name),
            if (identical(act, k)) tags$span(class = "badge bg-primary", "active")),
        tags$ul(class = "small text-muted mb-1 mt-1",
          tags$li(sprintf("%d layer(s): %s", length(st$layer_names),
                          paste(st$layer_names, collapse = ", "))),
          tags$li(sprintf("%d timepoints: %s", length(st$tp_levels),
                          paste(st$tp_levels, collapse = ", "))),
          tags$li(if (is.null(u)) "no subject declared - the profile is the unit"
                  else sprintf("%d independent units (%s)%s", nu, u,
                               if (!is.null(st$nest)) sprintf(", nested in %s", st$nest) else "")),
          tags$li(sprintf("%s features in %s", format(nrow(st$layers[[1]]$expr), big.mark = ","),
                          st$layer_names[1]))),
        div(class = "d-flex gap-2",
          actionButton(paste0("xs_use_", k), "Analyse this", class = "btn-sm btn-outline-primary"),
          actionButton(paste0("xs_drop_", k), "Remove", class = "btn-sm btn-outline-secondary")))
    })
    tagList(body,
      div(class = "d-flex gap-2 align-items-center",
        actionButton("byod_addslot", "Add another dataset slot", class = "btn-sm btn-outline-secondary"),
        tags$span(class = "small text-muted",
          sprintf("%s held in memory", format(structure(xs_ws_bytes(), class = "object_size"),
                                              units = "auto")))),

      local({
        host <- byod_host_memory_bytes()
        used <- xs_ws_bytes()
        if (!is.finite(host) || host <= 0 || used < 0.2 * host) return(NULL)
        div(class = "alert alert-warning py-2 px-3 small mt-2",
            sprintf("The workspace is holding %s, about %.0f%% of the %s this machine reports. Remove a dataset you are finished with before building another, or run TRACE-Omics on a larger machine.",
                    format(structure(used, class = "object_size"), units = "auto"),
                    100 * used / host,
                    format(structure(host, class = "object_size"), units = "auto")))
      }),
      if (length(xs_ws_filled()) >= 2)
        div(class = "alert alert-success py-2 px-3 small mt-2",
            tags$b("Two datasets are loaded. "),
            "Cross-study comparison can now run the same contrast separately in each and compare the answers.",
            div(class = "mt-1", actionButton("go_compare_upload2", "Compare them",
                                             class = "btn-sm btn-primary"))))
  })

  output$byod_ready <- renderText(if (is.null(byod_study())) "no" else "yes")
  outputOptions(output, "byod_ready", suspendWhenHidden = FALSE)

  byod_layer_d <- reactive({
    s <- byod_study(); req(s)
    ln <- input$byod_layer %||% s$layer_names[1]; if (!ln %in% s$layer_names) ln <- s$layer_names[1]
    s$layers[[ln]]
  })
  byod_subj <- reactive({ s <- byod_study(); if (is.null(s)) NULL else s$subject })

  observeEvent(input$byod_meta, {
    tryCatch({ byod_meta_raw(byod_read(input$byod_meta$datapath)); byod_msg(NULL) },
      error = function(e) byod_msg(list(type = "danger", txt = paste("Metadata read error:", conditionMessage(e)))))
  })

  BYOD_NLAYERS <- 3L
  byod_layer_file <- function(i) input[[paste0("byod_expr", i)]]

  for (i in seq_len(BYOD_NLAYERS)) local({
    j <- i
    output[[paste0("byod_layer_shown", j)]] <- renderText(
      if (is.null(byod_layer_file(j))) "no" else "yes")
    outputOptions(output, paste0("byod_layer_shown", j), suspendWhenHidden = FALSE)

    observeEvent(input[[paste0("byod_type", j)]], {
      p <- BYOD_OMICS_TYPES[[input[[paste0("byod_type", j)]]]]; req(p)
      updateTextInput(session, paste0("byod_name", j), value = p$name)
      updateTextInput(session, paste0("byod_unit", j), value = p$unit)
      updateSelectInput(session, paste0("byod_norm", j), selected = p$norm)
      updateSelectInput(session, paste0("byod_annotation", j), selected = if (isTRUE(p$genesets)) "human" else "none")
    }, ignoreInit = TRUE)

    output[[paste0("byod_mapping_hint", j)]] <- renderUI({
      helpText(annotation_resource_hint(
        input$byod_organism %||% "Human",
        input[[paste0("byod_idtype", j)]] %||% "symbol",
        input[[paste0("byod_idspace", j)]] %||% "native",
        input[[paste0("byod_annotation", j)]] %||% "human",
        input[[paste0("byod_mapmode", j)]] %||% "builtin"))
    })

    output[[paste0("byod_layer_hint", j)]] <- renderUI({
      ty <- input[[paste0("byod_type", j)]] %||% "other"
      p <- BYOD_OMICS_TYPES[[ty]]
      f <- byod_layer_file(j)
      if (is.null(f)) return(helpText(p$hint))
      rn <- tryCatch(as.character(byod_read(f$datapath)[[1]]), error = function(e) NULL)
      if (is.null(rn) || !length(rn)) return(helpText(p$hint))
      info <- byod_detect_ids(rn, ty)
      div(class = "small text-muted mt-1", tags$b("Detected: "), info$label,
        ". Check the ID type and source organism. Mapping coverage is calculated when the study is built.")
    })
  })

  byod_clean_opts <- reactive(list(
    trim_names    = isTRUE(input$clean_trim),
    strip_version = isTRUE(input$clean_strip),
    upper_names   = isTRUE(input$clean_upper),
    dup_rule      = input$clean_dup %||% "mean",
    min_value     = input$clean_minval %||% 0,
    min_samples   = input$clean_minsamp %||% 0,
    max_missing   = input$clean_maxna %||% 100,
    drop_constant = isTRUE(input$clean_const),
    drop_failed   = isTRUE(input$clean_failed),
    failed_frac   = input$clean_failedfrac %||% 0.01))

  output$byod_clean_preview <- renderUI({
    f <- byod_layer_file(1); if (is.null(f)) return(NULL)
    x <- tryCatch(byod_read(f$datapath), error = function(e) NULL); req(x)
    M <- suppressWarnings(vapply(x[, -1, drop = FALSE],
           function(col) as.numeric(as.character(col)), numeric(nrow(x))))
    if (!is.matrix(M)) M <- matrix(M, nrow = nrow(x))
    rownames(M) <- as.character(x[[1]]); colnames(M) <- names(x)[-1]
    cl <- tryCatch(byod_clean_matrix(M, byod_clean_opts(),
                     is_counts = identical(byod_diagnose(M)$scale, "raw counts")),
                   error = function(e) NULL)
    if (is.null(cl)) return(NULL)
    if (!cl$any_applied)
      return(div(class = "small text-muted", sprintf(
        "Layer 1 as uploaded: %s features x %d profiles. No cleaning step is active.",
        format(nrow(M), big.mark = ","), ncol(M))))
    kept <- nrow(cl$M); lost <- nrow(M) - kept
    div(class = if (kept < 0.25 * nrow(M)) "alert alert-warning py-2 px-3 small"
                else "alert alert-secondary py-2 px-3 small",
        tags$b(sprintf("Layer 1 preview: %s of %s features kept",
                       format(kept, big.mark = ","), format(nrow(M), big.mark = ","))),
        if (lost > 0) sprintf(" (%s removed)", format(lost, big.mark = ",")),
        tags$ul(class = "mb-0 mt-1", lapply(byod_clean_summary(cl), tags$li)))
  })

  output$byod_config <- renderUI({
    m <- byod_meta_raw()
    if (is.null(m)) return(helpText("Upload a metadata file to pick the time and pairing columns."))
    cols <- names(m)[-1]
    tagList(
      selectInput("byod_tpcol", "Time / order column", choices = cols,
                  selected = if ("Timepoint" %in% cols) "Timepoint" else cols[1]),
      selectInput("byod_subj", "Subject / pairing column (optional)", choices = c("(none)", cols), selected = "(none)"),
      selectInput("byod_nest", "Higher-level unit the subject nests inside (optional)",
                  choices = c("(none)", cols), selected = "(none)"),
      helpText(paste("If each subject belongs to a larger unit - two kidneys from one animal,",
                     "several biopsies from one donor - declare it here and contrasts are tested",
                     "across those units, not across subjects. Nothing is inferred: leave it blank",
                     "and the subject stays the unit.")),
      selectInput("byod_batch", "Batch column — ComBat correction (optional)", choices = c("(none)", cols), selected = "(none)"),

      textInput("byod_studyname", "Study name", value = "My study"))
  })

  byod_store_study <- function(s, slot = NULL, source = "upload") {
      k  <- slot %||% isolate(input$byod_target) %||% BYOD_SLOT_KEYS[1]
      s$dataset_id <- bundle_new_study_id(session, k)
      ws <- isolate(byod_ws())
      had <- !is.null(ws[[k]])
      ws[[k]] <- list(key = k, study = s, built = Sys.time())
      byod_ws(ws); byod_aslot(k); byod_study(s)

      log_event("load-study", source = source, name = s$name, slot = xs_slot_label(k), engine = "shiny",
                reference = s$reference_key %||% "",
                replaced = had, layers = paste(s$layer_names, collapse = ", "),
                timepoints = length(s$tp_levels),
                subject = s$subject %||% "(none declared)",
                nest = s$nest %||% "(none declared)",
                organism = s$organism %||% "(not declared)")
      byod_msg(list(type = "success",
        txt = sprintf("%s %s: %s - %d layer(s) [%s], %d timepoints%s. It is the active study.",
          if (had) "Replaced" else "Built into", xs_slot_label(k),
          s$name, length(s$layer_names), paste(s$layer_names, collapse = ", "), length(s$tp_levels),
          if (!is.null(s$subject))
            paste0(" \u00b7 ", byod_design_text(s$layers[[1]])) else "")))
  }

  byod_do_build <- function(exprs, meta, tp_col, units, name, tp_order = NULL, subject_col = NULL, nest_col = NULL,
                            norms = NULL, batch_col = NULL, id_maps = NULL,
                            cleans = NULL, omics = NULL, slot = NULL, organism = NULL, annotations = NULL) {
    withProgress(message = "Building study", detail = "clean + normalize + batch + PCA + UMAP…", value = .5, {
      s <- build_byod_study(exprs, meta, tp_col = tp_col, units = units, name = name, tp_order = tp_order,
                            subject_col = subject_col, norms = norms, batch_col = batch_col,
                            id_maps = id_maps, cleans = cleans, omics = omics, nest_col = nest_col,
                            organism = organism, annotations = annotations)
      byod_store_study(s, slot)
    })
  }

  observeEvent(input$byod_build, {
    req(input$byod_expr1, input$byod_meta)
    tryCatch({
    exprs <- list(); units <- list(); norms <- list(); idm <- list(); anns <- list()
    cleans <- list(); omics <- list(); cl <- byod_clean_opts()
    for (i in seq_len(BYOD_NLAYERS)) {
      f <- byod_layer_file(i); if (is.null(f)) next
      ty <- input[[paste0("byod_type", i)]] %||% "other"
      nm <- input[[paste0("byod_name", i)]]
      if (is.null(nm) || !nzchar(trimws(nm))) nm <- BYOD_OMICS_TYPES[[ty]]$name %||% paste0("Layer", i)

      while (nm %in% names(exprs)) nm <- paste0(nm, "_", i)
      exprs[[nm]]  <- byod_read(f$datapath)
      units[[nm]]  <- input[[paste0("byod_unit", i)]] %||% "expression"
      norms[[nm]]  <- input[[paste0("byod_norm", i)]] %||% "Auto-detect"
      idm[[nm]]    <- input[[paste0("byod_idmap", i)]] %||% "none"
      cleans[[nm]] <- cl
      omics[[nm]]  <- ty
      mf <- input[[paste0("byod_mapping", i)]]; gf <- input[[paste0("byod_gmt", i)]]
      anns[[nm]] <- list(mode = input[[paste0("byod_annotation", i)]] %||%
        if (input$byod_organism %in% names(ANNOTATION_SPECIES)) "human" else "none",
        id_type = input[[paste0("byod_idtype", i)]] %||% "symbol",
        id_species = input[[paste0("byod_idspace", i)]] %||% "native")
      choice <- input[[paste0("byod_mapmode", i)]] %||% "builtin"
      if (!choice %in% c("builtin", "custom")) stop("Unknown mapping source.")
      if (anns[[nm]]$mode == "human" && choice == "custom") {
        if (is.null(mf)) stop(sprintf("Layer '%s': upload the selected custom mapping or switch to built-in mapping.", nm))
        anns[[nm]]$mapping <- annotation_custom_map(mf$datapath, input[[paste0("byod_mapsource", i)]])
      }
      if (isTRUE(input[[paste0("byod_usegmt", i)]])) {
        if (is.null(gf)) stop(sprintf("Layer '%s': upload a GMT or uncheck Use a custom GMT.", nm))
        anns[[nm]]$gmt <- annotation_read_gmt(gf$datapath)
      }
    }
    validate(need(length(exprs) >= 1, "Upload at least one matrix."))
    subj <- input$byod_subj; if (is.null(subj) || subj == "(none)") subj <- NULL
    nest <- input$byod_nest; if (is.null(nest) || nest == "(none)") nest <- NULL
    bat <- input$byod_batch; if (is.null(bat) || bat == "(none)") bat <- NULL
    org <- input$byod_organism
    if (identical(org, "Other")) {
      org <- trimws(input$byod_custom_organism %||% "")
      if (!nzchar(org)) stop("Enter the scientific name for Other organism.")
    }
    if (is.null(org) || org == "(not declared)") org <- NULL
      byod_do_build(exprs, byod_read(input$byod_meta$datapath), input$byod_tpcol, units,
                    input$byod_studyname %||% "My study", subject_col = subj, norms = norms,
                    batch_col = bat, id_maps = idm, cleans = cleans, omics = omics, nest_col = nest,
                    organism = org, annotations = anns)
    },
      error = function(e) byod_msg(list(type = "danger", txt = sprintf(
        "Build failed, and nothing was changed - %s is still the active study. %s",
        isolate(byod_study()$name) %||% "no study", conditionMessage(e)))))
  })

  go <- function(target) nav_select("main_nav", target, session = session)
  observeEvent(input$go_atlas,           go("Kidney Perfusion Reference"))
  observeEvent(input$go_atlas2,          go("Kidney Perfusion Reference"))
  observeEvent(input$go_upload2,         go("Your data"))
  observeEvent(input$go_upload,          go("Your data"))
  observeEvent(input$go_analyse_atlas,   go("Analyse"))
  observeEvent(input$go_analyse_upload,  go("Analyse"))
  observeEvent(input$go_compare_atlas,   go("Cross-study comparison"))
  observeEvent(input$go_compare_upload,  go("Cross-study comparison"))
  observeEvent(input$csr_goto_upload,    go("Your data"))

  output$atlas_contents <- renderUI({
    rows <- lapply(names(STUDIES), function(k) {
      st <- STUDIES[[k]]

      lay <- lapply(st$layer_names, function(ln) {
        d <- st$layers[[ln]]
        tags$li(sprintf("%s: %s features (%s)", ln,
                        format(nrow(d$expr), big.mark = ","), d$unit))
      })

      d1 <- st$layers[[1]]; m1 <- d1$meta
      design <- list()
      if (!is.null(m1$Pig_Number) && !is.null(m1$WI_status)) {
        n_an <- length(unique(m1$Pig_Number)); n_kid <- length(unique(m1$Kidney_ID))
        target <- n_kid * length(st$tp_levels)
        design <- list(
          tags$li(sprintf("%d animals, %d kidneys — one warm-ischaemia kidney and one control kidney per animal",
                          n_an, n_kid)),
          tags$li("warm ischaemia is a within-animal contrast, and the first timepoint is an in vivo biopsy taken before retrieval"),
          tags$li(sprintf("%d biopsies by design (%d kidneys x %d timepoints)",
                          target, n_kid, length(st$tp_levels))))
      } else {
        nclin <- length(d1$clinical_cat %||% character(0)) + length(d1$clinical_num %||% character(0))
        if (nclin > 0)
          design <- list(tags$li(sprintf("%d donor and clinical variables linked to each kidney", nclin)))
      }
      tags$div(class = "mb-2",
        tags$b(st$name),
        tags$ul(class = "mb-1", lay,
          tags$li(sprintf("%d timepoints: %s", length(st$tp_levels),
                          paste(st$tp_levels, collapse = ", "))),
          design,
          if (!is.null(st$subject)) tags$li(paste("repeated sampling keyed by", st$subject))))
    })
    tagList(rows,
      div(class = "small text-muted",
          "The design is already declared for these cohorts, so no configuration is needed. ",
          "Load one and every module is immediately available."))
  })

  output$atlas_status <- renderUI({
    m <- byod_msg()
    if (is.null(m)) return(NULL)
    div(class = paste0("alert alert-", m$type, " py-2 px-3 mt-2 mb-0 small"), m$txt)
  })

  is_uploaded <- reactive({
    st <- byod_study()
    !is.null(st) && !bundle_study_id(st) %in% vapply(STUDIES, bundle_study_id, character(1))
  })
  output$upload_ready <- renderText(if (isTRUE(is_uploaded())) "yes" else "no")
  outputOptions(output, "upload_ready", suspendWhenHidden = FALSE)

  output$upload_guide <- renderUI({
    if (isTRUE(is_uploaded())) return(NULL)
    card(card_header("What the engine needs from you"),
      tags$ol(
        tags$li(tags$b("The matrices."), " Features in the first column, one column per sample. ",
                "Up to three layers, of any assay type; the second and third are optional. Give layers ",
                "the same sample IDs if they were measured on the same aliquot, and different IDs if ",
                "they were not — the engine checks, and withdraws sample-level cross-layer methods when ",
                "there is no overlap."),
        tags$li(tags$b("The unit of each layer."), " log2, log1p, linear counts, intensities, beta ",
                "values. This is not cosmetic: it decides which relative measures are meaningful, it is ",
                "the only source for a +1 offset, and it is what keeps a bounded 0–1 value from being ",
                "treated as a logarithm."),
        tags$li(tags$b("Species and identifiers."), " Declare the biological source and each layer's ID type. ",
                "Human gene sets use the layer's annotation; a custom GMT can instead use native feature IDs. ",
                "Expression rows retain their original identity after the chosen cleaning steps."),
        tags$li(tags$b("The design."), " Which column is time and in what order; which column identifies ",
                "the repeatedly sampled subject; whether subjects nest inside a higher-level unit such as ",
                "an animal; and which column, if any, is the treatment group. Everything else about your ",
                "metadata is derived from the data rather than declared \u2014 in particular whether a column ",
                "is a property of the subject or is measured again at every visit \u2014 and the build summary ",
                "lists what was derived, so no column disappears from a module without a stated reason."),
        tags$li(tags$b("Nothing else."), " The engine never guesses the design from column names or ",
                "sample order, and it refuses an analysis the declared design cannot support rather than ",
                "silently fitting a weaker one.")),
      div(class = "small text-muted",
          "Try the human kidney reference in the sidebar: mRNA alone or matched mRNA + protein, ",
          "with the same donor metadata and five timepoints. Both use the existing curated preprocessing."))
  })

  output$upload_built_summary <- renderUI({
    req(is_uploaded()); st <- byod_study(); req(st)
    lay <- lapply(st$layer_names, function(ln) {
      d <- st$layers[[ln]]; si <- cs_scale_info(d); ids <- d$ids
      tags$li(sprintf("%s: %s features x %d profiles — %s (declared '%s', %s)",
                      ln, format(nrow(d$expr), big.mark = ","), ncol(d$expr),
                      si$label, si$declared, si$source),

        if (!is.null(ids)) tags$div(class = "small text-muted",
          sprintf("identifiers: %s; %s of %s in Hallmark%s", ids$label,
                  format(ids$coverage_after$hit, big.mark = ","),
                  format(ids$coverage_after$n, big.mark = ","),
                  if (isTRUE(ids$map$applied)) paste0(" — ", ids$map$msg) else "")),

        if (!is.null(d$clean) && isTRUE(d$clean$any_applied))
          tags$div(class = "small text-muted",
                   "cleaning: ", paste(byod_clean_summary(d$clean), collapse = "; ")),
        tags$div(class = "small text-muted", annotation_summary(d)),
        if (!isTRUE(d$genesets_apply))
          tags$div(class = "small text-muted fst-italic",
                   "No configured gene set currently has at least three eligible features. Check the annotation settings or supply a GMT."))
    })
    tagList(
      downloadButton("byod_annotation_dl", "Download input data & annotation"),
      tags$ul(lay,
        tags$li(sprintf("%d timepoints: %s", length(st$tp_levels), paste(st$tp_levels, collapse = ", "))),
        if (!is.null(st$subject)) tags$li(paste("repeated sampling keyed by", st$subject))
             else tags$li("no subject column declared: samples are treated as independent"),
        if (length(st$layer_names) > 1)
          tags$li(sprintf("%d profiles and %s features shared between the layers",
                          length(st$matched_samples), format(length(st$shared_genes), big.mark = ",")))),

      {
        qcs <- Filter(Negate(is.null), lapply(st$layer_names, function(ln) {
          q <- st$layers[[ln]]$qc; if (is.null(q)) NULL else list(ln = ln, q = q) }))
        if (length(qcs)) tagList(lapply(qcs, function(z)
          div(class = "alert alert-warning py-2 px-3 small",
              tags$b(sprintf("%s — possible failed librar%s. ", z$ln,
                             if (length(z$q$samples) == 1) "y" else "ies")), z$q$msg)))
      },

      if (!is.null(st$excluded_untimed)) {
        eu <- st$excluded_untimed
        div(class = "alert alert-warning py-2 px-3 small",
            tags$b(sprintf("%d profile%s excluded.", eu$n, if (eu$n == 1) "" else "s")),
            sprintf(" Timepoint %s not among the levels you ordered, so %s cannot be placed in time. %s",
                    paste0("'", paste(eu$labels, collapse = "', '"), "'",
                           if (length(eu$labels) == 1) " is" else " are"),
                    if (eu$n == 1) "it" else "they",
                    "Add the level to the timepoint order if it belongs in the series."),
            tags$div(class = "text-muted mt-1",
                     paste("Excluded:", paste(head(eu$samples, 8), collapse = ", "),
                           if (eu$n > 8) sprintf("and %d more", eu$n - 8) else "")))
      },

      {
        sl  <- c(st$clinical_num %||% character(0), st$clinical_cat %||% character(0))
        vry <- st$subject_varying
        if (!length(sl) && !length(vry)) NULL else
        div(class = "mt-2 p-2 border rounded small",
          tags$b("How your other columns were classified"),
          tags$div(class = "text-muted mb-1",
            "You declared time, subject and nesting. Whether each remaining column is a ",
            "subject-level characteristic or a per-visit measurement is measured from the ",
            "data, by testing it for constancy within each subject."),
          tags$ul(class = "mb-1",
            tags$li(tags$b(sprintf("%d subject-level: ", length(sl))),
                    if (length(sl)) paste(sl, collapse = ", ") else "none",
                    tags$span(class = "text-muted",
                              " \u2014 offered to the subject-level modules, collapsed to one value per subject")),
            tags$li(tags$b(sprintf("%d change between visits: ", length(vry))),
                    if (length(vry)) paste(names(vry), collapse = ", ") else "none",
                    tags$span(class = "text-muted",
                              " \u2014 not offered there, because collapsing would pick one visit and call it the subject's value"))),
          if (length(vry)) tags$div(class = "text-muted",
            paste0("Disagreeing subjects: ",
                   paste(vapply(names(vry), function(k)
                     sprintf("%s (%d)", k, length(vry[[k]])), character(1)), collapse = "; "),
                   ". A per-visit column is still available everywhere that models time, ",
                   "and an association with one is tested within subject rather than across profiles.")))
      },
      div(class = "small text-muted",
          "This is now the same validated object as a built-in cohort."))
  })

  output$analyse_source_card <- renderUI({
    st <- byod_study()
    if (is.null(st)) return(div(class = "small text-muted", "No study loaded."))
    builtin <- bundle_study_id(st) %in% vapply(STUDIES, bundle_study_id, character(1))
    ws <- byod_ws(); act <- byod_aslot()
    others <- Filter(function(k) !is.null(ws[[k]]) && !identical(k, act), xs_slot_keys())
    tagList(
      tags$b("Active study"),
      div(class = "pill my-1", st$name),
      div(class = "small text-muted mb-2",
          if (builtin) "From the Kidney Perfusion Reference."
          else if (!is.null(st$reference_key)) "Workspace copy of the human kidney reference; not a new cohort."
          else "Your uploaded study."),
      if (length(others)) tagList(
        div(class = "small text-muted mb-1", "Run the same module on:"),
        div(class = "d-flex gap-1 flex-wrap mb-2", lapply(others, function(k)
          actionButton(paste0("xs_use_", k), sprintf("%s \u00b7 %s", xs_slot_label(k), ws[[k]]$study$name),
                       class = "btn-outline-primary btn-sm")))),
      div(class = "d-flex gap-2 flex-wrap",
          actionButton("go_atlas", "Switch cohort", class = "btn-outline-secondary btn-sm"),
          actionButton("go_upload", "Upload another", class = "btn-outline-secondary btn-sm")),
      hr(),
      tags$b("Active layer"),

      selectInput("byod_layer", NULL, choices = st$layer_names,
                  selected = isolate(input$byod_layer) %||% st$layer_names[1]))
  })
  output$byod_status2 <- renderUI({
    m <- byod_msg(); if (is.null(m)) return(NULL)
    div(class = paste0("alert alert-", m$type, " py-2 px-3 mb-0 small"), m$txt)
  })

  load_builtin_study <- function(key) {
    s <- STUDIES[[key]]
    if (is.null(s)) { byod_msg(list(type = "danger", txt = "Unknown study.")); return() }
    byod_meta_raw(data.frame(sample = rownames(s$meta), s$meta, check.names = FALSE))

    byod_aslot(NULL); byod_study(s)
    byod_msg(list(type = "success",

      txt = sprintf("Loaded %s: %d layer(s) [%s], %d timepoints%s.",
        s$name, length(s$layer_names), paste(s$layer_names, collapse = ", "),
        length(s$tp_levels),
        if (!is.null(s$subject)) paste0(", paired by ", s$subject) else "")))
  }
  observeEvent(input$builtin_load, {
    req(input$builtin_study)
    tryCatch(load_builtin_study(input$builtin_study),
      error = function(e) byod_msg(list(type = "danger", txt = conditionMessage(e))))
  })

  observe({
    if (is.null(isolate(byod_study()))) load_builtin_study(names(STUDY_LABELS)[1])
  })
  output$active_study_name <- renderText({ s <- byod_study(); if (is.null(s)) "—" else s$name })

  bundle_tick <- reactiveVal(0L)
  observe({
    analysis_log()
    if (identical(input$main_nav, "Reproducibility")) invalidateLater(4000, session)
    bundle_tick(isolate(bundle_tick()) + 1L)
  })

  bundle_inv <- reactive({ bundle_tick(); bundle_inventory(session) })

  output$bundle_summary <- renderUI({
    inv <- bundle_inv()
    st <- byod_study()
    if (!nrow(inv))
      return(div(class = "small text-muted",
                 "Nothing recorded yet. Every analysis you run is added here as it finishes."))

    ds <- if (is.null(st)) NULL else BUNDLE_SLUG(bundle_study_id(st))
    keep <- if (is.null(st)) rep(FALSE, nrow(inv)) else
      ifelse(inv$dataset == BUNDLE_PAIR_DATASET,
             vapply(inv$layer, bundle_pair_is_about, logical(1), dataset = bundle_study_id(st)),
             inv$dataset == ds)
    mine <- inv[keep, , drop = FALSE]; other <- inv[!keep, , drop = FALSE]
    recorded_layers <- bundle_recorded_layers(mine)
    tagList(
      div(class = "small text-muted mb-1",
          sprintf("%d analys%s recorded for %s, across %d layer%s (%s).",
                  length(unique(mine$analysis)), if (length(unique(mine$analysis)) == 1) "is" else "es",
                  st$name %||% "this study",
                  length(recorded_layers), if (length(recorded_layers) == 1) "" else "s",
                  format(structure(sum(mine$bytes), class = "object_size"), units = "auto"))),
      tags$ul(class = "small mb-1", lapply(sort(unique(mine$analysis)), function(a) {
        z <- mine[mine$analysis == a, , drop = FALSE]
        tags$li(tags$code(a), sprintf(" — %s (%d table%s, %d figure%s)",
                paste(z$layer, collapse = ", "), sum(z$n_tables),
                if (sum(z$n_tables) == 1) "" else "s", sum(z$n_figures),
                if (sum(z$n_figures) == 1) "" else "s")) })),
      if (nrow(other)) div(class = "small text-muted",
        sprintf("%d further analys%s are held for %s and will go into their own archive.",
                length(unique(other$analysis)), if (length(unique(other$analysis)) == 1) "is" else "es",
                paste(unique(other$dataset), collapse = ", "))))
  })

  output$byod_annotation_dl <- downloadHandler(
    filename = function() paste0("TRACE-Omics_input_annotation.", bundle_archive_format()),
    content = function(f) {
      st <- isolate(byod_study()); req(st)
      stage <- tempfile("trace_annotation_"); dir.create(stage)
      on.exit(unlink(stage, recursive = TRUE), add = TRUE)
      annotation_write_study(st, file.path(stage, "annotation"))
      bundle_write_archive(stage, f)
    })

  output$bundle_dl <- downloadHandler(
    filename = function() {
      st <- isolate(byod_study())

      sprintf("TRACE-Omics_%s_%s.%s", BUNDLE_SLUG(st$name %||% "session"),
              format(Sys.time(), "%Y%m%d_%H%M"), bundle_archive_format())
    },
    content = function(f) {
      st <- isolate(byod_study())
      withProgress(message = "Assembling the archive", detail = "rendering figures…", value = .4, {
        bundle_assemble(session, bundle_study_id(st) %||% "session", f,
                        log_records = isolate(analysis_log()),
                        design_text = tryCatch(byod_design_text(isolate(byod_layer_d())),
                                               error = function(e) NULL),

                        methods_text = { lr <- bundle_filter_log(isolate(analysis_log()),
                                                                 bundle_study_id(st) %||% "session")
                                         byod_methods_text(lr, unique(unlist(lapply(lr, function(r) r$params$study)))) },
                        env_text = capture.output(print(utils::sessionInfo())),
                        study_name = st$name, study = st,
                        related_studies = lapply(isolate(csr_registry()), `[[`, "study"),
                        include_comparison_inputs = isTRUE(isolate(input$bundle_comparison_inputs)))
      })
    })

  observeEvent(input$builtin_load, {
    s <- isolate(byod_study()); if (!is.null(s)) log_event("load-study", study = s$name, source = "built-in", engine = "shiny") })

  byod_load_human_example <- function(dual) {
    tryCatch({
      s <- byod_human_example(dual)
      byod_store_study(s, source = "reference-example")
    }, error = function(e) byod_msg(list(type = "danger", txt = conditionMessage(e))))
  }
  observeEvent(input$byod_example_human_rna, byod_load_human_example(FALSE))
  observeEvent(input$byod_example_human_multi, byod_load_human_example(TRUE))

  output$byod_status <- renderUI({
    m <- byod_msg(); if (is.null(m)) return(NULL)
    div(class = paste0("alert alert-", m$type, " mt-2 py-2"), style = "font-size:.85rem;", m$txt)
  })

  observeEvent(byod_study(), {
    s <- byod_study(); req(s)
    updateSelectInput(session, "byod_layer",  choices = s$layer_names, selected = s$layer_names[1])
    updateSelectInput(session, "byod_int_l1", choices = s$layer_names, selected = s$layer_names[1])
    updateSelectInput(session, "byod_int_l2", choices = s$layer_names,
                      selected = if (length(s$layer_names) >= 2) s$layer_names[2] else s$layer_names[1])
    updateSelectInput(session, "byod_conc_l1", choices = s$layer_names, selected = s$layer_names[1])
    updateSelectInput(session, "byod_conc_l2", choices = s$layer_names,
                      selected = if (length(s$layer_names) >= 2) s$layer_names[2] else s$layer_names[1])
    updateSelectInput(session, "byod_as_tp",  choices = c("All", s$tp_levels), selected = "All")
    updateSelectInput(session, "byod_diablo_out", choices = col_groups(s$layers[[1]]), selected = "Timepoint")

    clin_vars <- c(s$clinical_num, s$clinical_cat)
    if (length(clin_vars)) {
      cc_def <- clin_vars[which(clin_vars %in% c("Total duration CIT", "KDPI", "Donor age"))[1]]
      updateSelectInput(session, "byod_cc_var", choices = clin_vars, selected = cc_def %||% clin_vars[1])
      updateSelectInput(session, "byod_cc_tp",  choices = c("All (donor-averaged)" = "All", stats::setNames(s$tp_levels, s$tp_levels)), selected = "All")
    }
  })

  observeEvent(byod_layer_d(), {
    d <- byod_layer_d(); req(d)
    ac <- annotation_collections(d)
    for (id in c("byod_en_coll", "byod_pw_coll", "byod_cc_coll", "byod_gti_coll", "byod_tp_coll", "byod_de_ora_coll"))
      updateSelectInput(session, id, choices = ac, selected = ac[1])
    feats <- rownames(d$expr); dg <- if ("HAVCR1" %in% feats) "HAVCR1" else feats[1]
    updateSelectizeInput(session, "byod_gene",    choices = feats, selected = dg, server = TRUE)
    updateSelectizeInput(session, "byod_dr_gene", choices = feats, selected = dg, server = TRUE)
    updateSelectizeInput(session, "byod_as_feat", choices = feats, selected = dg, server = TRUE)
    gv <- group_vars(d); cg <- col_groups(d)
    updateSelectInput(session, "byod_dr_col", choices = c(gv, "Gene expression"), selected = "Timepoint")
    updateSelectInput(session, "byod_split",  choices = c("None", setdiff(cg, "Timepoint")), selected = "None")
    updateSelectInput(session, "byod_de_var", choices = gv, selected = "Timepoint")
    updateSelectInput(session, "byod_en_var", choices = gv, selected = "Timepoint")
    updateSelectInput(session, "byod_int_var", choices = gv, selected = "Timepoint")

    updateSelectInput(session, "byod_hm_by",
                      choices = c(setNames(byod_order_candidates(d), byod_order_candidates(d)),
                                  "cluster the samples" = "cluster"),
                      selected = if ("Timepoint" %in% cg) "Timepoint" else cg[1])
    updateSelectInput(session, "byod_ms_split", choices = c("None", setdiff(cg, "Timepoint")), selected = "None")
    updateSelectInput(session, "byod_pw_by",  choices = cg, selected = "Timepoint")
    updateSelectInput(session, "byod_pw_split", choices = c("None", setdiff(cg, "Timepoint")), selected = "None")
    updateSelectInput(session, "byod_as_var", choices = gv, selected = if (length(gv) > 1) gv[2] else gv[1])
    if (!nzchar(input$byod_hm_genes %||% ""))
      updateTextAreaInput(session, "byod_hm_genes", value = "HAVCR1, LCN2, VIM, SPP1, CCL2, IL6, FN1, TIMP1, CD3E, CD68")
  })

  output$byod_overview <- renderUI({
    s <- byod_study(); req(s)

    builtin <- bundle_study_id(s) %in% vapply(STUDIES, bundle_study_id, character(1))
    d1 <- s$layers[[1]]; m1 <- d1$meta
    design_head <- if (!builtin) NULL
      else if (!is.null(m1$Pig_Number) && !is.null(m1$Kidney_ID))
        sprintf("%d animals · %d kidneys · %d timepoints",
                length(unique(m1$Pig_Number)), length(unique(m1$Kidney_ID)), length(s$tp_levels))
      else if (!is.null(m1$Kidney_ID))
        sprintf("%d kidneys · %d timepoints",
                length(unique(m1$Kidney_ID)), length(s$tp_levels))
      else sprintf("%d timepoints", length(s$tp_levels))
    cards <- lapply(s$layer_names, function(ln) { d <- s$layers[[ln]]
      card(card_header(paste0(ln, "  ·  ", d$unit)),
        if (is.null(design_head)) tagList(
          div(class = "statnum", ncol(d$expr)), div(class = "statlab", "profiles"))
        else tagList(
          div(style = "font-size:1.05rem;font-weight:600;color:#1E3A5F;", design_head),
          div(class = "statlab", "study design")),
        div(style = "margin-top:6px;", tags$b(format(nrow(d$expr), big.mark = ",")), " features"),
        div(style = "margin-top:8px;font-size:.82rem;color:#556;",
          tags$b("Input scale: "), d$diag$scale, tags$br(),
          tags$b("Normalization: "), d$norm,
          sprintf("  (range %.1f–%.1f, %.0f%% zeros)", d$diag$min, d$diag$max, d$diag$pct_zero), tags$br(),
          tags$b("Batch correction: "), d$batch %||% "none")) })
    lc <- do.call(layout_columns, c(list(col_widths = if (length(cards) == 1) 12 else c(6, 6), fill = FALSE), cards))

    ln <- input$byod_layer %||% s$layer_names[1]
    if (!ln %in% s$layer_names) ln <- s$layer_names[1]
    da <- s$layers[[ln]]
    gv <- setdiff(group_vars(da), "Timepoint")
    tagList(lc,
      card(card_header(s$name),

        if (builtin) p(tags$b("Timepoints: "), paste(s$tp_levels, collapse = "  ·  "))
        else {
          tpc <- table(factor(da$meta$Timepoint, levels = s$tp_levels))
          p(tags$b(sprintf("Timepoints (%s): ", ln)),
            paste(names(tpc), as.integer(tpc), sep = " n=", collapse = "  ·  "))
        },
        p(tags$b("Pairing / subject: "), s$subject %||% "none"),
        p(tags$b("Extra grouping variables: "), if (length(gv)) paste(gv, collapse = ", ") else "none"),
        tags$details(tags$summary(style = "cursor:pointer;", "Preview sample metadata"),
          tags$pre(style = "font-size:.75rem;max-height:220px;overflow:auto;",
                   paste(capture.output(print(utils::head(da$meta, 8))), collapse = "\n")))))
  })

  output$byod_qc_plot <- renderPlotly({ d <- byod_layer_d(); req(d); byod_dist_plot(d) })

  byod_dr_base <- reactive({ d <- byod_layer_d(); req(d)
    byod_dimred(d, input$byod_dr_type %||% "UMAP", input$byod_dr_col %||% "Timepoint", input$byod_dr_gene) })
  output$byod_dr_kind <- renderText(byod_dr_base()$kind %||% "discrete")
  outputOptions(output, "byod_dr_kind", suspendWhenHidden = FALSE)
  byod_dr_levels <- reactive(byod_dr_base()$levels %||% character(0))
  byod_dr_export <- reactive({ o <- byod_dr_base()
    list(summary = o$summary, ylab = o$ax[2], ttl = o$ttl, splitby = o$scale_title,
         result = list(coordinates = o$summary),
         gg = function(s, yl, tt, co) gg_dimred_discrete(s, yl, tt, co, o$ax), prefix = "landscape") })
  byod_dr_cols <- apprServer("byod_dr_appr", byod_dr_levels, byod_dr_export)
  byod_dr_hm_export <- reactive({ o <- byod_dr_base()
    list(summary = o$summary, ttl = o$ttl, prefix = "landscape", diverging = FALSE,
         result = list(coordinates = o$summary),
         gg = function(mat, stops, tt, dv) gg_dimred_gradient(mat, stops, tt, dv, o$ax, o$scale_title)) })
  byod_dr_spec <- apprServerHeat("byod_dr_appr_heat", byod_dr_hm_export)
  output$byod_landscape <- renderPlotly({ d <- byod_layer_d(); req(d)
    if ((byod_dr_base()$kind %||% "discrete") == "gradient")
      byod_dimred(d, input$byod_dr_type %||% "UMAP", input$byod_dr_col %||% "Timepoint", input$byod_dr_gene,
                  scale_stops = byod_dr_spec()$stops)$plot
    else byod_dimred(d, input$byod_dr_type %||% "UMAP", input$byod_dr_col %||% "Timepoint", input$byod_dr_gene,
                     cols = byod_dr_cols())$plot })

  byod_ge_obj <- reactive({ d <- byod_layer_d(); req(d, input$byod_gene)
    byod_gene_tc(d, input$byod_gene, input$byod_split %||% "None") })
  byod_ge_levels <- reactive(byod_ge_obj()$splits %||% character(0))
  byod_ge_export <- reactive({ o <- byod_ge_obj()
    list(summary = o$summary, ylab = o$ylab, ttl = o$ttl,
         result = { d <- byod_layer_d()
                    v <- expr_vec(d, input$byod_gene)
                    if (is.null(v)) NULL
                    else list(value_per_profile = .per_sample(
                      d, stats::setNames(v, colnames(d$expr)), input$byod_gene)) },
         splitby = if (identical(input$byod_split, "None")) NULL else input$byod_split, gg = tc_gg, prefix = "trajectory") })
  byod_ge_cols <- apprServer("byod_ge_appr", byod_ge_levels, byod_ge_export)
  output$byod_gene_plot <- renderPlotly({ d <- byod_layer_d(); req(d, input$byod_gene)
    if (identical(input$byod_traj_mode, "gam"))
      byod_gene_gam(d, input$byod_gene, input$byod_split %||% "None", cols = byod_ge_cols())
    else byod_gene_tc(d, input$byod_gene, input$byod_split %||% "None", cols = byod_ge_cols())$plot })

  output$byod_de_groups <- renderUI({ d <- byod_layer_d(); req(d, input$byod_de_var)
    lv <- group_levels(d, input$byod_de_var); validate(need(length(lv) >= 2, "Need ≥2 groups."))
    layout_columns(col_widths = c(6, 6), fill = FALSE,
      selectInput("byod_de_g1", "Group 1 (reference)", lv, selected = lv[1]),
      selectInput("byod_de_g2", "Group 2", lv, selected = lv[length(lv)])) })
  byod_gen_layer <- reactiveVal(0L)
  observeEvent(byod_layer_d(), byod_gen_layer(isolate(byod_gen_layer()) + 1L), ignoreInit = TRUE)
  byod_gen_study <- reactiveVal(0L)
  observeEvent(byod_study(),   byod_gen_study(isolate(byod_gen_study()) + 1L), ignoreInit = TRUE)

  .outcome_of <- byod_outcome_of

  .log_module_outcome <- function(module, o, params = function(o) list()) {
    pr <- tryCatch(params(o), error = function(e) list())
    plan <- byod_outcome_plan(module, o, pr)
    for (rec in plan)
      do.call(log_event, c(list(rec$module), rec$params,
                           list(outcome = rec$outcome, engine = rec$engine)))
    invisible(!identical(plan[[1]]$outcome, "declined"))
  }

  resettable <- function(run_reactive, gen_reactive, valueFn, ignoreNULL = TRUE,
                         module = NULL, params = function(o) list()) {
    raw <- eventReactive(run_reactive(), {
      restore_rng <- trace_rng_scope(); on.exit(restore_rng(), add = TRUE)
      if (is.null(module)) return(valueFn())

      o <- tryCatch(valueFn(),
                    validation = function(cnd) {
                      log_declined(module, conditionMessage(cnd)); stop(cnd) })
      .log_module_outcome(module, o, params)
      o
    }, ignoreNULL = ignoreNULL)
    gen_at_run <- reactiveVal(-1L)
    seed_at_run <- reactiveVal(NA_integer_)
    observeEvent(run_reactive(), {
      gen_at_run(isolate(gen_reactive()))
      seed_at_run(isolate(analysis_seed()))
    })
    reactive({ req(identical(gen_at_run(), gen_reactive()),
                    identical(seed_at_run(), analysis_seed())); raw() })
  }

  byod_de_res <- resettable(reactive(input$byod_de_run), byod_gen_layer, function() { d <- byod_layer_d(); req(d, input$byod_de_g1, input$byod_de_g2)
    byod_run_limma(d, input$byod_de_var, input$byod_de_g1, input$byod_de_g2,
                   subject_col = byod_subj(), paired = isTRUE(input$byod_de_paired),
                   filter_constant = isTRUE(input$byod_de_filter_constant)) },
    module = "DE (limma)",
    params = function(o) list(layer = isolate(input$byod_layer), var = isolate(input$byod_de_var),
      g1 = isolate(input$byod_de_g1), g2 = isolate(input$byod_de_g2),
      paired = isTRUE(attr(o, "paired")), unit = attr(o, "unit") %||% "-",
      n_units = attr(o, "n_units") %||% NA, n_profiles = attr(o, "n") %||% NA,
      excluded_profiles = attr(o, "excluded_profiles") %||% 0L,
      excluded_subjects = paste(attr(o, "excluded_subjects"), collapse = ";"),
      filter_constant = isTRUE(attr(o, "filter_constant")),
      removed_features = length(attr(o, "removed_features")),
      block = attr(o, "block") %||% "-", fdr = isolate(input$byod_de_fdr),
      lfc = isolate(input$byod_de_lfc)))

  byod_de_lbl <- reactive({ o <- byod_de_res()
    nu <- attr(o, "n_units"); np <- attr(o, "n"); ul <- attr(o, "unit") %||% "unit"
    ns <- attr(o, "n_subjects"); sl <- attr(o, "subject_level")
    suffix <- if (is.null(nu) || is.na(nu)) sprintf(" (n=%d profiles)", np)
      else if (identical(ul, "profile")) sprintf(" (n=%d profiles, unpaired)", nu)
      else if (isTRUE(attr(o, "nested")) && !is.null(ns))

        sprintf(" (n=%d %s%s; %d paired %s%s nested inside them; %d profiles)",
                nu, ul, if (nu == 1L) "" else "s", ns, sl, if (ns == 1L) "" else "s", np)
      else sprintf(" (n=%d %s%s, %s; %d profiles)", nu, ul, if (nu == 1L) "" else "s",
                   if (isTRUE(attr(o, "paired"))) "paired" else "unpaired", np)
    ne <- attr(o, "excluded_profiles") %||% 0L
    exclusion <- if (ne > 0L) sprintf("; excluded %d profiles from %d incomplete subjects",
      ne, length(attr(o, "excluded_subjects"))) else ""
    paste0(input$byod_de_g2, " vs ", input$byod_de_g1, suffix, exclusion) })
  byod_de_levels <- reactive({ byod_de_res(); VOLCANO_LEVELS })
  output$byod_de_status <- renderUI({
    o <- byod_de_res(); req(o)
    if (!nrow(o)) return(div(class = "alert alert-warning py-2", attr(o, "declined") %||% "No result."))
    div(class = "small text-muted", byod_de_lbl(),
      if (length(attr(o, "excluded_subjects"))) tags$details(tags$summary("Excluded subject IDs"),
        paste(attr(o, "excluded_subjects"), collapse = ", ")),
      p(sprintf("Constant-feature filter: %s; %d features excluded.",
        if (isTRUE(attr(o, "filter_constant"))) "on" else "off", length(attr(o, "removed_features")))),
      if (length(attr(o, "constant_features"))) p(sprintf(
        "%d features are constant across the selected profiles; their IDs and removal status are included in the result archive.",
        length(attr(o, "constant_features")))))
  })
  byod_de_export <- reactive({ req(byod_de_res()); validate(need(nrow(byod_de_res()) > 0, "No results."))
    s <- plot_volcano(byod_de_res(), byod_de_lbl(), input$byod_de_fdr, input$byod_de_lfc, n_label = input$byod_de_nlab %||% 0)$summary
    list(summary = s, ylab = "-log10 P", ttl = byod_de_lbl(), splitby = NULL,
         csv = as.data.frame(s)[, c("Feature","logFC","AveExpr","t","P.Value","adj.P.Val","sig")],

         result = { o <- byod_de_res()
                    list(differential = o,
                         excluded_subjects = data.frame(subject = attr(o, "excluded_subjects") %||% character()),
                         constant_features = data.frame(
                           Feature = attr(o, "constant_features") %||% character(),
                           removed = (attr(o, "constant_features") %||% character()) %in% attr(o, "removed_features")),
                         contrast = data.frame(
                           variable = isolate(input$byod_de_var),
                           tested = isolate(input$byod_de_g2),
                           reference = isolate(input$byod_de_g1),
                           paired = isTRUE(attr(o, "paired")),
                           unit = attr(o, "unit") %||% NA_character_,
                           n_units = attr(o, "n_units") %||% NA_integer_,
                           n_profiles = attr(o, "n") %||% NA_integer_,
                           block = attr(o, "block") %||% NA_character_,
                           block_source = if (isTRUE(attr(o, "paired"))) "declared" else "none",
                           excluded_profiles = attr(o, "excluded_profiles") %||% 0L,
                           filter_constant = isTRUE(attr(o, "filter_constant")),
                           n_features_removed = length(attr(o, "removed_features")),
                           n_features_tested = nrow(o),
                           layer = isolate(input$byod_layer),
                           study = isolate(byod_study()$name),
                           fdr = isolate(input$byod_de_fdr),
                           lfc = isolate(input$byod_de_lfc),
                           stringsAsFactors = FALSE)) },
         gg = gg_volcano, prefix = paste0("volcano_", gsub("[^A-Za-z0-9]+", "_", byod_de_lbl()))) })
  byod_de_cols <- apprServer("byod_de_appr", byod_de_levels, byod_de_export)
  output$byod_volcano <- renderPlotly({ req(byod_de_res())
    validate(need(nrow(byod_de_res()) > 0, attr(byod_de_res(), "declined") %||% "No results."))
    plot_volcano(byod_de_res(), byod_de_lbl(), input$byod_de_fdr, input$byod_de_lfc,
                 cols = byod_de_cols(), n_label = input$byod_de_nlab %||% 0)$plot })
  output$byod_de_table <- renderDT({ req(byod_de_res())
    dt <- byod_de_res(); dt$logFC <- round(dt$logFC, 3); dt$AveExpr <- round(dt$AveExpr, 3)
    dt$t <- round(dt$t, 2); dt$P.Value <- signif(dt$P.Value, 3); dt$adj.P.Val <- signif(dt$adj.P.Val, 3)
    datatable(dt, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })
  output$byod_de_dl <- downloadHandler("byod_DE_results.csv", function(f) write.csv(byod_de_res(), f, row.names = FALSE))

  output$byod_en_groups <- renderUI({ d <- byod_layer_d(); req(d, input$byod_en_var)
    lv <- group_levels(d, input$byod_en_var); validate(need(length(lv) >= 2, "Need ≥2 groups."))
    layout_columns(col_widths = c(6, 6), fill = FALSE,
      selectInput("byod_en_g1", "Group 1 (reference)", lv, selected = lv[1]),
      selectInput("byod_en_g2", "Group 2", lv, selected = lv[length(lv)])) })
  byod_en_res <- resettable(reactive(input$byod_en_run), byod_gen_layer, function() { d <- byod_layer_d(); req(d, input$byod_en_coll)
    if (input$byod_en_mode == "gsea") { req(input$byod_en_g1, input$byod_en_g2)
      tt <- withProgress(message = "DE + GSEA…", byod_run_limma(d, input$byod_en_var, input$byod_en_g1,
              input$byod_en_g2, byod_subj(), isTRUE(input$byod_en_paired)))
      validate(need(nrow(tt) > 0, "Too few samples for this contrast."))
      list(kind = "gsea", res = run_gsea(setNames(tt$t, tt$Feature), input$byod_en_coll, d = d),
           lab = paste0(input$byod_en_coll, ": ", input$byod_en_g2, " vs ", input$byod_en_g1))
    } else { g <- match_feats(d, parse_genes(input$byod_en_genes))
      list(kind = "ora", res = run_ora(g, input$byod_en_coll, rownames(d$expr), d = d),
           lab = paste0(input$byod_en_coll, " · ORA (", length(g), " genes)")) } },
    module = "Enrichment",
    params = function(o) list(.engine = .mode_engine(isolate(input$byod_en_mode)),
      layer = isolate(input$byod_layer), mode = isolate(input$byod_en_mode),
      collection = isolate(input$byod_en_coll)))
  byod_en_obj <- reactive({ r <- byod_en_res(); if (r$kind == "gsea") plot_gsea(r$res, 20, r$lab) else plot_ora(r$res, 20, r$lab) })
  byod_en_levels <- reactive({ o <- byod_en_obj(); if (length(o$splits)) o$splits else character(0) })
  byod_en_export <- reactive({ r <- byod_en_res(); o <- byod_en_obj()
    list(summary = o$summary, ylab = o$ylab, ttl = o$ttl, splitby = NULL,
         gg = if (r$kind == "gsea") gg_gsea else gg_ora, csv = o$summary, prefix = paste0("enrichment_", r$kind)) })
  byod_en_cols <- apprServer("byod_en_appr", byod_en_levels, byod_en_export)
  output$byod_en_plot  <- renderPlotly({ req(byod_en_res()); r <- byod_en_res()
    if (r$kind == "gsea") plot_gsea(r$res, 20, r$lab, cols = byod_en_cols())$plot
    else plot_ora(r$res, 20, r$lab, cols = byod_en_cols())$plot })
  output$byod_en_table <- renderDT({ req(byod_en_res())
    datatable(byod_en_res()$res, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })

  output$byod_int_head <- renderUI({ s <- byod_study()
    if (is.null(s) || length(s$layer_names) < 2)
      div(class = "alert alert-info py-2", style = "font-size:.85rem;",
          "Load or upload a study with TWO layers (e.g. mRNA + Protein) to correlate their fold-changes.")
    else NULL })
  output$byod_int_groups <- renderUI({ s <- byod_study(); req(s)
    d <- s$layers[[input$byod_int_l1 %||% s$layer_names[1]]]; req(d, input$byod_int_var)
    lv <- group_levels(d, input$byod_int_var); validate(need(length(lv) >= 2, "Need ≥2 groups."))
    layout_columns(col_widths = c(3, 3, 6), fill = FALSE,
      selectInput("byod_int_g1", "Group 1 (reference)", lv, selected = lv[1]),
      selectInput("byod_int_g2", "Group 2", lv, selected = lv[length(lv)])) })
  byod_int_res <- resettable(reactive(input$byod_int_run), byod_gen_study, function() { s <- byod_study(); req(s, length(s$layer_names) >= 2)
    d1 <- s$layers[[input$byod_int_l1]]; d2 <- s$layers[[input$byod_int_l2]]
    req(d1, d2, input$byod_int_g1, input$byod_int_g2)
    withProgress(message = "Integrating layers", value = .5,
      byod_integrate(d1, d2, input$byod_int_var, input$byod_int_g1, input$byod_int_g2,
                     s$subject, isTRUE(input$byod_int_paired), input$byod_int_l1, input$byod_int_l2,
                     topn = as.integer(input$byod_int_topn %||% 10000))) },
    module = "Integration",
    params = function(o) list(x = isolate(input$byod_int_l1), y = isolate(input$byod_int_l2),
      var = isolate(input$byod_int_var), contrast = paste(isolate(input$byod_int_g2), "vs", isolate(input$byod_int_g1)),
      requested_pairing = isolate(isTRUE(input$byod_int_paired)), display_topn = isolate(input$byod_int_topn %||% 10000), .engine = c("limma", "stats")))
  byod_int_levels <- reactive(byod_int_res()$splits %||% character(0))
  byod_int_export <- reactive({ o <- byod_int_res(); s <- o$summary
    if (nrow(s)) { s$logFC_mRNA <- s$logFC_1; s$logFC_protein <- s$logFC_2 }
    list(summary = s, ylab = o$ylab, ttl = o$ttl, splitby = NULL, gg = gg_mp, csv = o$tab, prefix = "integration") })
  byod_int_cols <- apprServer("byod_int_appr", byod_int_levels, byod_int_export)
  output$byod_int_plot  <- renderPlotly({ req(byod_int_res()); byod_mp_plot(byod_int_res(), byod_int_cols()) })
  output$byod_int_table <- renderDT({ req(byod_int_res())
    datatable(byod_int_res()$tab, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })

  output$byod_conc_head <- renderUI({ s <- byod_study()
    if (is.null(s) || length(s$layer_names) < 2 || length(s$matched_samples) < 4)
      div(class = "alert alert-info py-2", style = "font-size:.85rem;",
          "Concordance needs TWO sample-matched layers (identical sample IDs in both). Try the ",
          tags$b("Matched kidney (mRNA+protein)"), " demo, or upload matched layers.")
    else div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
          sprintf("%d samples matched across %s; %d shared features — each feature in one layer is correlated with the same feature in the other, across samples.",
                  length(s$matched_samples), paste(s$layer_names, collapse = " + "), length(s$shared_genes))) })
  byod_conc_res <- resettable(reactive(input$byod_conc_run), byod_gen_study, function() {
    s <- byod_study(); req(s, length(s$layer_names) >= 2)
    l1 <- input$byod_conc_l1 %||% s$layer_names[1]
    l2 <- input$byod_conc_l2 %||% s$layer_names[2]
    req(l1 %in% s$layer_names, l2 %in% s$layer_names)
    validate(need(!identical(l1, l2), "Choose two different layers."))
    d1 <- s$layers[[l1]]; d2 <- s$layers[[l2]]

    withProgress(message = "Per-feature concordance", value = .5,
      byod_concordance(d1, d2,
                       shared_samples = intersect(colnames(d1$expr), colnames(d2$expr)),
                       shared_genes  = intersect(rownames(d1$expr), rownames(d2$expr)))) },
    module = "Concordance",
    params = function(o) list(study = isolate(byod_study()$name)))
  observeEvent(byod_conc_res(), { r <- byod_conc_res(); req(isTRUE(r$ok))
    updateSelectizeInput(session, "byod_conc_gene", choices = r$tab$gene, selected = r$tab$gene[1], server = TRUE)
  })
  output$byod_conc_note <- renderUI({ r <- byod_conc_res(); if (is.null(r)) return(NULL)
    if (!isTRUE(r$ok)) return(div(class = "alert alert-warning py-2 mt-2", style = "font-size:.85rem;", r$msg))
    div(class = "alert alert-success py-1 mt-2", style = "font-size:.8rem;",
        sprintf("median r = %.2f · %d concordant · %d discordant · %d weak positive", r$med_rho,
                sum(r$tab$class == "concordant"), sum(r$tab$class == "discordant"), sum(r$tab$class == "weak positive"))) })
  output$byod_conc_hist <- renderPlotly({ req(byod_conc_res()); byod_conc_hist(byod_conc_res()) })
  output$byod_conc_scatter <- renderPlotly({ r <- byod_conc_res(); s <- byod_study(); req(isTRUE(r$ok), input$byod_conc_gene)
    byod_conc_scatter(s$layers[[s$layer_names[1]]], s$layers[[s$layer_names[2]]], input$byod_conc_gene,
                      s$matched_samples, s$layer_names[1], s$layer_names[2]) })
  output$byod_conc_table <- renderDT({ r <- byod_conc_res(); req(isTRUE(r$ok))
    datatable(r$tab, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) })
  output$byod_conc_ora <- renderPlotly({ r <- byod_conc_res(); req(isTRUE(r$ok), input$byod_conc_coll)
    validate(need(length(r$discordant_genes) >= 3, "Too few discordant genes for enrichment"))
    d <- byod_study()$layers[[input$byod_conc_l1 %||% byod_study()$layer_names[1]]]
    plot_ora(run_ora(r$discordant_genes, input$byod_conc_coll, r$universe, d = d), 20, "Discordant genes")$plot })

  output$byod_diablo_head <- renderUI({ s <- byod_study()
    if (is.null(s) || length(s$layer_names) < 2 || length(s$matched_samples) < 8)
      div(class = "alert alert-info py-2", style = "font-size:.85rem;",
          "DIABLO needs TWO sample-matched layers with >= 8 matched samples. Try the ",
          tags$b("Matched kidney (mRNA+protein)"), " demo.")
    else div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
          sprintf("Multiblock sPLS-DA integrating %s across %d matched samples — an EXPLORATORY supervised projection of correlated cross-omics structure. Fitted on all samples without cross-validation or permutation, so read it as descriptive, not a validated predictive/discriminative signature.",
                  paste(s$layer_names[1:2], collapse = " + "), length(s$matched_samples))) })
  byod_diablo_res <- resettable(reactive(input$byod_diablo_run), byod_gen_study, function() { s <- byod_study(); req(s, length(s$layer_names) >= 2)
    withProgress(message = "Fitting DIABLO", detail = "block sPLS-DA…", value = .5,
      tryCatch(byod_diablo(s, s$layer_names[1], s$layer_names[2], input$byod_diablo_out %||% "Timepoint",
                           nfeat = input$byod_diablo_nfeat %||% 400, keepx = input$byod_diablo_keepx %||% 20),
               error = function(e) { byod_msg(list(type = "danger", txt = paste("DIABLO:", conditionMessage(e)))); NULL })) },
    module = "DIABLO",
    params = function(o) list(out = isolate(input$byod_diablo_out), nfeat = isolate(input$byod_diablo_nfeat),
      cross_validated = FALSE))
  output$byod_diablo_note <- renderUI({ r <- byod_diablo_res(); if (is.null(r)) return(NULL)
    div(class = "alert alert-success py-1", style = "font-size:.8rem;",
        sprintf("Fitted (exploratory, no CV): %d component(s) · %d samples · outcome '%s' (%d groups).", r$ncomp, r$n, r$outcome, nlevels(r$Y))) })
  output$byod_diablo_plot <- renderPlot({ r <- byod_diablo_res(); req(r); byod_diablo_plot(r, input$byod_diablo_which %||% "circos") })

  byod_diablo_cv_res <- resettable(reactive(input$byod_diablo_cv_run), byod_gen_study, function() {
    s <- byod_study(); req(s, length(s$layer_names) >= 2)
    withProgress(message = "Cross-validating DIABLO", value = 0.4,
      tryCatch(byod_diablo_cv(s, s$layer_names[1], s$layer_names[2],
                              input$byod_diablo_out %||% "Timepoint",
                              nfeat = input$byod_diablo_nfeat %||% 400,
                              keepx = input$byod_diablo_keepx %||% 20,
                              folds = input$byod_diablo_folds %||% 5,
                              repeats = input$byod_diablo_reps %||% 5),
               error = function(e) list(ok = FALSE, msg = conditionMessage(e))))
  },
    module = "DIABLO cross-validation",
    params = function(o) list(folds = o$folds %||% NA, repeats = o$repeats %||% NA, unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  output$byod_diablo_cv_note <- renderUI({ r <- byod_diablo_cv_res(); if (is.null(r)) return(NULL)
    if (!isTRUE(r$ok))
      return(div(class = "alert alert-warning py-1", style = "font-size:.8rem;",
                 r$msg %||% "Cross-validation unavailable."))
    rng <- if (length(r$per_repeat_balanced_accuracy) > 1)
      sprintf(" Balanced accuracy across the %d repeats ranged %.1f%%-%.1f%%.",
              length(r$per_repeat_balanced_accuracy),
              100 * min(r$per_repeat_balanced_accuracy),
              100 * max(r$per_repeat_balanced_accuracy)) else ""
    div(class = "alert alert-info py-1", style = "font-size:.8rem;",
        paste0(r$note, rng)) })
  output$byod_diablo_cv_tab <- renderDT({ r <- byod_diablo_cv_res(); req(isTRUE(r$ok))
    cm <- as.data.frame.matrix(r$confusion)
    cm <- cbind(observed = rownames(cm), cm, recall = sprintf("%.0f%%", 100 * r$per_class))
    datatable(cm, rownames = FALSE, options = list(dom = "t", paging = FALSE, ordering = FALSE),
              caption = sprintf("Confusion matrix pooled over %d folds x %d repeats; rows are the observed %s.",
                                r$folds, r$repeats, input$byod_diablo_out %||% "Timepoint")) })

  output$byod_hm_complex <- renderPlot({ trace_seed(); input$byod_hm_run; d <- byod_layer_d(); req(d)
    isolate(byod_cheatmap(d, match_feats(d, parse_genes(input$byod_hm_genes)),
                          input$byod_hm_by %||% "Timepoint", input$byod_hm_k %||% 3)) })

  byod_wg_res <- resettable(reactive(input$byod_wg_run), byod_gen_layer, function() { d <- byod_layer_d(); req(d)
    pw <- input$byod_wg_power; pw <- if (is.null(pw) || pw < 1) NULL else pw
    withProgress(message = "Building co-expression network", detail = "adjacency + modules…", value = .5,
      byod_coexpr(d, nGenes = input$byod_wg_ngenes %||% 2000, minModuleSize = input$byod_wg_minsize %||% 20, power = pw)) },
    module = "WGCNA",
    params = function(o) list(layer = isolate(input$byod_layer), ngenes = isolate(input$byod_wg_ngenes),
      power = isolate(input$byod_wg_power)))
  byod_wg_ok <- reactive({ r <- byod_wg_res(); if (is.null(r) || isFALSE(r$ok)) NULL else r })
  output$byod_wg_note <- renderUI({ r <- byod_wg_res(); if (is.null(r)) return(NULL)
    if (isFALSE(r$ok))
      return(div(class = "alert alert-warning py-2 px-3", style = "font-size:.85rem;",
                 tags$b("No network built. "), r$msg))
    div(class = "alert alert-secondary py-1", style = "font-size:.8rem;",
        sprintf("Engine: %s · soft power %s · %d modules from %d genes.", r$engine, r$power, ncol(r$ME), r$nGenes)) })
  byod_wg_hm_export <- reactive({ r <- byod_wg_ok(); req(r)
    list(summary = r$MT, gg = gg_heatmap, ttl = "Module-trait correlation", prefix = "module_trait", diverging = TRUE,

         result = list(module_trait_r = r$MT, module_trait_p = r$MTp,
                       modules = byod_module_summary(r),
                       feature_module = data.frame(feature = names(r$colors),
                                                   module = unname(r$colors),
                                                   stringsAsFactors = FALSE),
                       eigengenes = r$ME)) })
  byod_wg_spec <- apprServerHeat("byod_wg_appr", byod_wg_hm_export)
  output$byod_wg_plot  <- renderPlotly({ r <- byod_wg_ok(); req(r); byod_mt_heatmap(r, byod_wg_spec()$colorscale) })

  output$byod_wg_plot_ui <- renderUI({
    r <- byod_wg_ok(); req(r)
    plotlyOutput("byod_wg_plot", height = paste0(byod_mt_height(r), "px"))
  })
  output$byod_wg_table <- renderDT({ r <- byod_wg_ok(); req(r)
    datatable(byod_module_summary(r), rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) })

  observeEvent(byod_layer_d(), { d <- byod_layer_d(); req(d)
    lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
    if (length(lv) < 2) return()
    updateSelectInput(session, "byod_mc_ref", choices = lv, selected = lv[1])
    updateSelectizeInput(session, "byod_mc_tps", choices = lv, selected = lv[-1], server = TRUE) })

  observeEvent(input$byod_mc_ref, { d <- byod_layer_d(); req(d, input$byod_mc_ref)
    lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
    sel <- setdiff(input$byod_mc_tps %||% lv, input$byod_mc_ref)
    updateSelectizeInput(session, "byod_mc_tps", choices = lv,
                         selected = if (length(sel)) sel else setdiff(lv, input$byod_mc_ref),
                         server = TRUE) }, ignoreInit = TRUE)

  byod_mc_res <- resettable(reactive(input$byod_mc_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d, input$byod_mc_ref, input$byod_mc_tps)
    withProgress(message = "Multi-timepoint contrasts",
                 detail = sprintf("%d contrasts against %s...",
                                  length(setdiff(input$byod_mc_tps, input$byod_mc_ref)),
                                  input$byod_mc_ref),
                 value = .4,
      byod_multi_contrast(d, input$byod_mc_ref, input$byod_mc_tps,
                          fdr = input$byod_mc_fdr %||% .05,
                          lfc = input$byod_mc_lfc %||% 0)) },
    module = "Multi-timepoint",
    params = function(o) list(control = isolate(input$byod_mc_ctrl), unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  byod_mc_ok <- reactive({ r <- byod_mc_res(); if (is.null(r) || !isTRUE(r$ok)) NULL else r })
  output$byod_mc_note <- renderUI({ r <- byod_mc_res()
    if (is.null(r)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Declare a control timepoint and run every other timepoint against it in one pass. Each contrast is the same one the Differential tab computes, at the same resolved unit, so the numbers match; what this adds is one volcano per timepoint and the intersection structure across them."))
    if (!isTRUE(r$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", r$msg))
    ui <- list(tags$b(r$ttl), " ", r$note)
    if (length(r$declined))
      ui <- c(ui, list(tags$div(class = "mt-1", style = "font-size:.8rem;", tags$b("Not shown: "),
        paste(sprintf("%s (%s)", names(r$declined), unlist(r$declined)), collapse = "; "))))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;", ui) })
  output$byod_mc_volcano_ui <- renderUI({ r <- byod_mc_ok(); req(r)
    n <- length(r$test_tps); nr <- if (n <= 6) 1 else ceiling(n / ceiling(n / 2))
    plotlyOutput("byod_mc_volcano", height = paste0(150 + 260 * nr, "px")) })
  output$byod_mc_volcano <- renderPlotly({ r <- byod_mc_ok(); req(r); byod_multi_volcano_plot(r) })
  output$byod_mc_upset <- renderPlotly({ r <- byod_mc_ok(); req(r); byod_multi_upset_plot(r) })
  output$byod_mc_stats <- renderDT({ r <- byod_mc_ok(); req(r)
    st <- r$stats[, c("Timepoint", "n_sig", "n_up", "n_down", "n_units")]
    st$`only here` <- as.integer(r$n_only[as.character(st$Timepoint)])
    names(st) <- c("Timepoint", "significant", "up", "down", "n units", "only here")
    datatable(st, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE, dom = "t")) })
  output$byod_mc_core <- renderDT({ r <- byod_mc_ok(); req(r)
    if (!length(r$core)) return(datatable(data.frame(Feature = character(0)), rownames = FALSE))
    first <- r$results[[1]]
    tb <- first[match(r$core, first$Feature), c("Feature", "logFC", "adj.P.Val"), drop = FALSE]
    names(tb) <- c("Feature", sprintf("logFC at %s", r$test_tps[1]), "FDR there")
    tb[[2]] <- round(tb[[2]], 3); tb[[3]] <- signif(tb[[3]], 3)
    datatable(tb[order(tb[[3]]), ], rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE)) })

  output$byod_mc_vdl <- downloadHandler(
    filename = function() sprintf("volcanoes_vs_%s.zip",
                                  gsub("[^A-Za-z0-9._-]+", "_", byod_mc_ok()$reference)),
    content = function(f) {
      r <- byod_mc_ok(); req(r)
      td <- file.path(tempdir(), paste0("volc_", as.integer(Sys.time())))
      on.exit(unlink(td, recursive = TRUE), add = TRUE)
      withProgress(message = "Rendering volcanoes", detail = "each panel, then the composite...",
                   value = .4, {
        files <- byod_multi_volcano_export(r, td)
        writeLines(c(
          sprintf("Volcano plots, every timepoint against %s", r$reference),
          sprintf("FDR < %.3g%s", r$fdr,
                  if (r$lfc > 0) sprintf(", |log2FC| >= %g", r$lfc) else ""),
          sprintf("Unit: %s %s%s", r$n_units, r$unit, if (isTRUE(r$nested)) ", nested" else ""),
          "",
          "volcano_<timepoint>_vs_<control>  one panel on its own",
          "volcano_ALL_vs_<control>          the composite",
          "",
          "Every panel uses the SAME fold-change limits as the composite, so a",
          "single panel remains comparable with the others.", "",
          sprintf("Generated by TRACE-Omics on %s", Sys.Date())),
          file.path(td, "README.txt"))
        files <- c(files, file.path(td, "README.txt"))
      })
      owd <- setwd(td); on.exit(setwd(owd), add = TRUE)
      utils::zip(f, basename(files), flags = "-q9X")
    })

  output$byod_mc_dl <- downloadHandler("multi_timepoint_contrasts.csv", function(f) {
    r <- byod_mc_ok(); req(r)
    out <- do.call(rbind, lapply(names(r$results), function(tp) {
      x <- r$results[[tp]]
      data.frame(Timepoint = tp, Feature = x$Feature, logFC = x$logFC,
                 P.Value = x$P.Value, adj.P.Val = x$adj.P.Val, stringsAsFactors = FALSE) }))
    write.csv(out, f, row.names = FALSE) })

  observeEvent(byod_layer_d(), { d <- byod_layer_d(); req(d)
    lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
    if (length(lv) < 2) return()
    updateSelectInput(session, "byod_rc_ref", choices = lv, selected = lv[1])
    updateSelectizeInput(session, "byod_rc_tps", choices = lv, selected = lv[-1], server = TRUE)

    two <- names(d$meta)[vapply(names(d$meta), function(v) {
      z <- unique(as.character(d$meta[[v]])); z <- z[!is.na(z)]
      length(z) == 2L && !v %in% c("Timepoint", d$subject %||% "") }, logical(1))]
    updateSelectInput(session, "byod_rc_group", choices = c("none" = "", two)) })
  byod_rc_res <- resettable(reactive(input$byod_rc_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d, input$byod_rc_ref, input$byod_rc_tps)
    withProgress(message = "Response composition", detail = "contrasts per arm per timepoint...",
                 value = .4,
      byod_response_composition(d, input$byod_rc_ref, input$byod_rc_tps,
                                group_var = input$byod_rc_group)) },
    module = "Response composition",
    params = function(o) list(split = isolate(input$byod_rc_split), unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  byod_rc_ok <- reactive({ r <- byod_rc_res(); if (is.null(r) || !isTRUE(r$ok)) NULL else r })
  output$byod_rc_note <- renderUI({ r <- byod_rc_res()
    if (is.null(r)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Counting significant features answers 'how many'. This answers 'how much of the layer, and is it the same part': at every timepoint each feature is either changed from baseline or not, and with a two-level split, changed in both arms, in one, or in neither."))
    if (!isTRUE(r$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", r$msg))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;", tags$b(r$ttl), " ", r$note) })
  output$byod_rc_plot <- renderPlotly({ r <- byod_rc_ok(); req(r)
    byod_response_composition_plot(r, as_percent = identical(input$byod_rc_unit %||% "pct", "pct")) })
  output$byod_rc_table <- renderDT({ r <- byod_rc_ok(); req(r)
    tb <- r$tab; tb$pct <- round(tb$pct, 2)
    datatable(tb[, c("Timepoint", "category", "n", "pct")], rownames = FALSE,
              options = list(pageLength = 12, scrollX = TRUE, dom = "t")) })

  observeEvent(list(byod_layer_d(), byod_study()), {
    d <- byod_layer_d(); s <- byod_study(); req(d, s)
    lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
    v <- c(s$clinical_num, s$clinical_cat)
    if (!length(v)) v <- setdiff(names(d$meta), c("Timepoint", d$subject %||% "", "sample_id"))
    updateSelectizeInput(session, "byod_vs_vars", choices = v, selected = head(v, 8), server = TRUE)
    updateSelectizeInput(session, "byod_vs_tps", choices = lv, selected = lv, server = TRUE) })
  byod_vs_res <- resettable(reactive(input$byod_vs_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d, input$byod_vs_vars)
    how <- input$byod_vs_split %||% "route"
    sp <- if (identical(how, "route")) NULL
          else setNames(rep(list(how), length(input$byod_vs_vars)), input$byod_vs_vars)

    withProgress(message = "Variable screen",
                 detail = sprintf("%d variable%s x %d timepoint%s = %d cells",
                                  length(input$byod_vs_vars),
                                  if (length(input$byod_vs_vars) == 1) "" else "s",
                                  length(input$byod_vs_tps %||% 1),
                                  if (length(input$byod_vs_tps %||% 1) == 1) "" else "s",
                                  length(input$byod_vs_vars) * length(input$byod_vs_tps %||% 1)),
                 value = .4,
      byod_variable_screen(d, input$byod_vs_vars, input$byod_vs_tps, splits = sp,
                           min_per_group = as.integer(input$byod_vs_min_group %||% 3L),
                           min_n_cont    = as.integer(input$byod_vs_min_cont  %||% 8L))) },
    module = "Subject variables: screen",

    params = function(o) list(.engine = c(if ((o$n_continuous %||% 0) > 0) "stats",
                                        if ((o$n_categorical %||% 0) > 0) "limma"),
                              n_comparisons = o$n_comparisons %||% NA,
                              unit = o$unit %||% "-", n_units = o$n_units %||% NA,
                              continuous_handling = isolate(input$byod_vs_split %||% "route"),
                              n_continuous = o$n_continuous %||% NA,
                              n_categorical = o$n_categorical %||% NA,
                              min_subjects_continuous = isolate(input$byod_vs_min_cont %||% 8L),
                              min_per_group = isolate(input$byod_vs_min_group %||% 3L),
                              not_estimable = sum(is.na(o$tab$n_sig))))
  byod_vs_ok <- reactive({ r <- byod_vs_res(); if (is.null(r) || !isTRUE(r$ok)) NULL else r })
  output$byod_vs_note <- renderUI({ r <- byod_vs_res()
    if (is.null(r)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Each selected variable is analysed INSIDE each timepoint separately, by the route its measurement scale allows. A continuous variable is correlated with each feature on its own scale (Spearman, one subject per value, BH across features). A categorical variable is contrasted between its levels at the resolved unit. Splitting a continuous variable into groups is available above and is recorded when used; it is not the default, because a median split and an upper-quartile split answer different questions and both discard the ordering."))
    if (!isTRUE(r$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", r$msg))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;", tags$b(r$ttl), " ", r$note) })
  output$byod_vs_plot_ui <- renderUI({ r <- byod_vs_ok(); req(r)
    plotlyOutput("byod_vs_plot",
                 height = paste0(180 + 26 * min(40L, r$n_comparisons), "px")) })
  output$byod_vs_plot <- renderPlotly({ r <- byod_vs_ok(); req(r)
    byod_variable_screen_plot(r, input$byod_vs_value %||% "count") })
  output$byod_vs_table <- renderDT({ r <- byod_vs_ok(); req(r)
    tb <- r$tab; tb$pct <- round(tb$pct, 3)
    datatable(tb, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE)) })

  byod_tcon_res <- resettable(reactive(input$byod_tcon_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d)
    withProgress(message = "Trajectory consistency", detail = "leave-one-out correlations...", value = .4,
      byod_traj_consistency(d, min_tp = as.integer(input$byod_tcon_mintp %||% 3),
                            max_feat = as.integer(input$byod_tcon_nfeat %||% 4000))) },
    module = "Trajectory consistency",
    params = function(o) list(unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  byod_tcon_ok <- reactive({ r <- byod_tcon_res(); if (is.null(r) || !isTRUE(r$ok)) NULL else r })
  output$byod_tcon_note <- renderUI({ r <- byod_tcon_res()
    if (is.null(r)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Asks how reproducible each feature's time-course is across the subjects that produced it: every subject's trajectory is correlated with the mean of the others, leaving the subject (and anything nested with it) out of that mean. Run to compute."))
    if (!isTRUE(r$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", r$msg))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;",
        paste0(r$ttl, ". ", r$note)) })
  output$byod_tcon_plot <- renderPlotly({ r <- byod_tcon_ok(); req(r); byod_traj_consistency_plot(r) })
  output$byod_tcon_table <- renderDT({ r <- byod_tcon_ok(); req(r)
    tb <- utils::tail(r$tab, 200)[order(utils::tail(r$tab, 200)$mean_r), ]
    for (cn in c("mean_r", "sd_r", "min_r")) tb[[cn]] <- round(tb[[cn]], 3)
    datatable(tb, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })
  observeEvent(byod_tcon_res(), { r <- byod_tcon_ok(); req(r)
    updateSelectizeInput(session, "byod_tcon_feat", choices = r$tab$Feature,
                         selected = r$tab$Feature[1], server = TRUE) })
  output$byod_tcon_one <- renderPlotly({ r <- byod_tcon_ok(); d <- byod_layer_d()
    req(r, d, input$byod_tcon_feat)
    byod_gene_tc(d, input$byod_tcon_feat, r$subj_col)$plot })
  output$byod_tcon_dl <- downloadHandler("trajectory_consistency.csv", function(f) {
    r <- byod_tcon_ok(); req(r); write.csv(r$tab, f, row.names = FALSE) })

  observeEvent(byod_layer_d(), { d <- byod_layer_d(); req(d)
    lv <- levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels)))
    if (length(lv) < 3) return()
    updateSelectInput(session, "byod_bs_a", choices = lv, selected = lv[1])
    updateSelectInput(session, "byod_bs_b", choices = lv, selected = lv[2])
    updateSelectizeInput(session, "byod_bs_test", choices = lv,
                         selected = setdiff(lv, lv[1:2]), server = TRUE) })
  byod_bs_res <- resettable(reactive(input$byod_bs_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d, input$byod_bs_a, input$byod_bs_b, input$byod_bs_test)
    withProgress(message = "Baseline sensitivity",
                 detail = sprintf("two contrasts x %d timepoints...", length(input$byod_bs_test)),
                 value = .4,
      byod_baseline_sensitivity(d, input$byod_bs_a, input$byod_bs_b, input$byod_bs_test)) },
    module = "Baseline sensitivity",
    params = function(o) list(base1 = isolate(input$byod_bs_a), base2 = isolate(input$byod_bs_b), unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  byod_bs_ok <- reactive({ r <- byod_bs_res(); if (is.null(r) || !isTRUE(r$ok)) NULL else r })
  output$byod_bs_note <- renderUI({ r <- byod_bs_res()
    if (is.null(r)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Runs the SAME contrast against two different baselines, at every timepoint you select, and reports how much of the answer survives the choice. Both contrasts go through the ordinary differential-expression path, so each inherits the resolved unit and the same refusal rules."))
    if (!isTRUE(r$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", r$msg))
    ui <- list(tags$b(r$ttl), " ", r$note)
    if (length(r$declined))
      ui <- c(ui, list(tags$div(class = "mt-1", style = "font-size:.8rem;",
        tags$b("Not shown: "),
        paste(sprintf("%s (%s)", names(r$declined), unlist(r$declined)), collapse = "; "))))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;", ui) })

  output$byod_bs_plot_ui <- renderUI({ r <- byod_bs_ok(); req(r)
    nrow_panels <- ceiling(length(r$test_tps) / 4)
    plotlyOutput("byod_bs_plot", height = paste0(240 + 210 * nrow_panels, "px")) })
  output$byod_bs_plot <- renderPlotly({ r <- byod_bs_ok(); req(r); byod_baseline_plot(r) })
  output$byod_bs_overlap <- renderPlotly({ r <- byod_bs_ok(); req(r); byod_baseline_overlap_plot(r) })
  output$byod_bs_stats <- renderDT({ r <- byod_bs_ok(); req(r)
    st <- r$stats[, c("Timepoint", "r", "r2", "n_a", "n_b", "n_both", "pct_overlap", "n_units")]
    names(st) <- c("Timepoint", "r", "R2", paste0("sig ", r$baselines[1]),
                   paste0("sig ", r$baselines[2]), "shared", "% overlap", "n units")
    for (cn in c("r", "R2")) st[[cn]] <- round(st[[cn]], 4)
    st[["% overlap"]] <- round(st[["% overlap"]], 1)
    datatable(st, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE, dom = "t")) })
  output$byod_bs_table <- renderDT({ r <- byod_bs_ok(); req(r)
    tb <- r$tab[!r$tab$agreement %in% c("both", "neither"), , drop = FALSE]
    tb <- head(tb[order(pmin(tb$q_a, tb$q_b)), ], 500)
    for (cn in c("logFC_a", "logFC_b")) tb[[cn]] <- round(tb[[cn]], 3)
    for (cn in c("q_a", "q_b")) tb[[cn]] <- signif(tb[[cn]], 3)
    datatable(tb, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) })

  observeEvent(byod_study(), { s <- byod_study(); req(s)
    v <- c(s$clinical_num, s$clinical_cat)
    if (length(v)) updateSelectizeInput(session, "byod_pxc_vars", choices = v,
                                        selected = head(v, 6), server = TRUE) })
  byod_pxc_res <- resettable(reactive(input$byod_pxc_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d, input$byod_pxc_vars)
    g <- byod_gsva_rv()
    if (is.null(g) || is.null(g$mat))
      return(list(ok = FALSE, msg = "Compute pathway activity on the Pathways tab first; this grid tests those scores against subject characteristics."))
    withProgress(message = "Pathway x clinical", detail = "collapsing to subjects, then testing...", value = .4,
      byod_pathway_clinical(d, g$mat, input$byod_pxc_vars,
                            method = input$byod_pxc_method %||% "spearman")) },
    module = "Subject variables: pathways",
    params = function(o) list(unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  byod_pxc_ok <- reactive({ r <- byod_pxc_res(); if (is.null(r) || !isTRUE(r$ok)) NULL else r })
  output$byod_pxc_note <- renderUI({ r <- byod_pxc_res()
    if (is.null(r)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Tests every pathway score against every selected subject characteristic. Scores are collapsed to one value per subject first, so a donor variable is counted once and not once per profile, and the FDR is taken across the whole grid rather than per column."))
    if (!isTRUE(r$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", r$msg))
    directions <- unique(r$tab[r$tab$type == "two-group", c("variable", "tested", "reference"), drop = FALSE])
    div(class = "alert alert-success py-1", style = "font-size:.8rem;", paste0(r$ttl, ". ", r$note),
      if (nrow(directions)) tags$details(tags$summary("Binary effect directions"),
        tags$ul(lapply(seq_len(nrow(directions)), function(i) tags$li(paste0(directions$variable[i], ": ",
          directions$tested[i], " minus ", directions$reference[i]))))),
      if (length(r$skipped_inconsistent_variables)) p("Skipped variables that change within a subject: ",
        paste(r$skipped_inconsistent_variables, collapse = ", "))) })
  output$byod_pxc_plot <- renderPlotly({ r <- byod_pxc_ok(); req(r); byod_pathway_clinical_plot(r) })
  output$byod_pxc_table <- renderDT({ r <- byod_pxc_ok(); req(r)
    tb <- r$tab; tb$estimate <- round(tb$estimate, 3)
    tb$p <- signif(tb$p, 3); tb$fdr <- signif(tb$fdr, 3)
    datatable(tb, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) })

  byod_gsva_rv <- reactiveVal(NULL)

  observeEvent(byod_layer_d(), byod_gsva_rv(NULL), ignoreInit = TRUE)
  observeEvent(input$byod_pw_run, {
    d <- byod_layer_d(); req(d, input$byod_pw_coll)
    err <- NULL
    r <- withProgress(message = "Computing pathway activity", detail = paste0(input$byod_pw_method, "…"), value = .5,
      tryCatch(byod_pathway_activity(d, input$byod_pw_coll, input$byod_pw_method, input$byod_pw_by %||% "Timepoint"),
               error = function(e) { err <<- conditionMessage(e)
                 byod_msg(list(type = "danger", txt = paste("Pathway error:", err))); NULL }))

    if (is.null(r)) log_declined("Pathways", err %||% "no pathway could be scored on this layer")
    else .log_module_outcome("Pathways", r, function(o) list(
           layer = isolate(input$byod_layer), collection = o$collection,
           requested_method = o$requested_method %||% isolate(input$byod_pw_method),
           executed_method  = o$executed_method %||% o$method,
           by = o$by, n_sets = nrow(o$mat %||% matrix(nrow = 0))))
    byod_gsva_rv(r)
  })
  output$byod_pw_note <- renderUI({ r <- byod_gsva_rv(); if (is.null(r)) return(NULL)

    if (identical(r$ok, FALSE))
      return(div(class = "alert alert-warning py-2 px-3", style = "font-size:.82rem;",
                 tags$b("Cannot score this collection. "), r$msg))
    div(class = "alert alert-secondary py-1", style = "font-size:.8rem;",
        sprintf("Method: %s · %s · %d pathways scored across %d samples.",
                r$method %||% "?", r$collection %||% "?", nrow(r$mat), ncol(r$mat))) })
  outputOptions(output, "byod_pw_note", suspendWhenHidden = FALSE)
  output$byod_pw_pathui <- renderUI({ r <- byod_gsva_rv(); req(r)
    selectInput("byod_pw_path", "Pathway", choices = rownames(r$mat), selected = rownames(r$group_mat)[1]) })

  byod_pw_hm_export <- reactive({ r <- byod_gsva_rv(); req(r)
    if (identical(input$byod_pw_view %||% "group", "sample")) {
      d <- byod_layer_d(); req(d)
      list(summary = r$mat, prefix = "pathway_activity_per_sample", diverging = TRUE,

           result = list(scores_per_profile = r$mat, group_means = r$group_mat),
           ttl = sprintf("%s activity per sample", r$collection %||% "Pathway"),
           gg = function(mat, stops, tt, dv)
             byod_sample_heatmap_gg(mat, isolate(byod_layer_d())$meta,
                                    isolate(input$byod_pw_order) %||% "cluster",
                                    title = tt, value_name = r$method %||% "score",
                                    cluster_rows = isTRUE(isolate(input$byod_pw_rowclust)),
                                    stops = stops))
    } else {
      o <- byod_pathway_heatmap(r)
      list(summary = o$mat, gg = gg_heatmap, ttl = o$ttl,
           result = list(scores_per_profile = r$mat, group_means = r$group_mat),
           prefix = "pathway_activity", diverging = TRUE)
    } })
  byod_pw_hm_spec <- apprServerHeat("byod_pw_hm_appr", byod_pw_hm_export)

  observeEvent(byod_layer_d(), { d <- byod_layer_d(); req(d)
    cand <- byod_order_candidates(d)
    updateSelectInput(session, "byod_pw_order",
                      choices = c("cluster the samples" = "cluster", setNames(cand, cand),
                                  "leave as supplied" = "none"),
                      selected = if (length(cand)) cand[1] else "cluster") })
  output$byod_pw_heat_ui <- renderUI({
    r <- byod_gsva_rv(); req(r)
    h <- if (identical(input$byod_pw_view %||% "group", "sample")) 620 else 440
    plotlyOutput("byod_pw_heat", height = paste0(h, "px")) })
  output$byod_pw_heat <- renderPlotly({ r <- byod_gsva_rv(); req(r)
    if (identical(input$byod_pw_view %||% "group", "sample")) {
      d <- byod_layer_d(); req(d)
      g <- byod_sample_heatmap_gg(r$mat, d$meta, input$byod_pw_order %||% "cluster",
                                  title = sprintf("%s activity per sample", r$collection %||% "Pathway"),
                                  value_name = r$method %||% "score",
                                  cluster_rows = isTRUE(input$byod_pw_rowclust),

                                  stops = byod_pw_hm_spec()$stops)
      ggplotly(g) %>% ply()
    } else byod_pathway_heatmap(r, byod_pw_hm_spec()$colorscale)$plot })
  byod_pw_tc_obj <- reactive({ r <- byod_gsva_rv(); req(r, input$byod_pw_path)
    validate(need(input$byod_pw_path %in% rownames(r$mat), "Select a pathway"))
    byod_pathway_tc(byod_layer_d(), r, input$byod_pw_path, input$byod_pw_split %||% "None") })
  byod_pw_tc_levels <- reactive(byod_pw_tc_obj()$splits %||% character(0))
  byod_pw_tc_export <- reactive({ o <- byod_pw_tc_obj()
    list(summary = o$summary, ylab = o$ylab, ttl = o$ttl,
         result = { r <- byod_gsva_rv(); if (is.null(r)) NULL
                    else list(scores_per_profile = r$mat, group_means = r$group_mat) },
         splitby = if (identical(input$byod_pw_split, "None")) NULL else input$byod_pw_split, gg = tc_gg, prefix = "pathway_trajectory") })
  byod_pw_tc_cols <- apprServer("byod_pw_tc_appr", byod_pw_tc_levels, byod_pw_tc_export)
  output$byod_pw_tc <- renderPlotly({ r <- byod_gsva_rv(); req(r, input$byod_pw_path)
    validate(need(input$byod_pw_path %in% rownames(r$mat), "Select a pathway"))
    byod_pathway_tc(byod_layer_d(), r, input$byod_pw_path, input$byod_pw_split %||% "None", cols = byod_pw_tc_cols())$plot })

  byod_ms_res <- resettable(reactive(input$byod_ms_run), byod_gen_layer, function() { d <- byod_layer_d(); req(d)
    byod_module(d, match_feats(d, parse_genes(input$byod_ms_genes)), input$byod_ms_split %||% "None") }, ignoreNULL = FALSE,
    module = "Module score",
    params = function(o) list(layer = isolate(input$byod_layer), n_genes = o$n %||% NA,
      split = isolate(input$byod_ms_split)))
  byod_ms_levels <- reactive({ r <- byod_ms_res(); if (is.null(r)) character(0) else r$splits %||% character(0) })
  byod_ms_export <- reactive({ r <- byod_ms_res(); req(r)
    list(summary = r$summary, ylab = r$ylab, ttl = r$ttl,
         result = { d <- byod_layer_d()
                    g <- match_feats(d, parse_genes(input$byod_ms_genes))
                    v <- tryCatch(module_score(d, g), error = function(e) NULL)
                    if (is.null(v)) NULL else list(score_per_profile = .per_sample(d, v, "module_score"),
                                                   features_used = data.frame(feature = g, stringsAsFactors = FALSE)) },
         splitby = if (identical(input$byod_ms_split, "None")) NULL else input$byod_ms_split, gg = tc_gg, prefix = "module_score") })
  byod_ms_cols <- apprServer("byod_ms_appr", byod_ms_levels, byod_ms_export)
  output$byod_ms_plot <- renderPlotly({ r <- byod_ms_res(); req(r)
    tc_plotly(list(summary = r$summary, splits = r$splits), r$ylab, r$ttl, byod_ms_cols()) })

  byod_hm_obj <- resettable(reactive(input$byod_hm_run), byod_gen_layer, function() { d <- byod_layer_d(); req(d)
    byod_heatmap(d, match_feats(d, parse_genes(input$byod_hm_genes)),
                 if (identical(input$byod_hm_by, "cluster")) "Timepoint" else
                   (input$byod_hm_by %||% "Timepoint")) }, ignoreNULL = FALSE,
    module = "Heatmap",
    params = function(o) list(layer = isolate(input$byod_layer), n_features = nrow(o$mat %||% matrix(nrow = 0))))

  output$byod_hm_sample <- renderPlotly({ input$byod_hm_run
    d <- byod_layer_d(); req(d)
    isolate({
      gs <- match_feats(d, parse_genes(input$byod_hm_genes))
      gs <- intersect(gs, rownames(d$expr))
      if (length(gs) < 2) return(emptyplot("Enter at least two valid features."))
      Z <- t(scale(t(as.matrix(d$expr[gs, , drop = FALSE])))); Z[!is.finite(Z)] <- 0
      g <- byod_sample_heatmap_gg(Z, d$meta, input$byod_hm_by %||% "Timepoint",
                                  title = sprintf("%d features, row z-score", length(gs)),
                                  value_name = "z", max_rows = 120L)
      ggplotly(g) %>% ply()
    }) })
  byod_hm_export <- reactive({ o <- byod_hm_obj(); req(o$mat)
    list(summary = o$mat, gg = gg_heatmap, ttl = o$ttl, prefix = "heatmap", diverging = TRUE,
         result = list(values = o$mat)) })
  byod_hm_spec <- apprServerHeat("byod_hm_appr", byod_hm_export)
  output$byod_heatmap <- renderPlotly({ o <- byod_hm_obj(); req(o)
    if (is.null(o$mat)) return(o$plot)
    plot_ly(x = colnames(o$mat), y = rownames(o$mat), z = o$mat, type = "heatmap",
        colorscale = byod_hm_spec()$colorscale, zmid = 0, hovertemplate = "%{y} @ %{x}<br>z=%{z:.2f}<extra></extra>") %>%
      ply(xaxis = list(title = o$by), yaxis = list(title = ""), title = list(text = o$ttl, x = .02, font = list(size = 14))) })

  byod_as_obj <- reactive({ d <- byod_layer_d(); req(d, input$byod_as_feat, input$byod_as_var)
    byod_assoc(d, input$byod_as_feat, input$byod_as_var, input$byod_as_tp %||% "All") })
  byod_as_levels <- reactive(byod_as_obj()$splits %||% character(0))
  byod_as_export <- reactive({ o <- byod_as_obj()
    list(summary = o$summary, ylab = o$ylab, ttl = o$ttl, splitby = NULL, gg = o$gg, csv = o$csv, prefix = "association") })
  byod_as_cols <- apprServer("byod_as_appr", byod_as_levels, byod_as_export)
  output$byod_as_plot <- renderPlotly({ d <- byod_layer_d(); req(d, input$byod_as_feat, input$byod_as_var)
    byod_assoc(d, input$byod_as_feat, input$byod_as_var, input$byod_as_tp %||% "All", cols = byod_as_cols())$plot })

  output$byod_cc_has <- renderText({ s <- byod_study(); if (!is.null(s) && length(c(s$clinical_num, s$clinical_cat))) "yes" else "no" })
  outputOptions(output, "byod_cc_has", suspendWhenHidden = FALSE)
  byod_cc_res <- resettable(reactive(input$byod_cc_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d, input$byod_cc_var)
    withProgress(message = "Clinical correlation", detail = "per-feature across donors...", value = .5,
      byod_clin_corr(d, input$byod_cc_var, input$byod_cc_tp %||% "All", input$byod_cc_method %||% "spearman")) },
    module = "Clinical correlates",
    params = function(o) list(var = isolate(input$byod_cc_var),
      timepoint = isolate(input$byod_cc_tp) %||% "All", method = o$method %||% "moderated t-test",
      .engine = if (isTRUE(o$isnum)) "stats" else "limma",
      unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  output$byod_cc_note <- renderUI({ o <- byod_cc_res()
    if (is.null(o)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Pick a clinical variable + timepoint, then Correlate. Numeric variables use per-feature Spearman/Pearson across donors; categorical (e.g. DBD/DCD, HMP/SCS) use a donor-level moderated t-test. All statistics are at donor resolution (one value per kidney)."))
    if (!isTRUE(o$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", o$msg))
    nsig <- sum(is.finite(o$tab$adj.P.Val) & o$tab$adj.P.Val < (input$byod_cc_fdr %||% .05))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;",
      sprintf("%s -- %d features tested, %d significant at FDR < %s.", o$ttl, nrow(o$tab), nsig, input$byod_cc_fdr %||% .05)) })
  byod_cc_levels <- reactive({ byod_cc_res(); VOLCANO_LEVELS })
  byod_cc_export <- reactive({ o <- byod_cc_res(); req(isTRUE(o$ok))
    v <- byod_clin_volcano(o, input$byod_cc_fdr %||% .05, input$byod_cc_eff %||% 0)
    list(summary = v$summary, ylab = v$ylab, ttl = o$ttl, splitby = NULL, csv = o$tab, gg = gg_volcano,
         prefix = paste0("clin_corr_", gsub("[^A-Za-z0-9]+", "_", o$var))) })
  byod_cc_cols <- apprServer("byod_cc_appr", byod_cc_levels, byod_cc_export)
  output$byod_cc_volcano <- renderPlotly({ o <- byod_cc_res(); req(isTRUE(o$ok))
    byod_clin_volcano(o, input$byod_cc_fdr %||% .05, input$byod_cc_eff %||% 0, cols = byod_cc_cols(), n_label = 8)$plot })
  output$byod_cc_table <- renderDT({ o <- byod_cc_res(); req(isTRUE(o$ok))
    dt <- head(o$tab, 500); dt$logFC <- round(dt$logFC, 3); dt$AveExpr <- round(dt$AveExpr, 2)
    dt$t <- round(dt$t, 2); dt$P.Value <- signif(dt$P.Value, 3); dt$adj.P.Val <- signif(dt$adj.P.Val, 3)
    if (isTRUE(o$isnum)) names(dt)[names(dt) == "logFC"] <- if (o$method == "spearman") "rho" else "r"
    datatable(dt, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE)) })
  outputOptions(output, "byod_cc_table", suspendWhenHidden = FALSE)
  output$byod_cc_dl <- downloadHandler("clinical_correlates.csv", function(f) { o <- byod_cc_res(); req(isTRUE(o$ok)); write.csv(o$tab, f, row.names = FALSE) })
  observeEvent(byod_cc_res(), { o <- byod_cc_res(); req(isTRUE(o$ok))
    updateSelectizeInput(session, "byod_cc_feat", choices = o$tab$Feature, selected = o$tab$Feature[1], server = TRUE) })
  byod_cc_sc_obj <- reactive({ o <- byod_cc_res(); req(isTRUE(o$ok), input$byod_cc_feat)
    byod_clin_scatter(byod_layer_d(), input$byod_cc_feat, o$var, o$at_tp, o$method %||% "spearman") })
  byod_cc_sc_levels <- reactive(byod_cc_sc_obj()$splits %||% character(0))
  byod_cc_sc_export <- reactive({ x <- byod_cc_sc_obj()
    list(summary = x$summary, ylab = x$ylab, ttl = x$ttl, splitby = NULL, gg = x$gg, csv = x$csv, prefix = "clin_feature") })
  byod_cc_sc_cols <- apprServer("byod_cc_sc_appr", byod_cc_sc_levels, byod_cc_sc_export)
  output$byod_cc_scatter <- renderPlotly({ o <- byod_cc_res(); req(isTRUE(o$ok), input$byod_cc_feat)
    byod_clin_scatter(byod_layer_d(), input$byod_cc_feat, o$var, o$at_tp, o$method %||% "spearman", cols = byod_cc_sc_cols())$plot })
  byod_cc_en_obj <- reactive({ o <- byod_cc_res(); req(isTRUE(o$ok), input$byod_cc_coll)
    byod_clin_enrich(o, byod_layer_d(), input$byod_cc_coll, input$byod_cc_emode %||% "gsea", input$byod_cc_fdr %||% .05) })
  byod_cc_en_levels <- reactive(byod_cc_en_obj()$splits %||% character(0))
  byod_cc_en_export <- reactive({ x <- byod_cc_en_obj()
    list(summary = x$summary, ylab = x$ylab, ttl = x$ttl, splitby = NULL,
         gg = if ((input$byod_cc_emode %||% "gsea") == "ora") gg_ora else gg_gsea, csv = x$summary, prefix = "clin_enrichment") })
  byod_cc_en_cols <- apprServer("byod_cc_en_appr", byod_cc_en_levels, byod_cc_en_export)
  output$byod_cc_enrich <- renderPlotly({ o <- byod_cc_res(); req(isTRUE(o$ok), input$byod_cc_coll)
    byod_clin_enrich(o, byod_layer_d(), input$byod_cc_coll, input$byod_cc_emode %||% "gsea", input$byod_cc_fdr %||% .05, cols = byod_cc_en_cols())$plot })

  output$byod_gti_has <- renderText({ d <- byod_layer_d(); if (!is.null(d) && length(byod_gti_cands(d))) "yes" else "no" })
  outputOptions(output, "byod_gti_has", suspendWhenHidden = FALSE)
  observeEvent(byod_layer_d(), { d <- byod_layer_d(); req(d); cands <- byod_gti_cands(d)
    if (length(cands)) updateSelectInput(session, "byod_gti_var", choices = cands, selected = cands[1]) })
  byod_gti_res <- resettable(reactive(input$byod_gti_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d, input$byod_gti_var)
    withProgress(message = "Group x time interaction", detail = "limma + subject blocking...", value = .4,
      byod_gti(d, input$byod_gti_var)) },
    module = "Group x time interaction",
    params = function(o) list(var = isolate(input$byod_gti_var), paired = isTRUE(o$paired), unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  output$byod_gti_note <- renderUI({ o <- byod_gti_res()
    if (is.null(o)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Pick a 2-level group (e.g. warm-ischemia status) and Run. Tests which features' TIME-COURSE differs between the two groups (limma group x time interaction, subject-blocked)."))
    if (!isTRUE(o$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", o$msg))
    nsig <- sum(is.finite(o$tab$adj.P.Val) & o$tab$adj.P.Val < .05)
    div(class = "alert alert-success py-1", style = "font-size:.8rem;", sprintf("%s -- %d significant interaction features at FDR < 0.05.", o$ttl, nsig)) })
  output$byod_gti_table <- renderDT({ o <- byod_gti_res(); req(isTRUE(o$ok))
    dt <- head(o$tab, 500); dt[[o$stat]] <- round(dt[[o$stat]], 2); dt$P.Value <- signif(dt$P.Value, 3); dt$adj.P.Val <- signif(dt$adj.P.Val, 3)
    datatable(dt, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })
  outputOptions(output, "byod_gti_table", suspendWhenHidden = FALSE)
  observeEvent(byod_gti_res(), { o <- byod_gti_res(); req(isTRUE(o$ok))
    updateSelectizeInput(session, "byod_gti_gene", choices = o$tab$Feature, selected = o$tab$Feature[1], server = TRUE) })
  byod_gti_tc_obj <- reactive({ o <- byod_gti_res(); req(isTRUE(o$ok), input$byod_gti_gene)
    byod_gene_tc(byod_layer_d(), input$byod_gti_gene, o$group_var) })
  byod_gti_tc_levels <- reactive(byod_gti_tc_obj()$splits %||% character(0))
  byod_gti_tc_export <- reactive({ x <- byod_gti_tc_obj()
    list(summary = x$summary, ylab = x$ylab, ttl = x$ttl, splitby = isolate(byod_gti_res()$group_var), gg = tc_gg, csv = x$summary, prefix = "gti_trajectory") })
  byod_gti_tc_cols <- apprServer("byod_gti_tc_appr", byod_gti_tc_levels, byod_gti_tc_export)
  output$byod_gti_traj <- renderPlotly({ o <- byod_gti_res(); req(isTRUE(o$ok), input$byod_gti_gene)
    byod_gene_tc(byod_layer_d(), input$byod_gti_gene, o$group_var, cols = byod_gti_tc_cols())$plot })
  byod_gti_en_obj <- reactive({ o <- byod_gti_res(); req(isTRUE(o$ok), input$byod_gti_coll)
    byod_gti_enrich(o, byod_layer_d(), input$byod_gti_coll) })
  byod_gti_en_levels <- reactive(byod_gti_en_obj()$splits %||% character(0))
  byod_gti_en_export <- reactive({ x <- byod_gti_en_obj()
    list(summary = x$summary, ylab = x$ylab, ttl = x$ttl, splitby = NULL, gg = gg_ora, csv = x$summary, prefix = "gti_enrichment") })
  byod_gti_en_cols <- apprServer("byod_gti_en_appr", byod_gti_en_levels, byod_gti_en_export)
  output$byod_gti_enrich <- renderPlotly({ o <- byod_gti_res(); req(isTRUE(o$ok), input$byod_gti_coll)
    byod_gti_enrich(o, byod_layer_d(), input$byod_gti_coll, cols = byod_gti_en_cols())$plot })

  output$byod_tp3_has <- renderText({ d <- byod_layer_d(); if (!is.null(d) && length(d$tp_levels) >= 3) "yes" else "no" })
  outputOptions(output, "byod_tp3_has", suspendWhenHidden = FALSE)
  byod_tp_res <- resettable(reactive(input$byod_tp_run), byod_gen_layer, function() {
    d <- byod_layer_d(); req(d)
    withProgress(message = "Temporal programs", detail = "limma over time + fuzzy clustering...", value = .4,
      byod_temporal(d, input$byod_tp_mode %||% "factor", input$byod_tp_fdr %||% .05, input$byod_tp_k %||% 6,
                    core_cut = input$byod_tp_core %||% 0)) },
    module = "Temporal programs",
    params = function(o) list(mode = isolate(input$byod_tp_mode), k = isolate(input$byod_tp_k), unit = o$unit %||% "-", n_units = o$n_units %||% NA))
  output$byod_tp_note <- renderUI({ o <- byod_tp_res()
    if (is.null(o)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Detects features that change across timepoints (limma, subject-blocked) and clusters their z-scored profiles into temporal programs. Run to compute."))
    if (!isTRUE(o$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", o$msg))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;", paste0(o$ttl, ".", if (!is.null(o$note)) paste0(" ", o$note) else "")) })
  byod_tp_prof_levels <- reactive({ o <- byod_tp_res(); if (is.null(o) || !isTRUE(o$ok)) character(0) else byod_temporal_profiles(o)$splits })
  byod_tp_prof_export <- reactive({ o <- byod_tp_res(); req(isTRUE(o$ok)); pr <- byod_temporal_profiles(o)
    list(summary = pr$summary, ylab = pr$ylab, ttl = pr$ttl, splitby = "cluster", gg = tc_gg, csv = pr$summary, prefix = "temporal_programs") })
  byod_tp_prof_cols <- apprServer("byod_tp_prof_appr", byod_tp_prof_levels, byod_tp_prof_export)
  output$byod_tp_profiles <- renderPlotly({ o <- byod_tp_res(); req(isTRUE(o$ok))
    byod_temporal_profiles(o, cols = byod_tp_prof_cols())$plot })
  byod_tp_heat_export <- reactive({ o <- byod_tp_res(); req(isTRUE(o$ok)); hm <- byod_temporal_heat(o)
    list(summary = hm$mat, gg = gg_heatmap, ttl = hm$ttl, prefix = "temporal_heat", diverging = TRUE,
         result = list(feature_cluster = o$clusters, cluster_profiles = hm$mat,
                       standardised_values = o$z)) })
  byod_tp_heat_spec <- apprServerHeat("byod_tp_heat_appr", byod_tp_heat_export)
  output$byod_tp_heat <- renderPlotly({ o <- byod_tp_res(); req(isTRUE(o$ok))
    byod_temporal_heat(o, colorscale = byod_tp_heat_spec()$colorscale)$plot })
  observeEvent(byod_tp_res(), { o <- byod_tp_res(); req(isTRUE(o$ok))
    updateSelectInput(session, "byod_tp_clu", choices = sort(unique(o$clusters$cluster)), selected = sort(unique(o$clusters$cluster))[1]) })
  byod_tp_en_obj <- reactive({ o <- byod_tp_res(); req(isTRUE(o$ok), input$byod_tp_clu, input$byod_tp_coll)
    byod_temporal_enrich(o, byod_layer_d(), input$byod_tp_clu, input$byod_tp_coll) })
  byod_tp_en_levels <- reactive(byod_tp_en_obj()$splits %||% character(0))
  byod_tp_en_export <- reactive({ x <- byod_tp_en_obj()
    list(summary = x$summary, ylab = x$ylab, ttl = x$ttl, splitby = NULL, gg = gg_ora, csv = x$summary, prefix = "temporal_cluster_enrich") })
  byod_tp_en_cols <- apprServer("byod_tp_en_appr", byod_tp_en_levels, byod_tp_en_export)
  output$byod_tp_enrich <- renderPlotly({ o <- byod_tp_res(); req(isTRUE(o$ok), input$byod_tp_clu, input$byod_tp_coll)
    byod_temporal_enrich(o, byod_layer_d(), input$byod_tp_clu, input$byod_tp_coll, cols = byod_tp_en_cols())$plot })
  output$byod_tp_dl <- downloadHandler("temporal_programs.csv", function(f) { o <- byod_tp_res(); req(isTRUE(o$ok)); write.csv(o$tab, f, row.names = FALSE) })

  byod_de_top <- reactive({
    o <- byod_de_res(); req(is.data.frame(o), nrow(o) > 0)
    fdr <- input$byod_de_fdr %||% 0.05
    lfc <- input$byod_de_lfc %||% 0
    dir <- input$byod_de_follow_dir %||% "either"
    n   <- max(1L, as.integer(input$byod_de_follow_n %||% 30L))
    t <- o[is.finite(o$adj.P.Val) & o$adj.P.Val < fdr & abs(o$logFC) >= lfc, , drop = FALSE]
    if (identical(dir, "up"))   t <- t[t$logFC > 0, , drop = FALSE]
    if (identical(dir, "down")) t <- t[t$logFC < 0, , drop = FALSE]
    t <- t[order(t$adj.P.Val), , drop = FALSE]
    list(features = head(t$Feature, n), n_available = nrow(t), n_asked = n,
         fdr = fdr, lfc = lfc, direction = dir, label = byod_de_lbl())
  })

  output$byod_de_follow_note <- renderUI({
    tt <- tryCatch(byod_de_top(), error = function(e) NULL)
    if (is.null(tt)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Run a contrast above, then take its strongest features straight into a heatmap or an over-representation test. The list is taken from this result, so it cannot drift from the volcano you are looking at."))
    if (!length(tt$features)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;",
      sprintf("No feature passes FDR < %g and |log2FC| >= %g%s, so there is nothing to carry forward.",
              tt$fdr, tt$lfc,
              if (tt$direction == "either") "" else paste(" in the", tt$direction, "direction"))))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;",
        sprintf("%d of %s feature%s passing FDR < %g and |log2FC| >= %g%s, ranked by adjusted p.",
                length(tt$features), format(tt$n_available, big.mark = ","),
                if (tt$n_available == 1) "" else "s", tt$fdr, tt$lfc,
                if (tt$direction == "either") "" else paste(",", tt$direction, "only")))
  })

  byod_de_hm_obj <- resettable(reactive(input$byod_de_hm_run), byod_gen_layer, function() {
      d <- byod_layer_d(); tt <- byod_de_top(); req(d, length(tt$features) >= 2)
      o <- byod_heatmap(d, tt$features, input$byod_de_hm_by %||% "Timepoint")
      o$ttl <- sprintf("%d features from %s (row z-score)", length(tt$features), tt$label)
      o },
    module = "Heatmap",
    params = function(o) list(layer = isolate(input$byod_layer), source = "differential result",
                              n_features = nrow(o$mat %||% matrix(nrow = 0)),
                              ranked_by = "adjusted p", direction = isolate(input$byod_de_follow_dir)))
  output$byod_de_hm_plot <- renderPlotly({ o <- byod_de_hm_obj(); req(o)
    if (is.null(o$mat)) return(o$plot)
    ggplotly(gg_heatmap(o$mat, byod_de_hm_spec()$stops, o$ttl, TRUE)) })
  byod_de_hm_export <- reactive({ o <- byod_de_hm_obj(); req(o$mat)
    list(summary = o$mat, gg = gg_heatmap, ttl = o$ttl, diverging = TRUE,
         result = list(values = o$mat, differential = byod_de_res()),
         prefix = "differential_top_heatmap") })
  byod_de_hm_spec <- apprServerHeat("byod_de_hm_appr", byod_de_hm_export)

  byod_de_ora_res <- resettable(reactive(input$byod_de_ora_run), byod_gen_layer, function() {
      d <- byod_layer_d(); tt <- byod_de_top(); req(d, length(tt$features) >= 3)
      coll <- input$byod_de_ora_coll %||% names(GENESETS)[1]
      withProgress(message = "Over-representation...", value = .5,
        list(res = run_ora(tt$features, coll, rownames(d$expr), d = d),
             lab = sprintf("%s - top %d of %s by adjusted p", coll, length(tt$features),
                           format(tt$n_available, big.mark = ",")),
             n_in = length(tt$features), collection = coll,
             universe = nrow(d$expr))) },
    module = "Enrichment",
    params = function(o) list(layer = isolate(input$byod_layer), mode = "ora",
                              source = "differential result", collection = o$collection,
                              n_features = o$n_in, universe = o$universe))

  byod_de_ora_cut <- reactive({
    list(cut = max(1e-12, min(1, as.numeric(input$byod_de_ora_cut %||% 0.05))),
         alpha = max(0.05, min(1, as.numeric(input$byod_de_ora_alpha %||% 0.35))),
         fade = isTRUE(input$byod_de_ora_fade)) })

  byod_de_ora_obj <- reactive({ r <- byod_de_ora_res(); req(r); cf <- byod_de_ora_cut()
    plot_ora(r$res, as.integer(input$byod_de_ora_ntop %||% 20L), r$lab,
             sig_cut = cf$cut, faded_alpha = cf$alpha, fade = cf$fade) })
  output$byod_de_ora_plot <- renderPlotly({ r <- byod_de_ora_res(); req(r)
    cf <- byod_de_ora_cut()
    plot_ora(r$res, as.integer(input$byod_de_ora_ntop %||% 20L), r$lab,
             cols = byod_de_ora_cols(), sig_cut = cf$cut,
             faded_alpha = cf$alpha, fade = cf$fade)$plot })
  output$byod_de_ora_table <- renderDT({ r <- byod_de_ora_res(); req(r)
    validate(need(nrow(r$res) > 0, "No set is over-represented in this list."))
    dt <- r$res; dt$p <- signif(dt$p, 3); dt$padj <- signif(dt$padj, 3)
    datatable(dt, rownames = FALSE, options = list(pageLength = 8, scrollX = TRUE)) })
  byod_de_ora_levels <- reactive(byod_de_ora_obj()$splits %||% character(0))
  byod_de_ora_export <- reactive({ o <- byod_de_ora_obj(); req(o); cf <- byod_de_ora_cut()
    list(summary = o$summary, ylab = o$ylab, ttl = o$ttl, splitby = NULL,
         csv = o$summary, prefix = "differential_top_ora",

         gg = function(sm, yl, tt, co) gg_ora(sm, yl, tt, co,
                sig_cut = cf$cut, faded_alpha = cf$alpha, fade = cf$fade)) })
  byod_de_ora_cols <- apprServer("byod_de_ora_appr", byod_de_ora_levels, byod_de_ora_export)

  output$byod_de_top_dl <- downloadHandler("differential_top_features.csv", function(f) {
    tt <- byod_de_top()
    utils::write.csv(data.frame(rank = seq_along(tt$features), Feature = tt$features,
                                contrast = tt$label, fdr_threshold = tt$fdr,
                                lfc_threshold = tt$lfc, direction = tt$direction),
                     f, row.names = FALSE) })

  byod_gti_top <- reactive({
    o <- byod_gti_res(); req(isTRUE(o$ok))
    fdr <- input$byod_gti_fdr %||% 0.05
    n   <- max(2L, as.integer(input$byod_gti_follow_n %||% 30L))
    t <- o$tab[is.finite(o$tab$adj.P.Val) & o$tab$adj.P.Val < fdr, , drop = FALSE]
    t <- t[order(t$adj.P.Val), , drop = FALSE]
    list(features = head(t$Feature, n), n_available = nrow(t), fdr = fdr,
         group_var = o$group_var, unit = o$unit, n_units = o$n_units,
         paired = isTRUE(o$paired))
  })

  output$byod_gti_follow_note <- renderUI({
    tt <- tryCatch(byod_gti_top(), error = function(e) NULL)
    if (is.null(tt)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Run the interaction above, then look at the features themselves. Arranged group by timepoint, the two arms sit side by side at every point."))
    if (!length(tt$features)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;",
      sprintf("No feature reaches FDR < %g, so there is nothing to carry forward.", tt$fdr)))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;",
        sprintf("%d of %s feature%s with an interaction at FDR < %g, ranked by adjusted p, from %s %s.",
                length(tt$features), format(tt$n_available, big.mark = ","),
                if (tt$n_available == 1) "" else "s", tt$fdr,
                tt$n_units %||% NA, tt$unit %||% "units"))
  })

  byod_gti_hm_obj <- resettable(reactive(input$byod_gti_hm_run), byod_gen_layer, function() {
      d <- byod_layer_d(); tt <- byod_gti_top(); req(d, length(tt$features) >= 2)
      by <- input$byod_gti_hm_by %||% "group_time"
      if (identical(by, "group_time")) {

        tp <- factor(d$meta$Timepoint, levels = d$tp_levels)
        gf <- factor(as.character(d$meta[[tt$group_var]]))

        lab <- paste(as.character(gf), as.character(tp), sep = " @ ")
        lv  <- as.vector(outer(levels(gf), levels(droplevels(tp)),
                               function(a, b) paste(a, b, sep = " @ ")))
        d2 <- d; d2$meta[[".gti_group_time"]] <- factor(lab, levels = lv[lv %in% lab])
        o <- byod_heatmap(d2, tt$features, ".gti_group_time")
      } else o <- byod_heatmap(d, tt$features, "Timepoint")
      o$ttl <- sprintf("%d interaction features - %s x time, %s, n = %s %s (row z-score)",
                       length(tt$features), tt$group_var,
                       if (tt$paired) "within-unit paired" else "subject-blocked",
                       tt$n_units %||% NA, tt$unit %||% "units")
      o },
    module = "Heatmap",
    params = function(o) list(layer = isolate(input$byod_layer),
                              source = "group x time interaction",
                              n_features = nrow(o$mat %||% matrix(nrow = 0)),
                              columns = isolate(input$byod_gti_hm_by), ranked_by = "adjusted p"))
  output$byod_gti_hm_plot <- renderPlotly({ o <- byod_gti_hm_obj(); req(o)
    if (is.null(o$mat)) return(o$plot)
    ggplotly(gg_heatmap(o$mat, byod_gti_hm_spec()$stops, o$ttl, TRUE)) })
  byod_gti_hm_export <- reactive({ o <- byod_gti_hm_obj(); req(o$mat)
    list(summary = o$mat, gg = gg_heatmap, ttl = o$ttl, diverging = TRUE,
         result = { g <- byod_gti_res()
                    list(values = o$mat, interaction = if (isTRUE(g$ok)) g$tab else NULL) },
         prefix = "interaction_top_heatmap") })
  byod_gti_hm_spec <- apprServerHeat("byod_gti_hm_appr", byod_gti_hm_export)

  output$byod_gti_top_dl <- downloadHandler("interaction_top_features.csv", function(f) {
    tt <- byod_gti_top()
    utils::write.csv(data.frame(rank = seq_along(tt$features), Feature = tt$features,
                                group_var = tt$group_var, fdr_threshold = tt$fdr,
                                unit = tt$unit, n_units = tt$n_units),
                     f, row.names = FALSE) })

  byod_qc_export <- reactive({ d <- byod_layer_d(); req(d)
    list(summary = data.frame(), ylab = d$unit, ttl = paste(d$omic %||% "layer", "- per-sample distribution"),
         result = list(per_profile = {
           q <- t(apply(d$expr, 2, function(z) stats::quantile(z, c(0, .25, .5, .75, 1), na.rm = TRUE)))
           colnames(q) <- c("min", "q25", "median", "q75", "max")
           data.frame(sample = colnames(d$expr), as.data.frame(q),
                      mean = colMeans(d$expr, na.rm = TRUE),
                      n_missing = colSums(is.na(d$expr)),
                      d$meta[colnames(d$expr), , drop = FALSE],
                      check.names = FALSE, stringsAsFactors = FALSE) }),
         splitby = "Timepoint", prefix = "qc_distribution",
         gg = function(s, yl, tt, co) byod_dist_gg(byod_layer_d(), co)) })
  byod_qc_levels <- reactive({ d <- byod_layer_d()
    if (is.null(d)) character(0) else levels(droplevels(factor(d$meta$Timepoint, levels = d$tp_levels))) })
  byod_qc_cols <- apprServer("byod_qc_appr", byod_qc_levels, byod_qc_export)

  byod_gti_ov <- reactive({ o <- byod_gti_res(); req(isTRUE(o$ok))
    byod_gti_overview(o, byod_layer_d(), fdr = input$byod_gti_fdr %||% 0.05) })
  output$byod_gti_overview <- renderPlotly({ ov <- byod_gti_ov(); req(isTRUE(ov$ok))
    ggplotly(byod_gti_overview_gg(ov, byod_gti_ov_cols())) })
  byod_gti_ov_levels <- reactive({ ov <- byod_gti_ov()
    if (!isTRUE(ov$ok)) character(0) else unique(as.character(ov$summary$sig)) })
  byod_gti_ov_export <- reactive({ ov <- byod_gti_ov(); req(isTRUE(ov$ok))
    list(summary = ov$summary, ylab = ov$ylab, ttl = ov$ttl, splitby = NULL,
         csv = ov$summary, prefix = "gti_overview",
         gg = function(s, yl, tt, co) byod_gti_overview_gg(byod_gti_ov(), co)) })
  byod_gti_ov_cols <- apprServer("byod_gti_ov_appr", byod_gti_ov_levels, byod_gti_ov_export)

  byod_en_curve_obj <- reactive({ r <- byod_en_res(); req(r$kind == "gsea", input$byod_en_curve_path)
    d <- byod_layer_d(); req(d)
    tt <- withProgress(message = "Ranking for the curve...",
      byod_run_limma(d, input$byod_en_var, input$byod_en_g1, input$byod_en_g2,
                     byod_subj(), isTRUE(input$byod_en_paired)))
    byod_gsea_curve(setNames(tt$t, tt$Feature), input$byod_en_coll, input$byod_en_curve_path, d = byod_layer_d()) })
  output$byod_en_curve <- renderPlotly({ o <- byod_en_curve_obj(); req(isTRUE(o$ok))
    ggplotly(byod_gsea_curve_gg(o)) })
  output$byod_en_curve_pick <- renderUI({ r <- byod_en_res(); req(r$kind == "gsea")
    ch <- head(r$res$pathway[order(r$res$padj)], 60)
    selectizeInput("byod_en_curve_path", "Set to draw", choices = ch, selected = ch[1]) })
  byod_en_curve_export <- reactive({ o <- byod_en_curve_obj(); req(isTRUE(o$ok))
    list(summary = o$curve, ylab = "running enrichment score", ttl = o$ttl, splitby = NULL,
         csv = o$curve, prefix = "gsea_curve",
         gg = function(s, yl, tt, co) byod_gsea_curve_gg(byod_en_curve_obj())) })
  byod_en_curve_cols <- apprServer("byod_en_curve_appr", reactive(character(0)), byod_en_curve_export)

  byod_conc_hist_export <- reactive({ o <- byod_conc_res(); req(o)
    list(summary = o$tab, ylab = "genes", ttl = o$ttl %||% "transcript-protein concordance",
         splitby = NULL, csv = o$tab, prefix = "concordance_histogram",
         gg = function(s, yl, tt, co) byod_conc_hist_gg(byod_conc_res())) })
  byod_conc_hist_cols <- apprServer("byod_conc_hist_appr",
    reactive(names(BYOD_CONC_COLS)), byod_conc_hist_export)

  byod_pconc_res <- resettable(reactive(input$byod_pconc_run), byod_gen_study, function() {
      s <- byod_study(); req(s, length(s$layer_names) >= 2)
      l1 <- s$layer_names[1]; l2 <- s$layer_names[2]
      coll <- input$byod_pconc_coll %||% "Hallmark"
      validate(need(annotation_compatible(s$layers[[l1]], s$layers[[l2]], coll),
        "The layers need compatible gene-set annotation. Custom GMT requires identical set definitions."))
      withProgress(message = "Scoring pathways in both layers...", value = .4, {
        a <- byod_pathway_activity(s$layers[[l1]], coll, "GSVA", "Timepoint")
        b <- byod_pathway_activity(s$layers[[l2]], coll, "GSVA", "Timepoint")
        byod_pathway_concordance(a, b, l1, l2) }) },
    module = "Concordance",
    params = function(o) list(level = "pathway", n_pathways = o$n_pathways %||% NA,
                              n_samples = o$n_samples %||% NA))
  output$byod_pconc_note <- renderUI({ o <- byod_pconc_res()
    if (is.null(o)) return(div(class = "alert alert-secondary py-2", style = "font-size:.85rem;",
      "Gene-level agreement answers whether one transcript predicts its protein. This asks the same of a whole pathway, which is the level at which the two layers are usually compared."))
    if (!isTRUE(o$ok)) return(div(class = "alert alert-warning py-2", style = "font-size:.85rem;", o$msg))
    div(class = "alert alert-success py-1", style = "font-size:.8rem;", o$ttl) })
  output$byod_pconc_plot <- renderPlotly({ o <- byod_pconc_res(); req(isTRUE(o$ok))
    ggplotly(byod_pathway_concordance_gg(o)) })
  byod_pconc_export <- reactive({ o <- byod_pconc_res(); req(isTRUE(o$ok))
    list(summary = o$tab, ylab = "pathways", ttl = o$ttl, splitby = NULL, csv = o$tab,
         prefix = "pathway_concordance",
         gg = function(s, yl, tt, co) byod_pathway_concordance_gg(byod_pconc_res())) })
  byod_pconc_cols <- apprServer("byod_pconc_appr", reactive(character(0)), byod_pconc_export)

  output$byod_diablo_sample <- renderPlotly({ o <- byod_diablo_res(); req(o)
    ggplotly(byod_diablo_sample_gg(o)) })
  byod_diablo_sample_export <- reactive({ o <- byod_diablo_res(); req(o)
    list(summary = data.frame(), ylab = "consensus component 2",
         ttl = "Multiblock projection - descriptive fit, not cross-validated",
         splitby = isolate(input$byod_diablo_out), prefix = "diablo_projection",
         gg = function(s, yl, tt, co) byod_diablo_sample_gg(byod_diablo_res())) })
  byod_diablo_sample_levels <- reactive({ o <- byod_diablo_res()
    if (is.null(o)) character(0) else unique(as.character(o$Y %||% character(0))) })
  byod_diablo_sample_cols <- apprServer("byod_diablo_appr",
    byod_diablo_sample_levels, byod_diablo_sample_export)

  byod_diablo_load <- reactive({ o <- byod_diablo_res(); req(o)
    byod_diablo_loadings(o, comp = as.integer(input$byod_diablo_comp %||% 1L),
                         top_n = as.integer(input$byod_diablo_ntop %||% 20L)) })
  output$byod_diablo_loadings <- renderPlotly({ o <- byod_diablo_load(); req(isTRUE(o$ok))
    ggplotly(byod_diablo_loadings_gg(o)) })
  byod_diablo_load_export <- reactive({ o <- byod_diablo_load(); req(isTRUE(o$ok))
    list(summary = o$tab, ylab = "loading", ttl = o$ttl, splitby = NULL, csv = o$tab,
         prefix = "diablo_loadings",
         gg = function(s, yl, tt, co) byod_diablo_loadings_gg(byod_diablo_load())) })
  byod_diablo_load_cols <- apprServer("byod_diablo_load_appr",
    reactive(character(0)), byod_diablo_load_export)

  output$byod_diablo_cv_plot <- renderPlotly({ r <- byod_diablo_cv_res(); req(isTRUE(r$ok))
    ggplotly(byod_cv_confusion_gg(r)) })
  byod_diablo_cv_export <- reactive({ r <- byod_diablo_cv_res(); req(isTRUE(r$ok))
    list(summary = as.data.frame(as.table(as.matrix(r$confusion))),
         ylab = "observed", ttl = "cross-validated confusion matrix", splitby = NULL,
         csv = as.data.frame(as.table(as.matrix(r$confusion))), prefix = "diablo_cv_confusion",
         gg = function(s, yl, tt, co) byod_cv_confusion_gg(byod_diablo_cv_res())) })
  byod_diablo_cv_cols <- apprServer("byod_diablo_cv_appr",
    reactive(character(0)), byod_diablo_cv_export)

  byod_vs_export <- reactive({ r <- byod_vs_ok(); req(r)
    list(summary = r$tab, ylab = "features", ttl = r$ttl, splitby = NULL, csv = r$tab,
         prefix = "variable_screen",
         gg = function(sm, yl, tt, co)
           byod_variable_screen_gg(byod_vs_ok(), input$byod_vs_value %||% "count")) })
  byod_vs_cols <- apprServer("byod_vs_appr", reactive(character(0)), byod_vs_export)

  byod_pxc_export <- reactive({ r <- byod_pxc_ok(); req(r)
    list(summary = r$tab, ylab = "pathways", ttl = r$ttl, splitby = NULL, csv = r$tab,
         prefix = "pathway_clinical",
         gg = function(sm, yl, tt, co) byod_pathway_clinical_gg(byod_pxc_ok())) })
  byod_pxc_cols <- apprServer("byod_pxc_appr", reactive(character(0)), byod_pxc_export)

  byod_rc_export <- reactive({ r <- byod_rc_ok(); req(r)
    list(summary = r$tab %||% r$long, ylab = "share of the layer", ttl = r$ttl,
         splitby = isolate(input$byod_rc_split), csv = r$tab %||% r$long,
         prefix = "response_composition",
         gg = function(sm, yl, tt, co) byod_response_composition_gg(
                byod_rc_ok(), as_percent = identical(input$byod_rc_unit %||% "pct", "pct"))) })
  byod_rc_cols <- apprServer("byod_rc_appr",
    reactive({ r <- byod_rc_ok(); if (is.null(r)) character(0)
               else unique(as.character(r$tab$category %||% r$long$category %||% character(0))) }),
    byod_rc_export)

  byod_tcon_export <- reactive({ r <- byod_tcon_ok(); req(r)
    list(summary = r$tab, ylab = "held-out r", ttl = r$ttl, splitby = NULL, csv = r$tab,
         prefix = "trajectory_consistency",
         gg = function(sm, yl, tt, co) byod_traj_consistency_gg(byod_tcon_ok())) })

  byod_tcon_cols <- apprServer("byod_tcon_appr", reactive(character(0)), byod_tcon_export)

  byod_mc_upset_export <- reactive({ r <- byod_mc_ok(); req(r)
    list(summary = r$upset %||% r$tab, ylab = "features", ttl = r$ttl, splitby = NULL,
         csv = r$upset %||% r$tab, prefix = "multi_timepoint_upset",
         gg = function(sm, yl, tt, co) byod_multi_upset_gg(byod_mc_ok())) })
  byod_mc_upset_cols <- apprServer("byod_mc_upset_appr", reactive(character(0)), byod_mc_upset_export)

  byod_mc_volc_export <- reactive({ r <- byod_mc_ok(); req(r)
    list(summary = r$tab, ylab = "-log10 P", ttl = r$ttl, splitby = NULL, csv = r$tab,
         prefix = "multi_timepoint_volcanoes",
         gg = function(sm, yl, tt, co) byod_multi_volcano_gg(byod_mc_ok(),
                label_n = as.integer(input$byod_mc_nlab %||% 0L))) })
  byod_mc_volc_cols <- apprServer("byod_mc_volc_appr", reactive(VOLCANO_LEVELS), byod_mc_volc_export)

  byod_bs_export <- reactive({ r <- byod_bs_ok(); req(r)
    list(summary = r$tab, ylab = "effect size", ttl = r$ttl, splitby = NULL, csv = r$tab,
         prefix = "baseline_sensitivity",
         gg = function(sm, yl, tt, co) byod_baseline_gg(byod_bs_ok())) })
  byod_bs_cols <- apprServer("byod_bs_appr",
    reactive({ r <- byod_bs_ok(); if (is.null(r)) character(0)
               else unique(as.character(r$tab$agreement %||% character(0))) }),
    byod_bs_export)

  byod_bs_ov_export <- reactive({ r <- byod_bs_ok(); req(r)
    list(summary = r$summary %||% r$tab, ylab = "agreement", ttl = r$ttl, splitby = NULL,
         csv = r$summary %||% r$tab, prefix = "baseline_overlap",
         gg = function(sm, yl, tt, co) byod_baseline_overlap_gg(byod_bs_ok())) })
  byod_bs_ov_cols <- apprServer("byod_bs_ov_appr", reactive(character(0)), byod_bs_ov_export)

  .xs_side_params <- function(o) {
    one <- function(sd, tag) {
      if (is.null(sd)) return(stats::setNames(list("-", "-", "-", NA), paste0(tag, c("_contrast", "_unit", "_pairing", "_n_units"))))
      executed <- isTRUE(sd$ok)
      c(stats::setNames(list(
        sd$text %||% "-",
        sd$unit %||% "not fitted",
        if (!executed) "not fitted"
        else if (isTRUE(sd$nested)) "nested, within-unit paired"
        else if (isTRUE(sd$paired)) "paired within subject" else "unpaired",
        sd$n_units %||% NA_integer_),
        paste0(tag, c("_contrast", "_unit", "_pairing", "_n_units"))),
        stats::setNames(list(sd$spec$paired, sd$declared_subject %||% NA_character_,
          sd$declared_nest %||% NA_character_, if (executed) "fitted" else "not fitted"),
          paste0(tag, c("_requested_pairing", "_declared_subject", "_declared_nest", "_execution"))),
        stats::setNames(list(sd$excluded_profiles %||% 0L, paste(sd$excluded_subjects, collapse = ";")),
          paste0(tag, c("_excluded_profiles", "_excluded_subjects"))))
    }
    c(one(o$sideA, "A"), one(o$sideB, "B"),
      list(model = if (isTRUE(o$sideA$ok) && isTRUE(o$sideB$ok))
        "limma moderated t, fitted separately in each study" else "see per-side execution status"))
  }

  csr_sA <- reactive({ reg <- csr_registry(); req(input$csr_A); req(reg[[input$csr_A]]); reg[[input$csr_A]]$study })
  csr_sB <- reactive({ reg <- csr_registry(); req(input$csr_B); req(reg[[input$csr_B]]); reg[[input$csr_B]]$study })

  xs_gen <- reactiveVal(0L)
  observeEvent(list(csr_idA(), csr_idB(), input$csr_layA, input$csr_layB, input$csr_1to1),
               xs_gen(isolate(xs_gen()) + 1L), ignoreInit = TRUE)

  xs_bridge <- reactive({ dA <- csr_dA(); dB <- csr_dB(); req(dA, dB)
    xs_feature_bridge(dA, dB, isTRUE(input$csr_1to1)) })

  output$xs_tp_map_ui <- renderUI({
    dA <- csr_dA(); dB <- csr_dB(); req(dA, dB)
    a <- cs_timepoints(dA); b <- cs_timepoints(dB)
    tagList(
      div(class = "small text-muted mb-1",
          sprintf("For each timepoint of %s, name the timepoint of %s it corresponds to.",
                  csr_labA(), csr_labB())),
      lapply(seq_along(a), function(i)
        selectInput(paste0("xs_tp_", i), a[i], choices = c("(none)", b),
                    selected = if (a[i] %in% b) a[i] else "(none)")))
  })
  xs_manual_map <- reactive({
    dA <- csr_dA(); req(dA); a <- cs_timepoints(dA)
    v <- vapply(seq_along(a), function(i) input[[paste0("xs_tp_", i)]] %||% "(none)", character(1))
    v[v == "(none)"] <- ""
    setNames(v, a)
  })
  xs_align <- reactive({
    dA <- csr_dA(); dB <- csr_dB(); req(dA, dB)
    m <- input$xs_tp_mode %||% "label"
    xs_align_timepoints(dA, dB, m, manual = if (identical(m, "manual")) xs_manual_map() else NULL)
  })
  output$xs_align_note <- renderUI({
    al <- xs_align()
    cls <- if (al$n >= 2) "alert-secondary" else "alert-warning"
    div(class = paste("alert py-2 px-3 small mt-2", cls),
        tags$b(sprintf("%d timepoint%s matched. ", al$n, if (al$n == 1) "" else "s")),
        al$note %||% sprintf("Matched: %s.", paste(al$label, collapse = ", ")))
  })

  xs_design_r <- reactive({ dA <- csr_dA(); dB <- csr_dB(); req(dA, dB)
    xs_design_table(csr_sA(), csr_sB(), dA, dB, csr_labA(), csr_labB(),
                    bridge = xs_bridge(), align = xs_align()) })
  xs_verdict_r <- reactive(xs_design_verdict(xs_design_r()))

  output$xs_pair_verdict <- renderUI({
    v <- tryCatch(xs_verdict_r(), error = function(e) NULL); req(v)
    div(class = paste0("alert alert-", v$type, " py-2 px-3 small mb-0 mt-2"), v$txt)
  })
  output$xs_verdict_banner <- renderUI({
    v <- xs_verdict_r(); br <- xs_bridge()
    tagList(
      div(class = paste0("alert alert-", v$type, " py-2 px-3 small"), v$txt),
      div(class = "small text-muted mb-2", br$note))
  })
  output$xs_design_tab <- renderDT({
    d <- xs_design_r()
    datatable(d, rownames = FALSE, options = list(pageLength = 12, dom = "t", ordering = FALSE)) %>%
      formatStyle("Verdict", target = "row",
                  backgroundColor = styleEqual(
                    c("blocking", "feature-level blocked", "time-contrast blocked", "check", "ok"),
                    c("#fbeaea", "#fbeaea", "#fdf5e6", "#fdf5e6", "#eef6f6")))
  })
  output$xs_design_dl <- downloadHandler(
    filename = function() sprintf("design_comparison_%s_vs_%s.csv",
                                  gsub("[^A-Za-z0-9]+", "_", csr_labA()),
                                  gsub("[^A-Za-z0-9]+", "_", csr_labB())),
    content = function(f) write.csv(xs_design_r(), f, row.names = FALSE))
  output$xs_design_plot <- renderPlot({ xs_design_gg(xs_design_r(), csr_labA(), csr_labB()) })
  xs_design_export <- reactive({ d <- xs_design_r()
    list(summary = d, csv = d, ylab = "", ttl = "design comparison", splitby = NULL,
         prefix = "design_comparison",
         gg = function(sm, yl, tt, co) xs_design_gg(isolate(xs_design_r()),
                                                    isolate(csr_labA()), isolate(csr_labB()))) })
  xs_design_cols <- apprServer("xs_design_appr", reactive(character(0)), xs_design_export)

  local({
    for (side in c("A", "B")) local({
      sd <- side
      dfun <- if (sd == "A") function() csr_dA() else function() csr_dB()
      observeEvent(list(dfun(), input[[paste0("csr_lay", sd)]]), {
        d <- dfun(); req(d)
        gv <- group_vars(d)
        cur <- isolate(input[[paste0("xs_var", sd)]])
        updateSelectInput(session, paste0("xs_var", sd), choices = gv,
                          selected = if (!is.null(cur) && cur %in% gv) cur else "Timepoint")
        updateSelectInput(session, paste0("xs_tp", sd), choices = c("All", cs_timepoints(d)),
                          selected = "All")
      })
      observeEvent(list(dfun(), input[[paste0("xs_var", sd)]]), {
        d <- dfun(); req(d)
        v <- input[[paste0("xs_var", sd)]]
        if (is.null(v) || !v %in% group_vars(d)) v <- "Timepoint"
        lv <- group_levels(d, v)

        updateSelectInput(session, paste0("xs_test", sd), choices = lv, selected = tail(lv, 1))
        updateSelectInput(session, paste0("xs_ref", sd), choices = lv, selected = head(lv, 1))
      })
    })
  })

  observeEvent(input$xs_mirror, {
    al <- xs_align(); dB <- csr_dB(); req(dB)
    if (!identical(input$xs_varA, "Timepoint"))
      return(byod_msg(list(type = "warning",
        txt = "The timepoint mapping can only fill in a contrast that is itself between timepoints. Declare study B's contrast directly.")))
    m <- setNames(al$B, al$A)
    g1 <- m[[input$xs_testA %||% ""]] %||% NA_character_
    g2 <- m[[input$xs_refA %||% ""]] %||% NA_character_
    if (is.na(g1) || is.na(g2))
      return(byod_msg(list(type = "warning", txt = sprintf(
        "The current mapping does not cover both of study A's levels (%s, %s). Declare them on the Side by side tab first.",
        input$xs_testA, input$xs_refA))))
    updateSelectInput(session, "xs_varB", selected = "Timepoint")
    updateSelectInput(session, "xs_testB", choices = cs_timepoints(dB), selected = g1)
    updateSelectInput(session, "xs_refB", choices = cs_timepoints(dB), selected = g2)
    byod_msg(list(type = "success", txt = sprintf(
      "Study B's contrast set to %s vs %s, from the declared mapping.", g1, g2)))
  })

  local({ for (side in c("A", "B")) local({
    sd <- side
    output[[paste0("xs_dir", sd)]] <- renderUI({
      te <- input[[paste0("xs_test", sd)]]; rf <- input[[paste0("xs_ref", sd)]]
      if (is.null(te) || is.null(rf) || identical(te, rf)) return(NULL)
      div(class = "small text-muted mb-1",
          sprintf("A positive effect estimate means higher in %s than in %s, on the declared analysis scale.", te, rf))
    })
  }) })

  xs_specA <- reactive(xs_spec(input$xs_varA, input$xs_testA, input$xs_refA,
                               input$xs_tpA %||% "All", isTRUE(input$xs_paired)))
  xs_specB <- reactive(xs_spec(input$xs_varB, input$xs_testB, input$xs_refB,
                               input$xs_tpB %||% "All", isTRUE(input$xs_paired)))

  xs_de_res <- resettable(reactive(input$xs_de_run), xs_gen, function() {
    dA <- csr_dA(); dB <- csr_dB()
    req(dA, dB, input$xs_testA, input$xs_refA, input$xs_testB, input$xs_refB)
    validate(
      need(!identical(input$xs_testA, input$xs_refA),
           sprintf("%s's two levels are the same; a contrast needs two different levels.", csr_labA())),
      need(!identical(input$xs_testB, input$xs_refB),
           sprintf("%s's two levels are the same; a contrast needs two different levels.", csr_labB())))
    np <- if (isTRUE(input$xs_perm)) as.integer(input$xs_nperm %||% 200L) else 0L
    withProgress(message = "Fitting both studies", value = .25,
      detail = if (np > 0) sprintf("then %d permutations of %s", np, csr_labB()) else NULL, {
      xs_de_compare(dA, dB, xs_specA(), xs_specB(), csr_labA(), csr_labB(),
                    fdr = input$xs_fdr %||% 0.05, one_to_one_only = isTRUE(input$csr_1to1),
                    n_perm = np, seed = trace_seed())
    })
  }, module = "cross-dataset differential concordance",
  params = function(o) c(
    list(study_A = isolate(csr_labA()), study_B = isolate(csr_labB()),
         study_id_A = isolate(csr_idA()), study_id_B = isolate(csr_idB())),
    .xs_side_params(o),
    list(feature_rule = isolate(xs_bridge()$rule),
         shared_features = if (isTRUE(o$ok)) o$n_shared else NA,
         r = if (isTRUE(o$ok)) round(o$r, 3) else NA,
         r_discovery_set = if (isTRUE(o$ok) && is.finite(o$r_sig)) round(o$r_sig, 3) else NA,
         replication = if (isTRUE(o$ok) && is.finite(o$rep_rate)) sprintf("%.0f%%", 100 * o$rep_rate) else "-",
         permutations = if (isTRUE(o$ok) && !is.null(o$perm) && isTRUE(o$perm$ok)) o$perm$n else 0L,
         permutation_scheme = if (isTRUE(o$ok) && isTRUE(o$perm$ok)) o$perm$scheme else "-",
         permutation_seed = trace_seed(),
         fdr = isolate(input$xs_fdr))))

  xs_de_ok <- reactive({ o <- xs_de_res(); req(o); o })
  xs_agree_lv <- reactive({ o <- xs_de_res()
    if (is.null(o) || !isTRUE(o$ok)) character(0) else xs_agree_levels(o$labA, o$labB) })

  output$xs_de_head <- renderUI({
    o <- xs_de_res(); req(o)
    if (!isTRUE(o$ok))
      return(div(class = "alert alert-warning py-2 px-3 small", tags$b("Not run. "), o$msg))
    div(class = "alert alert-secondary py-2 px-3 small", xs_de_text(o))
  })
  output$xs_de_scatter <- renderPlotly({ o <- xs_de_res(); req(o)
    if (!isTRUE(o$ok)) return(emptyplot(o$msg))
    xs_de_scatter(o, xs_de_sc_cols()) })
  output$xs_de_quad <- renderPlot({ o <- xs_de_res(); req(o)
    xs_de_quadrant_gg(o, xs_de_quad_cols()) })
  output$xs_de_null <- renderPlotly({ o <- xs_de_res(); req(o); xs_de_null_plot(o) })
  output$xs_de_tab <- renderDT({ o <- xs_de_res(); req(o, isTRUE(o$ok))
    d <- o$tab[, c("Feature", "logFC_A", "logFC_B", "FDR_shared_A", "FDR_shared_B", "Agreement")]
    names(d)[2:5] <- c(paste(o$labA, "effect"), paste(o$labB, "effect"),
                       paste(o$labA, "FDR"), paste(o$labB, "FDR"))
    datatable(d, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE)) %>%
      formatRound(2:3, 3) %>% formatSignif(4:5, 3) })
  output$xs_de_dl <- downloadHandler(
    filename = function() "cross_dataset_differential_concordance.csv",
    content = function(f) { o <- xs_de_res(); req(o, isTRUE(o$ok)); write.csv(o$tab, f, row.names = FALSE) })
  output$xs_de_dl_txt <- downloadHandler(
    filename = function() "cross_dataset_comparison_summary.txt",
    content = function(f) { o <- xs_de_res(); req(o)
      writeLines(c(sprintf("TRACE-Omics cross-dataset comparison  (%s)", format(Sys.time())),
                   "", strwrap(xs_de_text(o), 96), "",
                   "Feature correspondence:", strwrap(xs_bridge()$note, 96), "",
                   "Timepoint correspondence:", strwrap(xs_align()$note %||% "-", 96)), f) })

  xs_de_sc_export <- reactive({ o <- xs_de_ok()
    list(summary = o$tab, csv = o$tab, ylab = "log2 FC", splitby = NULL,
         result = list(shared_features = o$tab,
                       differential_A = o$sideA$tt, differential_B = o$sideB$tt),
         ttl = sprintf("%s vs %s", o$labA, o$labB), prefix = "cross_dataset_concordance",
         gg = function(sm, yl, tt, co) xs_de_scatter_gg(isolate(xs_de_ok()), co,
                label_n = as.integer(isolate(input$xs_label_n) %||% 0L))) })
  xs_de_sc_cols <- apprServer("xs_de_sc_appr", xs_agree_lv, xs_de_sc_export)

  xs_de_quad_export <- reactive({ o <- xs_de_ok()
    list(summary = o$tab, csv = as.data.frame(o$quad), ylab = "features", splitby = NULL,
         result = list(shared_features = o$tab, contingency = as.data.frame(o$quad)),
         ttl = "agreement", prefix = "cross_dataset_agreement",
         gg = function(sm, yl, tt, co) xs_de_quadrant_gg(isolate(xs_de_ok()), co)) })
  xs_de_quad_cols <- apprServer("xs_de_quad_appr", xs_agree_lv, xs_de_quad_export)

  xs_de_null_export <- reactive({ o <- xs_de_ok()
    nulls <- if (!is.null(o$perm) && isTRUE(o$perm$ok)) data.frame(null_r = o$perm$null) else data.frame()
    list(summary = nulls, csv = nulls, ylab = "permutations", splitby = NULL,
         result = if (nrow(nulls)) list(null_draws = nulls) else NULL,
         ttl = "permutation null", prefix = "cross_dataset_null",
         gg = function(sm, yl, tt, co) xs_de_null_gg(isolate(xs_de_ok()))) })
  xs_de_null_cols <- apprServer("xs_de_null_appr", reactive(character(0)), xs_de_null_export)

  observeEvent(list(csr_dA(), csr_dB()), {
    a <- csr_dA(); b <- csr_dB()
    ac <- intersect(annotation_collections(a), annotation_collections(b))
    ac <- ac[vapply(ac, function(k) annotation_compatible(a, b, k), logical(1))]
    updateSelectInput(session, "xs_pw_coll", choices = ac, selected = ac[1])
  })
  observeEvent(list(byod_study(), input$byod_conc_l1), {
    s <- byod_study(); req(s)
    d <- s$layers[[input$byod_conc_l1 %||% s$layer_names[1]]]; req(d)
    ac <- annotation_collections(d)
    updateSelectInput(session, "byod_conc_coll", choices = ac, selected = ac[1])
    ac <- if (length(s$layers) < 2L) character() else
      intersect(annotation_collections(s$layers[[1]]), annotation_collections(s$layers[[2]]))
    if (length(ac)) ac <- ac[vapply(ac, function(k) annotation_compatible(s$layers[[1]], s$layers[[2]], k), logical(1))]
    updateSelectInput(session, "byod_pconc_coll", choices = ac, selected = ac[1])
  })
  output$xs_pw_contrast_echo <- renderUI({
    div(class = "small text-muted mb-2",
        tags$b(csr_labA()), ": ", xs_spec_text(xs_specA()), br(),
        tags$b(csr_labB()), ": ", xs_spec_text(xs_specB()))
  })
  xs_pw_res <- resettable(reactive(input$xs_pw_run), xs_gen, function() {
    dA <- csr_dA(); dB <- csr_dB()
    req(dA, dB, input$xs_testA, input$xs_refA, input$xs_testB, input$xs_refB, input$xs_pw_coll)
    withProgress(message = "Scoring both studies", value = .3, {
      xs_pathway_compare(dA, dB, xs_specA(), xs_specB(), input$xs_pw_coll,
                         input$xs_pw_mode %||% "gsea", csr_labA(), csr_labB(),
                         fdr = input$xs_fdr %||% 0.05)
    })
  }, module = "cross-dataset pathway concordance",
  params = function(o) c(
    list(.engine = .mode_engine(isolate(input$xs_pw_mode)),
         annotation = paste(isolate(csr_labA()), annotation_summary(isolate(csr_dA())),
                            isolate(csr_labB()), annotation_summary(isolate(csr_dB()))),
         study_A = isolate(csr_labA()), study_B = isolate(csr_labB()),
         study_id_A = isolate(csr_idA()), study_id_B = isolate(csr_idB())),
    .xs_side_params(o),
    list(collection = isolate(input$xs_pw_coll), scored_by = isolate(input$xs_pw_mode),
         sets = if (isTRUE(o$ok)) o$n else NA,
         r = if (isTRUE(o$ok)) round(o$r, 3) else NA)))

  xs_pw_ok <- reactive({ o <- xs_pw_res(); req(o); o })
  xs_pw_lv <- reactive({ o <- xs_pw_res()
    if (is.null(o) || !isTRUE(o$ok)) character(0) else xs_agree_levels(o$labA, o$labB) })

  output$xs_pw_head <- renderUI({
    o <- xs_pw_res(); req(o)
    if (!isTRUE(o$ok))
      return(div(class = "alert alert-warning py-2 px-3 small", tags$b("Not run. "), o$msg))
    div(class = "alert alert-secondary py-2 px-3 small", xs_pw_text(o))
  })
  output$xs_pw_scatter <- renderPlotly({ o <- xs_pw_res(); req(o)
    if (!isTRUE(o$ok)) return(emptyplot(o$msg))
    xs_pw_scatter(o, xs_pw_sc_cols()) })
  output$xs_pw_bar <- renderPlot({ o <- xs_pw_res(); req(o)
    xs_pw_bar_gg(o, as.integer(input$xs_pw_top %||% 20L), xs_pw_bar_cols()) })
  output$xs_pw_tab <- renderDT({ o <- xs_pw_res(); req(o, isTRUE(o$ok))
    d <- o$tab[, c("pathway", o$sa, o$sb, "padj_A", "padj_B", "Agreement")]
    names(d) <- c("Gene set", paste(o$labA, o$score_name), paste(o$labB, o$score_name),
                  paste(o$labA, "FDR"), paste(o$labB, "FDR"), "Agreement")
    datatable(d, rownames = FALSE, options = list(pageLength = 12, scrollX = TRUE)) %>%
      formatRound(2:3, 2) %>% formatSignif(4:5, 3) })
  output$xs_pw_dl <- downloadHandler(
    filename = function() "cross_dataset_pathway_concordance.csv",
    content = function(f) { o <- xs_pw_res(); req(o, isTRUE(o$ok)); write.csv(o$tab, f, row.names = FALSE) })

  xs_pw_sc_export <- reactive({ o <- xs_pw_ok()
    list(summary = o$tab, csv = o$tab, ylab = o$score_name, splitby = NULL,
         result = list(gene_sets = o$tab,
                       differential_A = o$sideA$tt, differential_B = o$sideB$tt),
         ttl = o$collection, prefix = "cross_dataset_pathway_concordance",
         gg = function(sm, yl, tt, co) xs_pw_scatter_gg(isolate(xs_pw_ok()), co)) })
  xs_pw_sc_cols <- apprServer("xs_pw_sc_appr", xs_pw_lv, xs_pw_sc_export)

  xs_pw_bar_export <- reactive({ o <- xs_pw_ok()
    list(summary = o$tab, csv = o$tab, ylab = o$score_name, splitby = NULL,
         result = list(gene_sets = o$tab),
         ttl = o$collection, prefix = "cross_dataset_pathway_top",
         gg = function(sm, yl, tt, co) xs_pw_bar_gg(isolate(xs_pw_ok()),
                as.integer(isolate(input$xs_pw_top) %||% 20L), co)) })
  xs_pw_bar_cols <- apprServer("xs_pw_bar_appr",
    reactive({ o <- xs_pw_res(); if (is.null(o) || !isTRUE(o$ok)) character(0)
               else xs_two_labels(o$labA, o$labB) }),
    xs_pw_bar_export)

  observeEvent(input$go_compare_upload2, go("Cross-study comparison"))
}
