trace_home <-
function (reference_button, upload_button) 
{
    reference_button$children <- list("Open kidney reference")
    upload_button$children <- list("Upload a study")
    t <- htmltools::tags
    t$div(class = "trace-home", t$div(class = "trace-home-top", t$section(class = "trace-home-intro", t$p(class = "trace-kicker", 
        "TRACE-Omics / Analysis workspace"), t$h1("Longitudinal multi-omics analysis"), t$p("Analyse molecular changes across repeated samples, compare treatment trajectories, ", 
        "and examine their association with clinical variables. Use the kidney cohorts below ", "or bring your own processed data."), 
        t$div(class = "trace-entry-actions", reference_button, upload_button)), t$aside(class = "trace-home-aside", 
        `aria-label` = "Data and outputs", t$p(class = "trace-kicker", "Working with your data"), t$dl(t$dt("Input"), 
            t$dd("Processed molecular matrices and sample metadata (CSV or TSV)."), t$dt("Study design"), 
            t$dd("Timepoints, subjects, paired organs and groups."), t$dt("Download"), t$dd("Figures, full result tables, plotted data and analysis records.")))), 
        t$section(`aria-labelledby` = "trace-reference-title", t$div(class = "trace-section-heading", t$h2(id = "trace-reference-title", 
            "Kidney perfusion reference"), t$span("Two built-in longitudinal designs")), t$div(class = "trace-cohorts", 
            t$article(class = "trace-cohort", t$p(class = "trace-kicker", "Porcine / paired experiment"), 
                t$h3("Two kidneys, different treatments"), t$p("Warm ischaemia and control are compared within the same animal."), 
                t$dl(t$dt("Study size"), t$dd("30 animals · 60 kidneys"), t$dt("Sampling"), t$dd("Six timepoints, including an in vivo baseline"), 
                  t$dt("Molecular layers"), t$dd("Transcriptome and proteome"))), t$article(class = "trace-cohort", 
                t$p(class = "trace-kicker", "Human / clinical cohort"), t$h3("Trajectories in their clinical context"), 
                t$p("Study continuous clinical associations or compare donor-defined groups."), t$dl(t$dt("Study size"), 
                  t$dd("24 donor kidneys"), t$dt("Sampling"), t$dd("Five timepoints · 33 clinical variables"), 
                  t$dt("Molecular layers"), t$dd("Transcriptome and proteome"))))), t$section(`aria-labelledby` = "trace-tools-title", 
            t$div(class = "trace-section-heading", t$h2(id = "trace-tools-title", "Analysis modules"), t$span("Choose a module after loading a study")), 
            t$div(class = "trace-tools", t$div(class = "trace-tool", t$h3("Time courses and paired comparisons"), 
                t$p("Differential analysis, group × time, temporal programmes and baseline sensitivity.")), 
                t$div(class = "trace-tool", t$h3("Clinical variables"), t$p("Continuous associations, stratified trajectories and variable-by-time screening.")), 
                t$div(class = "trace-tool", t$h3("Pathways and molecular layers"), t$p("Enrichment, pathway activity, cross-layer concordance and DIABLO.")), 
                t$div(class = "trace-tool", t$h3("Comparison between studies"), t$p("Fit a contrast in each study, then compare feature effects or pathway responses.")))), 
        t$footer(class = "trace-home-footer", t$span("R / Shiny · Processed-data analysis"), t$span("Upload → Declare the design → Analyse → Download")))
}

trace_navigation <- function(node) {
  if (inherits(node, 'shiny.tag')) {
    labels <- c('Kidney Perfusion Reference' = 'Kidney reference', 'Cross-study comparison' = 'Compare studies')
    value <- node$attribs[['data-value']]
    if (identical(node$name, 'a') && length(value) == 1L && value %in% names(labels))
      node$children <- list(unname(labels[value]))
    node$children <- lapply(node$children, trace_navigation)
  } else if (is.list(node) && !inherits(node, 'html_dependency')) {
    for (i in seq_along(node)) node[i] <- list(trace_navigation(node[[i]]))
  }
  node
}
