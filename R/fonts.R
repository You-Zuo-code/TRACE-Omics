trace_font <- function(requested = "", cjk = FALSE, available = NULL) {
  if (is.null(available)) available <- unique(systemfonts::system_fonts()$family)
  preferred <- if (cjk)
    c("Noto Sans CJK SC", "Noto Sans CJK TC", "Noto Sans CJK JP", "Noto Sans CJK KR",
      "Microsoft YaHei", "Yu Gothic", "Malgun Gothic", "Songti SC", "PingFang SC") else
    c("Arial", "Helvetica", "Liberation Sans", "Noto Sans", "DejaVu Sans")
  choices <- unique(c(requested[nzchar(requested)], preferred))
  if (cjk) {
    choices <- choices[choices %in% available]
    # An explicit Latin font must not disable fallback for uploaded Unicode labels.
    choices <- choices[vapply(choices, function(f)
      all(systemfonts::glyph_info(c("\u4f9b", "\u4f53"), family = f)$index != 0), logical(1))]
  }
  matched <- match(tolower(choices), tolower(available))
  matched <- matched[!is.na(matched)]
  if (length(matched)) return(available[matched[1]])
  if (cjk) stop("No CJK font is available. Install Noto Sans CJK or set TRACE_PLOT_FONT to a font covering your labels.", call. = FALSE)
  "sans"
}
