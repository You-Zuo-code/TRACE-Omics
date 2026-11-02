library(testthat)

test_that("contrast levels follow each selected study even when the variable is unchanged", {
  root <- normalizePath(if (file.exists("global.R")) "." else "../..")
  wd <- getwd(); setwd(root); on.exit(setwd(wd), add = TRUE)
  if (!exists("STUDIES")) suppressWarnings(suppressMessages(source("global.R")))
  server <- source("server.R", local = TRUE)$value

  make_study <- function(name, times) {
    s <- STUDIES$nmp_hum
    s$name <- name
    s$dataset_id <- paste0("controls-", name)
    s$layer_names <- "mRNA"
    s$layers <- s$layers["mRNA"]
    d <- s$layers$mRNA
    selected <- as.character(d$meta$Timepoint) %in% head(d$tp_levels, length(times))
    d$expr <- d$expr[seq_len(80), selected, drop = FALSE]
    d$meta <- d$meta[selected, , drop = FALSE]
    d$meta$Timepoint <- factor(times[match(as.character(d$meta$Timepoint),
                                         head(d$tp_levels, length(times)))], levels = times)
    d$tp_levels <- times
    d$clinical_cat <- d$clinical_num <- character(0)
    s$layers$mRNA <- d
    s$meta <- d$meta
    s$tp_levels <- times
    s
  }
  a <- make_study("A", c("D0", "D3", "D7"))
  b <- make_study("B", c("D0", "D3", "D7"))

  shiny::testServer(server, {
    sent <- list()
    session$sendInputMessage <- function(inputId, message) {
      sent[[length(sent) + 1L]] <<- list(id = inputId, message = message)
    }
    last_message <- function(id) {
      messages <- Filter(function(x) identical(x$id, id), sent)
      if (length(messages)) tail(messages, 1)[[1]]$message else NULL
    }
    session$setInputs(csr_A = "nmp_pig", csr_B = "nmp_hum",
      csr_layA = "mRNA", csr_layB = "mRNA",
      xs_varA = "Timepoint", xs_varB = "Timepoint",
      xs_testA = "6h", xs_refA = "in_vivo", xs_testB = "6h", xs_refB = "preNMP")
    byod_ws(list(ds1 = list(study = a), ds2 = list(study = b)))
    session$flushReact()
    sent <- list()
    session$setInputs(csr_A = "ds1", csr_B = "ds2")
    for (side in c("A", "B")) {
      test <- last_message(paste0("xs_test", side))
      ref <- last_message(paste0("xs_ref", side))
      expect_identical(test$value, "D7")
      expect_identical(ref$value, "D0")
      expect_true(is.character(test$options) && grepl('value="D3"', test$options, fixed = TRUE))
      expect_false(any(grepl("in_vivo|preNMP|6h", c(test$options, ref$options))))
    }

    sent <- list()
    byod_ws(list(ds1 = list(study = a), ds2 = list(study = make_study("B-new", c("V0", "V1", "V2")))))
    session$flushReact()
    expect_identical(last_message("xs_testB")$value, "V2")
    expect_identical(last_message("xs_refB")$value, "V0")
    expect_identical(last_message("xs_testA")$value, "D7")

    sent <- list()
    session$setInputs(csr_A = "nmp_pig")
    expect_identical(last_message("xs_testA")$value, "6h")
    expect_identical(last_message("xs_refA")$value, "in_vivo")
    session$setInputs(xs_varA = "WI_status")
    expect_identical(last_message("xs_testA")$value, tail(group_levels(csr_dA(), "WI_status"), 1))
    expect_identical(last_message("xs_refA")$value, head(group_levels(csr_dA(), "WI_status"), 1))
  })
})
