library(testthat)

test_that("two uploaded files build, compare, export and survive replacement independently", {
  root <- normalizePath(if(file.exists("global.R")) "." else "../..")
  wd <- getwd();setwd(root);on.exit(setwd(wd),add=TRUE)
  if(!exists("STUDIES")) suppressWarnings(suppressMessages(source("global.R")))
  server <- source("server.R",local=TRUE)$value
  dir <- tempfile("two-upload-");dir.create(dir)
  on.exit(unlink(dir,recursive=TRUE),add=TRUE)
  make_files <- function(prefix,seed) {
    set.seed(seed)
    meta <- expand.grid(Timepoint=c("t0","t1","t2","t3"),Subject=paste0("s",1:8),stringsAsFactors=FALSE)
    meta$sample <- paste(prefix,meta$Subject,meta$Timepoint,sep="_")
    x <- matrix(rnorm(100*nrow(meta),8,.3),100,nrow(meta))
    x[1:20,] <- x[1:20,] + rep(as.integer(factor(meta$Timepoint))-1,each=20)
    x <- data.frame(Feature=paste0("GENE",1:100),x,check.names=FALSE)
    names(x)[-1] <- meta$sample
    ex <- file.path(dir,paste0(prefix,"_expr.csv")); me <- file.path(dir,paste0(prefix,"_meta.csv"))
    write.csv(x,ex,row.names=FALSE);write.csv(meta[,c("sample","Timepoint","Subject")],me,row.names=FALSE)
    file <- function(p) data.frame(name=basename(p),size=file.info(p)$size,type="text/csv",datapath=p)
    list(expr=file(ex),meta=file(me))
  }
  fa<-make_files("A",1);fb<-make_files("B",2)
  shiny::testServer(server, {
    session$setInputs(byod_target="ds1",byod_expr1=fa$expr,byod_meta=fa$meta,
      byod_type1="transcriptomics",byod_name1="RNA",byod_unit1="log2 intensity",byod_norm1="None",
      byod_idmap1="none",byod_tpcol="Timepoint",byod_subj="Subject",byod_nest="(none)",
      byod_batch="(none)",byod_organism="(not declared)",byod_studyname="\u7814\u7a76 (paired)")
    session$setInputs(byod_build=1);session$flushReact()
    a<-byod_ws()$ds1$study
    expect_false(is.null(a));expect_equal(ncol(a$layers$RNA$expr),32L)
    session$setInputs(byod_target="ds2",byod_expr1=fb$expr,byod_meta=fb$meta)
    session$setInputs(byod_tpcol="Timepoint",byod_subj="Subject",byod_studyname="\u7814\u7a76 (paired)",byod_build=2)
    session$flushReact();b<-byod_ws()$ds2$study
    expect_false(is.null(b));expect_identical(byod_ws()$ds1$study$dataset_id,a$dataset_id)
    expect_false(identical(a$dataset_id,b$dataset_id))
    session$setInputs(csr_A="ds1",csr_B="ds2",csr_layA="RNA",csr_layB="RNA",
      xs_varA="Timepoint",xs_varB="Timepoint",xs_testA="t3",xs_refA="t0",xs_testB="t3",xs_refB="t0",
      xs_paired=TRUE,xs_perm=FALSE,xs_fdr=.05,csr_1to1=FALSE)
    session$setInputs(xs_de_run=1);session$flushReact()
    result<-xs_de_res();expect_true(result$ok);expect_equal(result$n_shared,100L)
    expect_equal(result$sideA$n_units,8);expect_equal(result$sideB$n_units,8)
    pair<-bundle_pair_folder(a$dataset_id,b$dataset_id)
    expect_true(bundle_capture(session,"xs_de_sc_appr",xs_de_sc_export(),xs_de_scatter_gg(result),BUNDLE_PAIR_DATASET,pair))
    zipfile<-file.path(dir,"A_results.zip")
    expect_true(bundle_assemble(session,a$dataset_id,zipfile,log_records=analysis_log(),study_name=a$name,
      study=a,related_studies=list(b)))
    entries<-unzip(zipfile,list=TRUE)$Name
    expect_true(any(grepl("\\.csv$",entries)));expect_true(any(grepl("\\.pdf$",entries)))
    actual <- entries[grepl("result_shared_features.csv$", entries)]
    expect_length(actual, 1L); expect_equal(nrow(read.csv(unz(zipfile, actual))), 100L)
    related <- entries[grepl("comparison_annotation/", entries)]
    expect_gt(length(related), 0L)
    expect_false(any(grepl("input_expression|input_metadata", related)))
    expect_true(any(grepl("^annotation/.*/input_expression.csv.gz$", entries)))
    consented <- file.path(dir,"A_with_comparison_inputs.zip")
    expect_true(bundle_assemble(session,a$dataset_id,consented,study=a,related_studies=list(b),
      include_comparison_inputs=TRUE))
    included <- unzip(consented,list=TRUE)$Name
    expect_true(any(grepl("^comparison_annotation/.*/input_expression.csv.gz$", included)))
    expect_true(any(grepl("^comparison_annotation/.*/input_metadata.csv$", included)))
    previous<-xs_gen()
    session$setInputs(byod_target="ds1",byod_build=3);session$flushReact()
    expect_gt(xs_gen(),previous)
    expect_identical(byod_ws()$ds2$study$dataset_id,b$dataset_id)
    expect_false(identical(byod_ws()$ds1$study$dataset_id,a$dataset_id))
  })
})
