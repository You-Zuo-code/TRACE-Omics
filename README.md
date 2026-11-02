# TRACE-Omics

**Temporal and Repeated-measures Analysis of Cross-layer Experiments**

TRACE-Omics is an R/Shiny application for longitudinal multi-omics analysis. It supports repeated samples, paired organs, clinical variables and matched molecular layers. Uploaded studies can be analysed individually or compared using separately fitted contrasts and pathway responses.

The built-in reference contains porcine and human kidney perfusion data, each with mRNA and protein measurements. Uploads use processed expression matrices, not raw reads or mass spectra.

## Install and run

TRACE-Omics runs on your own computer. You need [R](https://cran.r-project.org/) 4.5 or later and a web browser. RStudio is optional. Python, Docker and a database server are not required.

Installation needs internet access and about 1.4 GB of disk space. You only need to install once.

### Windows

1. Install R from [CRAN for Windows](https://cran.r-project.org/bin/windows/base/).
2. On this page, click **Code > Download ZIP**. Extract the ZIP completely, for example to `C:\work`.
3. Open R or RStudio and run the lines below. When asked, select the folder that directly contains `install.R`.

   ```r
   setwd(choose.dir(caption = "Select the folder containing install.R"))
   source("install.R")
   ```

4. When `Installation complete` appears, start the application:

   ```r
   source("run.R")
   ```

### macOS

1. Install R from [CRAN for macOS](https://cran.r-project.org/bin/macosx/). Choose **arm64** for Apple silicon (M1 or later) or **x86_64** for Intel, and check that the release supports your macOS version.
2. On this page, click **Code > Download ZIP** and extract it.
3. Open Terminal and install. Replace the path with your extracted folder.

   ```sh
   cd ~/Downloads/TRACE-Omics-main
   Rscript install.R
   ```

4. When `Installation complete` appears, start the application from the same folder:

   ```sh
   Rscript run.R
   ```

If you prefer RStudio, choose **Session > Set Working Directory > Choose Directory**, select the folder containing `install.R`, then run `source("install.R")`, followed by `source("run.R")`.

### Linux

Linux needs system libraries and builds packages from source. See [INSTALL.md](INSTALL.md).

### Open the application

Started from R or RStudio, the application normally opens in your browser. Started from a terminal, open <http://127.0.0.1:3838/> yourself.

This address is the same on every computer. It means "this computer", so it is not a public website: the application runs in your own R session and your data stay on your machine. The page works only while R is running, so leave R open. To stop, press **Ctrl+C** in a terminal, or **Esc** or **Stop** in R or RStudio.

Next time, skip the installation and start with `Rscript run.R` or `source("run.R")` from the same folder. For troubleshooting, see [INSTALL.md](INSTALL.md).

### Installation details

The default installer uses CRAN and the Bioconductor release selected for your R version. It accepts newer R releases without an upper-version restriction. Packages are installed into the project-local `.library`; existing system libraries are not replaced. If interrupted, rerun it.

To check your installation, stop the application and run `Rscript tests/run.R`, or `source("tests/run.R")` in R.

### Reproduce the reference environment

If a newer environment has dependency problems, use R 4.5.x (tested: 4.5.1) and the optional lockfile installer in a separate copy of this repository:

```sh
Rscript install.R --locked
```

In RStudio:

```r
Sys.setenv(TRACE_INSTALL_MODE = "locked")
tryCatch(source("install.R"), finally = Sys.unsetenv("TRACE_INSTALL_MODE"))
```

Only this optional mode requires R 4.5.x. It restores and checks the versions in `renv.lock`. The default installation does not require a downgrade. After a successful installation, start with `Rscript run.R` or `source("run.R")` as above.

### System dependencies

Windows usually uses binary packages; source builds need [Rtools](https://cran.r-project.org/bin/windows/Rtools/) matching the R release. macOS source builds may need Xcode Command Line Tools and a matching Fortran compiler. See the [R installation manual](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).

On Debian/Ubuntu:

```sh
sudo apt-get update
sudo apt-get install build-essential gfortran cmake pkg-config \
  libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
  libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
  libjpeg-dev libtiff-dev libwebp-dev libicu-dev libmagick++-dev \
  libgsl-dev libglpk-dev libgmp-dev fonts-liberation fonts-noto-cjk
```

## Use

1. Load a cohort from **Kidney reference**, or open **Your data**.
2. Upload one to three expression matrices and a metadata table. The first matrix column contains feature IDs; remaining columns contain samples. The first metadata column contains the same sample IDs.
3. Declare the time variable and order, subject, optional nesting unit, measurement scale and normalization.
4. Open **Analyse** for quality checks, contrasts, trajectories, clinical associations or multi-omics integration.
5. For cross-study analysis, load studies into separate slots and select each contrast in **Compare studies**.

CSV/TSV uploads may be gzip-compressed. The upload examples contain the same human reference cohort as either mRNA alone or matched mRNA and protein; they are not independent studies.

The random seed defaults to 1234. Change it in **Your data**, click **Apply**, then rerun affected analyses. Each captured result retains its own recorded seed until replaced.

Feature and sample filters are off by default. Optional upload filters act on the supplied values before normalization. In **Differential**, **Exclude constant features** removes features with exactly zero variance across the profiles selected for that contrast, before model fitting and BH adjustment. It is off by default and does not change the input matrix or other analyses. The result archive records the setting, removed IDs and tested-feature count. It is not a low-expression or near-zero-variance filter; zero residual variances can still arise in fitted models.

Declare the source organism and feature ID type for each layer. Human annotation and offline pig, mouse and rat ortholog mappings are included. Other species or identifier systems can use a custom mapping or GMT. See [annotation formats](data/annotation.md) for coverage and matching rules.

## Export

Download figures, plotted data or complete result tables from each analysis. **Download results for the active study (.zip)** on the **Reproducibility** page packages the latest captured results, inputs, annotations and run settings. A session can contain captures made with different settings; rerun and recapture panels when changing the seed or analysis parameters. Download results before closing the session or replacing a study.

Cross-study captures contain derived results from both studies. Input matrices and metadata from comparison studies are excluded by default; the archive checkbox explicitly enables their inclusion. Review the archive before sharing it.

Paired contrasts use only declared subject or higher-level pairing IDs. For independent clinical groups, turn **Paired** off. Incomplete pairs and their excluded profile counts are reported with the contrast; the application does not choose a metadata column as a substitute pairing variable.

Exports prefer installed Arial or Helvetica, with Liberation Sans or another installed sans-serif font as a fallback. Set `TRACE_PLOT_FONT` before launch to request a specific installed font. PDF text remains editable text; font substitution in a graphics editor still depends on that editor and its installed fonts.

The server listens locally by default. Set `TRACE_PORT` to change the port.

## Licence and citation

Code is licensed under [MIT](LICENSE). Data and annotation resources have separate [sources and terms](data/README.md), including [gene-set terms](data/genesets_LICENSE.txt). Cite the software and the studies used in your analysis.
