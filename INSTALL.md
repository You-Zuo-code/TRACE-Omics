# Installation

TRACE-Omics runs locally as an R/Shiny application. Install R, install the R packages, then start the application. RStudio is optional. No Python environment, Docker image, Node.js installation or database server is needed to run the application.

## Requirements and versions

- R 4.5 or later, a web browser, internet access during installation, and a writable local project folder.
- The default installer accepts newer R versions and selects packages through CRAN and Bioconductor. It does not force R 4.5.1 or use the lockfile.
- R 4.5.1 / Bioconductor 3.21 is the recorded reference environment. The optional locked installer is described below.
- Storage and memory requirements depend on the input matrices and analysis. The source download is much smaller than the installed package library; compilation also needs temporary space. No minimum-RAM guarantee has been established.

## Get the source

Download **Code > Download ZIP** from the repository, then extract it completely. Alternatively, download a tagged source ZIP from **Releases**. Work inside the extracted directory containing `install.R`, `run.R`, `DESCRIPTION` and `data/atlas.rds`, not inside the ZIP or its parent directory.

Keep the `R/`, `tools/`, `tests/`, `data/` and `www/` directories. Do not download only the three Shiny files. Keep the bundled data and resource licence notices with the source.

## Windows 10 / Windows 11

1. Install 64-bit R from [CRAN for Windows](https://cran.r-project.org/bin/windows/base/). Use an R release supported by your Windows version. On Windows ARM, distinguish an x64 R installation running under emulation from native ARM R.
2. Optionally install [RStudio Desktop](https://posit.co/download/rstudio-desktop/). RStudio is an editor, not a replacement for R. Check the active version in its R console with `R.version.string`.
3. Windows binary packages do not require a compiler. If a dependency must be built from source, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) as specified for your R release. Rtools45 currently covers R 4.5.0 and later, including R 4.6. Check the official table for other R versions. Do not install an arbitrary GCC toolchain.
4. Extract the project into a writable folder, for example `C:/work/TRACE-Omics`. Avoid installing it inside `Program Files` or a read-only/shared folder.

In the R or RStudio console, select the extracted project directory:

```r
app_dir <- choose.dir(caption = "Select the TRACE-Omics folder containing install.R")
if (is.na(app_dir)) stop("No folder selected")
setwd(app_dir)
stopifnot(file.exists("install.R"), file.exists("data/atlas.rds"))
source("install.R")
```

Wait for `Installation complete`. Then run:

```r
source("run.R")
```

For PowerShell, if `Rscript` is not on PATH, use its installed full path. This example selects R 4.5.1; replace it with your installed version rather than installing that version just for this command:

```powershell
Set-Location "C:\work\TRACE-Omics"
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" install.R
```

After installation succeeds:

```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" run.R
```

If several R versions are installed, start the intended R executable or select it in RStudio's R-version settings. A newer R version is allowed by the default installer.

## macOS

1. Install R from [CRAN for macOS](https://cran.r-project.org/bin/macosx/). Choose **arm64** for Apple silicon or **x86_64** for Intel, and check that the selected release supports your macOS version.
2. Optionally install RStudio Desktop for the same architecture.
3. CRAN binaries usually avoid compilation. For packages requiring source builds, install Apple's Command Line Tools from Terminal:

```sh
xcode-select --install
```

For packages containing Fortran, use the compiler specified for your R build on the [R for macOS tools page](https://mac.r-project.org/tools/). Xcode does not include that compiler. Do not mix a Homebrew-built R/package toolchain with CRAN R binaries without managing that toolchain consistently.

In Terminal, replace the example directory with the extracted folder:

```sh
cd "/path/to/TRACE-Omics"
Rscript --version
Rscript install.R
```

After installation succeeds:

```sh
Rscript run.R
```

Alternatively, in the R/RStudio console choose `install.R` in the file dialog:

```r
setwd(dirname(file.choose()))
stopifnot(file.exists("install.R"), file.exists("data/atlas.rds"))
source("install.R")
```

Then run `source("run.R")` after installation completes.

## Ubuntu / Debian Linux

Use R 4.5 or later. Distribution repositories may supply an older R version, so check `Rscript --version` first. Follow the [CRAN Ubuntu instructions](https://cran.r-project.org/bin/linux/ubuntu/) or [CRAN Debian instructions](https://cran.r-project.org/bin/linux/debian/) for your distribution. Do not add an Ubuntu repository to Debian.

For **Ubuntu 24.04 (noble)**, the following sets up CRAN's current R repository. It installs the current repository version, not necessarily R 4.5.1. These are system-administration commands requiring sudo; do not use them on a managed machine without its administrator's agreement.

```sh
sudo apt-get update
sudo apt-get install --no-install-recommends ca-certificates curl gnupg software-properties-common
TRACE_CRAN_KEY="$(mktemp -t trace-cran-key.XXXXXX)"
curl -fsSL https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc -o "$TRACE_CRAN_KEY"
gpg --show-keys --with-fingerprint "$TRACE_CRAN_KEY"
```

Check the fingerprint against CRAN's official page before continuing. At preparation time it is `E298 A3A8 25C0 D65D FD57 CBB6 5171 6619 E084 DAB9`. Stop if the download fails or the fingerprint differs. Only after verification, in the same terminal, install the key and add the repository below. For other Ubuntu releases, use their codename and supported-version instructions from CRAN, not `noble`.

```sh
sudo install -m 644 "$TRACE_CRAN_KEY" /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu noble-cran40/"
sudo apt-get update
sudo apt-get install --no-install-recommends r-base r-base-dev
Rscript --version
```

Install build and graphics dependencies on Debian/Ubuntu:

```sh
sudo apt-get install build-essential gfortran cmake pkg-config \
  libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev \
  libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
  libjpeg-dev libtiff-dev libwebp-dev libicu-dev libmagick++-dev \
  libgsl-dev libglpk-dev libgmp-dev fonts-liberation fonts-noto-cjk
```

These are the expected system dependencies. A package release may need additional system libraries; use the first package configuration error to identify any missing library. Other Linux distributions use different package names and package managers.

CRAN does not provide binary R packages for Linux, so the default installer compiles them from source. Expect this to take much longer than on Windows or macOS. As your normal user, not through `sudo Rscript`, install the project packages:

```sh
cd "/path/to/TRACE-Omics"
Rscript install.R
```

After installation succeeds:

```sh
Rscript run.R
```

In a Linux VM, open the URL in the guest browser. For a remote host, keep the application bound to loopback and use an SSH tunnel from your own computer, for example `ssh -L 3838:127.0.0.1:3838 user@server`; then browse to the local URL. Do not expose an unauthenticated development server to the internet. A managed Shiny service is a separate deployment task.

## Open, check and stop

Open **http://127.0.0.1:3838/** on the computer running R. This address is the same on every computer: it means "this computer", not a public website, and your data stay on your machine. The terminal should show `Listening on ...`; it remains occupied while the server runs. That is normal. In RStudio an interactive launch can open the browser automatically.

Load a human or porcine cohort from **Kidney reference**, open **Overview**, run a small analysis, and download a result to confirm the complete workflow on your machine. For automated checks, stop the server first and run from the project directory:

```sh
Rscript tests/run.R
```

Or in an R console: `source("tests/run.R")`. Retain any warnings and the first error when reporting a failure. A passing installation alone is not a check of every analysis on your uploaded data.

Stop a terminal launch with **Ctrl+C**. Use **Stop** or **Esc** in the R console. Download session results before stopping. Later launches need only `Rscript run.R` or `source("run.R")`, not another installation.

## Package library and reference installation

Packages are installed into `.library/<R-version>-<platform>` in the project. The installer puts that library first and can reuse packages already available in other libraries. It does not replace system libraries. Keep `.library` locally; do not commit it to GitHub or copy it between operating systems/CPU architectures. After changing R version, rerun `install.R` for that version.

The default mode uses [BiocManager](https://bioconductor.org/install/) to choose packages for the active R. A recent package update can introduce incompatibilities; allowing a version is not a promise that all future combinations will work.

For the recorded dependency versions, use **a separate copy of the source with R 4.5.x**:

```sh
Rscript install.R --locked
```

In an R console:

```r
Sys.setenv(TRACE_INSTALL_MODE = "locked")
tryCatch(source("install.R"), finally = Sys.unsetenv("TRACE_INSTALL_MODE"))
```

Only locked mode requires R 4.5.x. It restores and verifies `renv.lock`; it is not a way to force R 4.6 or a different Bioconductor release to use those old packages. The lockfile does not install operating-system libraries, and future availability of historical downloads is outside this repository's control.

## Troubleshooting

| Symptom | What to check |
|---|---|
| `Run install.R from the TRACE-Omics directory` | Run `getwd()` and `list.files()`. The working directory must contain `DESCRIPTION` and `tools/install_helpers.R`. |
| `Rscript` is not recognised | Use the R/RStudio console route or the full path to the Rscript executable. |
| R 4.5.x requested despite using a newer R | Remove `--locked`, or run `Sys.unsetenv("TRACE_INSTALL_MODE")` before `source("install.R")`. Confirm the active version with `R.version.string`. |
| `org.Hs.eg.db` download/install failure | Read the first error above the final failure summary. Check network access, free disk space and permissions. Large annotation packages may take time to download and unpack. Rerun the installer after resolving the error; do not mix an unrelated package version into locked mode. |
| Download cannot reach the default Bioconductor mirror | In R, run `Sys.setenv(TRACE_BIOC_MIRROR = "https://bioconductor.org")`, then `source("install.R")`. Do not disable TLS verification. |
| Compilation fails | Check Rtools on Windows, the matching macOS compiler tools, or the missing Linux `-dev` package named in the error. |
| An installation directory is locked | First stop other R installations using that project library. Only if no install is running, inspect and remove the specific stale `00LOCK...` directory inside this project's `.library`, then retry. Do not delete the entire library. |
| Port 3838 is in use | In R, run `Sys.setenv(TRACE_PORT = "3839")` before `source("run.R")`; open `http://127.0.0.1:3839/`. |
| Empty page or connection refused | Check that R is still running and shows a listening address. The browser and server must refer to the same machine/forwarded port. Read the R error before restarting. |
| Missing or substituted export font | Install a suitable font. Exports prefer Arial/Helvetica and use installed fallbacks such as Liberation Sans. `Sys.setenv(TRACE_PLOT_FONT = "Liberation Sans")` before launch requests an installed font; it does not install one. |

For an issue report, include your OS and CPU architecture, `R.version.string`, `R.version$platform`, `sessionInfo()`, current or locked installation mode, the release/commit used, and the first relevant error. The installer also saves `installed-packages.csv` inside its project library. Redact personal paths, credentials and unpublished patient/sample data before sharing logs or screenshots.
