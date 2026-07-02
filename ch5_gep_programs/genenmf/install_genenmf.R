#!/usr/bin/env Rscript
# install_genenmf.R
#
# Ch.5 (GEP programs) — install GeneNMF (carmonalab/GeneNMF) into the user
# library, a one-off prerequisite for run_genenmf_per_cancer.R.
#
# We use `remotes::install_github` rather than devtools: devtools drags in `fs`,
# which fails to compile in the Rstudio apptainer image (no libuv headers:
# "uv.h: No such file"). `remotes` is pure-R and has no such dependency.
# USE_BUNDLED_LIBUV=1 is set as a fallback in case any dependency still needs fs.
#
# RcppML (GeneNMF's NMF backend) is already installed in the user library.
#
# Run once inside the apptainer image before submitting the pipeline:
#   sbatch --job-name genenmf_install --time 0-01:00:00 --mem 16G --cpus-per-task 4 \
#          --output genenmf_install-%J.out --error genenmf_install-%J.err \
#          ~/rstudio/rstudio.sbatch --script genenmf/install_genenmf.R

Sys.setenv(USE_BUNDLED_LIBUV = "1")

userlib <- Sys.getenv("R_USER_LIBRARY", unset = "")
if (nzchar(userlib)) {
  .libPaths(c(userlib, .libPaths()))
  dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
}
options(repos = c(CRAN = "https://cloud.r-project.org"))
install_lib <- if (nzchar(userlib)) userlib else .libPaths()[1]

if (!requireNamespace("remotes", quietly = TRUE)) {
  # remotes is pure R (only base-package deps) — installs without compilation.
  install.packages("remotes", lib = install_lib)
}

# upgrade = "never": don't rebuild already-present deps (avoids the fs/libuv trap).
remotes::install_github("carmonalab/GeneNMF", lib = install_lib, upgrade = "never")

suppressPackageStartupMessages(library(GeneNMF))
cat("GeneNMF installed OK. Version: ",
    as.character(utils::packageVersion("GeneNMF")), "\n", sep = "")
