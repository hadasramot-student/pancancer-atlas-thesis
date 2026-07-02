#!/usr/bin/env Rscript
# export_hallmark_gmt.R
#
# Ch.5 (GEP programs) — Spectra pipeline step C. Export MSigDB Hallmark
# (collection H) gene sets to a JSON of {set_name: [gene_symbols]} so the
# Python prior-builder (build_spectra_priors.py) can assemble the Spectra
# nested gene_set_dictionary. Uses the same msigdbr source as
# ../../ch4_subclustering/gsea_hallmarks.R.
#
# Usage:
#   Rscript export_hallmark_gmt.R --out ch5_gep_programs/spectra/priors/hallmark.json

userlib <- Sys.getenv("R_USER_LIBRARY", unset = "")
if (nzchar(userlib)) .libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(optparse)
  library(msigdbr)
  library(jsonlite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--out", type = "character",
              default = "ch5_gep_programs/spectra/priors/hallmark.json",
              help = "Output JSON path"),
  make_option("--species", type = "character", default = "Homo sapiens")
)))

message("Fetching MSigDB Hallmark (H) for ", opt$species)
h <- msigdbr(species = opt$species, collection = "H")
pathways <- split(h$gene_symbol, h$gs_name)
pathways <- lapply(pathways, function(x) sort(unique(x)))
message("Loaded ", length(pathways), " Hallmark gene sets")

dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write_json(pathways, opt$out, auto_unbox = FALSE, pretty = TRUE)
message("Wrote ", opt$out)
