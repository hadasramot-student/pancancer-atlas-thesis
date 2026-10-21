#!/usr/bin/env Rscript
# export_breast_chunks_to_mtx.R
#
# Exports raw counts from Breast Epithelial chunks directly to .mtx files,
# one per chunk, after subsetting to common genes.
# Python will stack them into a single matrix.

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

CHUNK_DIR <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/tmp_merge/Epithelial"
OUT_DIR   <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast/epithelial_mtx_chunks"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

chunk_files <- sort(list.files(CHUNK_DIR, pattern = "^chunk_.*\\.RDS$", full.names = TRUE))
log_msg(sprintf("Found %d chunk files", length(chunk_files)))

# ── Step 1: find common genes ──────────────────────────────────────────────
log_msg("Scanning gene names across all chunks...")
gene_lists <- lapply(chunk_files, function(f) {
  obj <- readRDS(f)
  g <- rownames(obj)
  rm(obj); gc()
  g
})
common_genes <- Reduce(intersect, gene_lists)
log_msg(sprintf("Common genes: %d", length(common_genes)))
rm(gene_lists); gc()

# write gene list once
write.table(data.frame(gene = common_genes),
            file.path(OUT_DIR, "genes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# ── Step 2: export each chunk as mtx ──────────────────────────────────────
all_barcodes <- c()

for (i in seq_along(chunk_files)) {
  f <- chunk_files[[i]]
  log_msg(sprintf("  Chunk %d/%d: %s", i, length(chunk_files), basename(f)))
  obj <- readRDS(f)
  obj <- obj[common_genes, ]

  # get counts — handle both joined and split layers
  lyrs <- Layers(obj[["RNA"]])
  count_lyrs <- grep("^counts", lyrs, value = TRUE)
  if (length(count_lyrs) == 1) {
    counts <- LayerData(obj, assay = "RNA", layer = count_lyrs[1])
  } else {
    # multiple split layers — cbind them (each split layer is small enough)
    counts <- do.call(cbind, lapply(count_lyrs, function(l) {
      LayerData(obj, assay = "RNA", layer = l)
    }))
  }

  log_msg(sprintf("  %d cells x %d genes", ncol(counts), nrow(counts)))

  chunk_out <- file.path(OUT_DIR, sprintf("chunk_%03d", i))
  dir.create(chunk_out, showWarnings = FALSE)
  writeMM(counts, file.path(chunk_out, "matrix.mtx"))
  write.table(data.frame(barcode = colnames(counts)),
              file.path(chunk_out, "barcodes.tsv"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)

  all_barcodes <- c(all_barcodes, colnames(counts))
  rm(obj, counts); gc()
  log_msg(sprintf("  Saved chunk %d", i))
}

# write combined barcode list
write.table(data.frame(barcode = all_barcodes),
            file.path(OUT_DIR, "barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

log_msg(sprintf("Done! Total barcodes: %d", length(all_barcodes)))
log_msg(sprintf("Output: %s", OUT_DIR))
