#!/usr/bin/env Rscript
log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

RDS     <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Brain/Epithelial_Brain_merged_full_scanvi_final.RDS"
OUT_DIR <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Brain/brain_mtx"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

log_msg("Loading Brain Malignant RDS...")
obj <- readRDS(RDS)
log_msg(sprintf("Cells: %d  Genes: %d", ncol(obj), nrow(obj)))

# Extract counts
lyrs <- Layers(obj[["RNA"]])
count_lyrs <- grep("^counts", lyrs, value = TRUE)
if (length(count_lyrs) == 0) {
  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
} else {
  counts <- LayerData(obj, assay = "RNA", layer = count_lyrs[1])
}

sample_vals <- counts@x[seq_len(min(200, length(counts@x)))]
log_msg(sprintf("Count range sample: %.3f - %.3f (integers: %s)",
                min(sample_vals), max(sample_vals),
                all(sample_vals == floor(sample_vals))))

log_msg("Writing mtx...")
writeMM(counts, file.path(OUT_DIR, "matrix.mtx"))
write.table(data.frame(gene = rownames(counts)),
            file.path(OUT_DIR, "genes.tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(data.frame(barcode = colnames(counts)),
            file.path(OUT_DIR, "barcodes.tsv"), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.csv(obj@meta.data, file.path(OUT_DIR, "metadata.csv"), row.names=TRUE)

log_msg(sprintf("Done. Output: %s", OUT_DIR))
