#!/usr/bin/env Rscript
# build_breast_final_seurat.R
#
# Merges 13 Breast Epithelial chunk RDS files (split layers, no JoinLayers),
# attaches scANVI latent + labels, saves final RDS.
#
# Key: Seurat v5 merge() keeps split layers — avoids 32-bit overflow entirely.

log_msg <- function(msg) cat(sprintf("[%s] %s\n", Sys.time(), msg))

userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

BASE      <- "/ofir_grp/users/hadasra/temp_objects/merged_by_query/Breast"
CHUNK_DIR <- file.path(BASE, "tmp_merge/Epithelial")
LATENT    <- file.path(BASE, "Epithelial_Breast_merged_full_scanvi_scanvi_latent.csv")
LABELS    <- file.path(BASE, "Epithelial_Breast_merged_full_scanvi_scanvi_labels.csv")
OUT_RDS   <- file.path(BASE, "Epithelial_Breast_merged_full_scanvi_final.RDS")

chunk_files <- sort(list.files(CHUNK_DIR, pattern = "chunk_.*\\.RDS$", full.names = TRUE))
log_msg(sprintf("Found %d chunk files", length(chunk_files)))

# ── Load all chunks ────────────────────────────────────────────────────────
log_msg("Loading all chunks...")
obj_list <- lapply(seq_along(chunk_files), function(i) {
  log_msg(sprintf("  Chunk %d/%d: %s", i, length(chunk_files), basename(chunk_files[i])))
  readRDS(chunk_files[i])
})

# ── Merge (split layers — no JoinLayers = no 32-bit overflow) ─────────────
log_msg("Merging chunks (split layers kept)...")
obj <- merge(obj_list[[1]], y = obj_list[-1])
rm(obj_list); gc()
log_msg(sprintf("Merged: %d cells x %d genes", ncol(obj), nrow(obj)))
log_msg(sprintf("Layers: %s", paste(Layers(obj), collapse=", ")))

# ── Attach scANVI latent embedding ────────────────────────────────────────
log_msg("Reading latent CSV...")
latent <- read.csv(LATENT, row.names = 1, check.names = FALSE)
log_msg(sprintf("Latent: %d cells x %d dims", nrow(latent), ncol(latent)))

common <- intersect(colnames(obj), rownames(latent))
log_msg(sprintf("Cells matched: %d / %d", length(common), ncol(obj)))

latent_mat <- as.matrix(latent[common, ])
colnames(latent_mat) <- paste0("SCANVI_", seq_len(ncol(latent_mat)))
rownames(latent_mat) <- common

obj[["scanvi"]] <- CreateDimReducObject(
  embeddings = latent_mat,
  key        = "SCANVI_",
  assay      = "RNA"
)
log_msg("scANVI reduction added.")

# ── Attach scANVI labels ──────────────────────────────────────────────────
log_msg("Reading labels CSV...")
labels    <- read.csv(LABELS, row.names = 1, check.names = FALSE)
label_col <- colnames(labels)[1]
log_msg(sprintf("Label column: %s", label_col))

common2 <- intersect(colnames(obj), rownames(labels))
obj$scanvi_pred <- NA_character_
obj$scanvi_pred[common2] <- as.character(labels[common2, label_col])
log_msg(sprintf("scanvi_pred assigned to %d / %d cells", length(common2), ncol(obj)))
log_msg("Top predictions:")
print(sort(table(obj$scanvi_pred), decreasing = TRUE)[1:10])

# ── Save ──────────────────────────────────────────────────────────────────
log_msg(sprintf("Saving to: %s", OUT_RDS))
saveRDS(obj, OUT_RDS)
log_msg("Done.")
