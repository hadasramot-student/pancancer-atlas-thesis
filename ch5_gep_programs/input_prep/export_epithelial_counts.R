library(Seurat)
library(Matrix)

base_dir <- "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/merged_by_query"
out_dir  <- "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/epithelial_cnmf_input"
dir.create(out_dir, showWarnings = FALSE)

# Use _scanvi_final where available, else _merged
rds_files <- c(
  Colorectal   = file.path(base_dir, "Colorectal/Epithelial_Colorectal_merged.RDS"),
  Head_and_Neck= file.path(base_dir, "Head_and_Neck/Epithelial_Head_and_Neck_merged.RDS"),
  Kidney       = file.path(base_dir, "Kidney/Epithelial_Kidney_merged_full_scanvi_final.RDS"),
  Liver_Biliary= file.path(base_dir, "Liver_Biliary/Epithelial_Liver_Biliary_merged_full_scanvi_final.RDS"),
  Lung         = file.path(base_dir, "Lung/Epithelial_Lung_merged_full_scanvi_final.RDS"),
  Ovarian      = file.path(base_dir, "Ovarian/Epithelial_Ovarian_merged_full_scanvi_final.RDS"),
  Pancreas     = file.path(base_dir, "Pancreas/Epithelial_Pancreas_merged_full_scanvi_final.RDS"),
  Prostate     = file.path(base_dir, "Prostate/Epithelial_Prostate_merged_full_scanvi_final.RDS"),
  Sarcoma      = file.path(base_dir, "Sarcoma/Epithelial_Sarcoma_merged_full_scanvi_final.RDS"),
  Skin         = file.path(base_dir, "Skin/Epithelial_Skin_merged_full_scanvi_final.RDS")
)

all_counts <- list()
all_meta   <- list()

for (cancer_type in names(rds_files)) {
  path <- rds_files[[cancer_type]]
  message("Loading: ", cancer_type, " (", basename(path), ")")
  obj <- readRDS(path)

  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  message("  Cells: ", ncol(counts), "  Genes: ", nrow(counts))

  # Verify these are raw integer counts
  sample_vals <- counts@x[1:min(100, length(counts@x))]
  if (!all(sample_vals == floor(sample_vals))) {
    warning("  WARNING: counts slot does not appear to contain integers for ", cancer_type)
  } else {
    message("  Counts look like raw integers: OK")
  }

  meta <- data.frame(
    cell_barcode = colnames(counts),
    cancer_type  = cancer_type,
    row.names    = colnames(counts)
  )
  # Carry over useful obs columns if present
  for (col in c("SRX_accession", "orig.ident", "seurat_clusters")) {
    if (col %in% colnames(obj@meta.data)) {
      meta[[col]] <- obj@meta.data[colnames(counts), col]
    }
  }

  all_counts[[cancer_type]] <- counts
  all_meta[[cancer_type]]   <- meta
  rm(obj); gc()
}

# Find common genes across all cancer types
message("\nFinding common genes...")
gene_lists <- lapply(all_counts, rownames)
common_genes <- Reduce(intersect, gene_lists)
message("Common genes: ", length(common_genes))

# Subset to common genes and cbind
message("Merging count matrices...")
merged_counts <- do.call(cbind, lapply(all_counts, function(m) m[common_genes, ]))
merged_meta   <- do.call(rbind, all_meta)

message("Final merged matrix: ", nrow(merged_counts), " genes x ", ncol(merged_counts), " cells")

# Write outputs
message("Writing to disk...")
writeMM(merged_counts, file.path(out_dir, "counts.mtx"))
write.table(data.frame(gene = rownames(merged_counts)),
            file.path(out_dir, "genes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(data.frame(barcode = colnames(merged_counts)),
            file.path(out_dir, "barcodes.tsv"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.csv(merged_meta, file.path(out_dir, "cell_metadata.csv"), row.names = FALSE)

message("\nDone! Files written to: ", out_dir)
message("  counts.mtx      - sparse count matrix (genes x cells)")
message("  genes.tsv       - gene names")
message("  barcodes.tsv    - cell barcodes")
message("  cell_metadata.csv - cell metadata (cancer_type, sample, etc.)")
