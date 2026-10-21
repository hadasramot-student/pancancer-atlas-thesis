# 00_probe.R  --  read-only reconnaissance for the 200k scANVI Figure 1/2 build.
# No sampling, no modelling. Answers the open structural questions:
#   * metadata table: needed columns, join key, categorical levels
#   * project_name -> *_step4.RDS coverage (every atlas project must have a file)
#   * step4 RDS structure: class, barcodes == metadata cell_id, counts assay,
#     gene-name consistency across repos
options(width = 200, stringsAsFactors = FALSE)
TO   <- "/ofir_grp/users/hadasra/temp_objects"
OUT  <- file.path(TO, "figure1_2_scanvi", "work")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
lg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

DUPS <- c("WU_SZ_ET_AL_2021_P15_filtered","GRIFFITHS_JI_ET_AL_2021_P18_filtered",
          "GONZALEZ_H_ET_AL_2022_P31_filtered","SAVAS_P_ET_AL_2018_P35_filtered",
          "GAO_R_ET_AL_2021_P17_filtered")

## ---------------------------------------------------------------- metadata ----
f <- file.path(TO, "all_cells.gapfix_20260727_164657.RDS")
lg("reading ", f)
md <- as.data.frame(readRDS(f))
lg("raw dim: ", paste(dim(md), collapse = " x "))
md <- md[!(md$project_name %in% DUPS), , drop = FALSE]
lg("post-dedup dim: ", paste(dim(md), collapse = " x "))

need <- c("cell_id","project_name","cancer_type","Consensus_Cell_Type",
          "maj3_cell","maj2_cell","maj2_cluster","maj3_cluster",
          "scimilarity_celltype","SRX_accession")
lg("needed columns present: ",
   paste(need[need %in% colnames(md)], collapse = ", "))
lg("needed columns MISSING: ",
   paste(need[!need %in% colnames(md)], collapse = ", "))

cat("\n--- cell_id head ---\n");        print(utils::head(md$cell_id, 5))
cat("--- cell_id dup? ---\n");          print(anyDuplicated(md$cell_id))
cat("\n--- cancer_type ---\n");         print(sort(table(md$cancer_type), decreasing = TRUE))
cat("\n--- Consensus_Cell_Type ---\n"); print(sort(table(md$Consensus_Cell_Type), decreasing = TRUE))
cat("\n--- maj3_cell (top 20) ---\n");  print(utils::head(sort(table(md$maj3_cell,  useNA="ifany"), decreasing = TRUE), 20))
cat("\n--- maj2_cell (top 20) ---\n");  print(utils::head(sort(table(md$maj2_cell,  useNA="ifany"), decreasing = TRUE), 20))
cat("\n--- maj2_cluster (top 20) ---\n");print(utils::head(sort(table(md$maj2_cluster,useNA="ifany"),decreasing = TRUE), 20))
cat("\n--- maj3_cluster (top 20) ---\n");print(utils::head(sort(table(md$maj3_cluster,useNA="ifany"),decreasing = TRUE), 20))
cat("\n--- SRX_accession NA count / non-NA count ---\n")
print(c(NA_ = sum(is.na(md$SRX_accession)), nonNA = sum(!is.na(md$SRX_accession))))
cat("\n--- cell_name / cell_type non-NA (3CA marker) ---\n")
print(c(cell_type_nonNA = sum(!is.na(md$cell_type)),
        cell_name_nonNA = if ("cell_name" %in% colnames(md)) sum(!is.na(md$cell_name)) else NA))

## agreement-level cross-check: does maj3_cell imply Consensus, etc.
cat("\n--- agreement flag tabulation (spec logic) ---\n")
agr <- ifelse(!is.na(md$maj3_cell), "Full",
       ifelse(!is.na(md$maj2_cell) | !is.na(md$maj2_cluster), "Partial", NA))
agr[md$Consensus_Cell_Type == "other"] <- "Other/Unknown"
agr[is.na(agr)] <- "UNASSIGNED_check"
print(table(agr, useNA = "ifany"))
cat("(cross) agr vs Consensus==other:\n")
print(table(agr, is_other = md$Consensus_Cell_Type == "other"))

## ------------------------------------------------- project -> step4 file ------
proj <- as.data.frame(table(md$project_name), responseName = "n_atlas")
colnames(proj)[1] <- "project_name"
proj <- proj[order(-proj$n_atlas), ]
lg("n atlas projects (post-dedup): ", nrow(proj))
step4 <- list.files(TO, pattern = "_step4\\.RDS$", full.names = FALSE)
step4_stem <- sub("_step4\\.RDS$", "", step4)
proj$has_step4 <- proj$project_name %in% step4_stem
cat("\n--- atlas projects WITHOUT a *_step4.RDS (BLOCKER if any) ---\n")
print(proj[!proj$has_step4, ])
cat("\n--- *_step4.RDS not matched to any atlas project ---\n")
print(setdiff(step4_stem, proj$project_name))
write.csv(proj, file.path(OUT, "atlas_project_step4_coverage.csv"), row.names = FALSE)

## ------------------------------------------------- step4 RDS structure -------
suppressPackageStartupMessages(library(Seurat))
pick <- c(`3CA` = "Azizi2018_Breast")
sc_srp <- step4_stem[grepl("^SRP", step4_stem)][1]
bs_flt <- step4_stem[grepl("_filtered$", step4_stem)][1]
pick <- c(pick, ScBaseCount = sc_srp, Breastscape = bs_flt)
lg("inspecting step4 objects: ", paste(pick, collapse = ", "))

rn_list <- list()
for (nm in names(pick)) {
  p <- file.path(TO, paste0(pick[[nm]], "_step4.RDS"))
  cat("\n=================== ", nm, " : ", pick[[nm]], " ===================\n")
  if (!file.exists(p)) { cat("  MISSING FILE\n"); next }
  cat("  file size (GB): ", round(file.info(p)$size/1e9, 2), "\n")
  o <- readRDS(p)
  cat("  class: ", paste(class(o), collapse = ","), "\n")
  cat("  dim (genes x cells): ", paste(dim(o), collapse = " x "), "\n")
  cat("  colnames[1:3]: ", paste(utils::head(colnames(o), 3), collapse = " | "), "\n")
  cat("  rownames[1:6]: ", paste(utils::head(rownames(o), 6), collapse = " | "), "\n")
  if ("RNA" %in% Assays(o)) {
    a <- o[["RNA"]]
    cat("  RNA assay class: ", paste(class(a), collapse = ","), "\n")
    cat("  RNA layers: ", paste(tryCatch(SeuratObject::Layers(a), error=function(e) "n/a"), collapse = ", "), "\n")
    cm <- tryCatch(SeuratObject::GetAssayData(o, assay = "RNA", layer = "counts"),
                   error = function(e) tryCatch(SeuratObject::LayerData(a, "counts"), error=function(e2) NULL))
    if (!is.null(cm) && length(cm@x) > 0) {
      cat("  counts class: ", paste(class(cm), collapse = ","),
          " | nnz: ", length(cm@x),
          " | integer-valued: ", isTRUE(all(cm@x == floor(cm@x))),
          " | max: ", max(cm@x), "\n")
    } else cat("  counts layer: EMPTY or unreadable\n")
  }
  rn_list[[nm]] <- rownames(o)
  ## barcode <-> metadata cell_id cross-check
  ids_md <- md$cell_id[md$project_name == pick[[nm]]]
  cat("  n cells in obj: ", ncol(o), " | n cell_id in atlas md for this project: ", length(ids_md), "\n")
  cat("  setequal(colnames, atlas cell_id): ", setequal(colnames(o), ids_md), "\n")
  cat("  in obj not md (up to 3): ", paste(utils::head(setdiff(colnames(o), ids_md), 3), collapse = " | "), "\n")
  cat("  in md not obj (up to 3): ", paste(utils::head(setdiff(ids_md, colnames(o)), 3), collapse = " | "), "\n")
  rm(o); gc(verbose = FALSE)
}

cat("\n--- gene-name consistency across repos ---\n")
if (length(rn_list) >= 2) {
  nms <- names(rn_list)
  for (i in seq_along(rn_list)) cat(sprintf("  %s: %d genes\n", nms[i], length(rn_list[[i]])))
  common <- Reduce(intersect, rn_list)
  cat("  intersection size: ", length(common), "\n")
  cat("  identical(3CA, ScBaseCount): ", identical(rn_list[[1]], rn_list[[2]]), "\n")
  if (length(rn_list) >= 3)
    cat("  identical(3CA, Breastscape): ", identical(rn_list[[1]], rn_list[[3]]), "\n")
}
lg("DONE")
