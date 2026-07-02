# merge_by_cell_type_from_all_cells.R
#
# Ch.4 (Subclustering) — lineage-splitting step. Streams through the
# atlas-wide all_cells metadata table (built by rbind_tsvs.R) and, for each
# distinct Consensus_Cell_Type (Ch.2 consensus vote), pulls the matching
# cells out of each per-object *_step4.RDS and merges them into one
# lineage-wide Seurat object — the input to Ch.4 clustering
# (clustering_jobs/run_clustering.R).
#
# Inputs:  ATLAS_DATA_DIR/all_cells.RDS (atlas-wide metadata, see usage note
#          below), ATLAS_DATA_DIR/<project>_step4.RDS per source object
# Outputs: ATLAS_DATA_DIR/merged_by_cell_type/<Consensus_Cell_Type>_merged.RDS
#
# Usage:
#   ATLAS_DATA_DIR=/path/to/temp_objects Rscript merge_by_cell_type_from_all_cells.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

merge_by_celltype_streaming <- function(all_cells,
                                        input_dir,
                                        output_dir,
                                        celltype_col = "Consensus_Cell_Type",
                                        project_col  = "project_name",
                                        cellid_col   = c("cell_name", "cell_id"),
                                        file_suffix  = "_step4.RDS",
                                        assay        = NULL,
                                        diet         = TRUE,
                                        keep_layers  = c("counts"),
                                        verbose      = TRUE) {
  
  # dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # pick first existing cell-id column
  cellid_col <- cellid_col[cellid_col %in% names(all_cells)]
  if (length(cellid_col) == 0) stop("None of cellid_col exist in all_cells (tried: cell_name, cell_id).")
  cellid_col <- cellid_col[2]
  
  req <- c(celltype_col, project_col, cellid_col)
  missing <- setdiff(req, names(all_cells))
  if (length(missing) > 0) stop("Missing required columns in all_cells: ", paste(missing, collapse = ", "))
  
  meta <- all_cells %>%
    transmute(
      celltype = as.character(.data[[celltype_col]]),
      project  = as.character(.data[[project_col]]),
      cell     = as.character(.data[[cellid_col]])
    ) %>%
    filter(!is.na(celltype), celltype != "",
           !is.na(project), project != "",
           !is.na(cell), cell != "")
  
  celltypes <- sort(unique(meta$celltype))
  
  read_obj <- function(path) {
    tryCatch(readRDS(path), error = function(e) {
      message("Failed to read: ", path, " | ", conditionMessage(e))
      NULL
    })
  }
  
  diet_obj <- function(obj, assay = NULL) {
    if (!diet) return(obj)
    
    if (!is.null(assay) && assay %in% Seurat::Assays(obj)) {
      Seurat::DefaultAssay(obj) <- assay
    }
    
    tryCatch(
      Seurat::DietSeurat(
        obj,
        assays = Seurat::DefaultAssay(obj),
        counts = "counts" %in% keep_layers,
        data   = "data"   %in% keep_layers,
        scale.data = "scale.data" %in% keep_layers,
        dimreducs = NULL,
        graphs = NULL
      ),
      error = function(e) obj
    )
  }
  
  # ct <- "Fibroblasts"
  for (ct in celltypes) {
    print(paste0("Starting: ", ct))
    ct_meta <- meta %>% filter(celltype == ct)
    proj_list <- split(ct_meta$cell, ct_meta$project)
    
    if (verbose) message("\n=== Cell type: ", ct, " | projects: ", length(proj_list), " ===")
    
    merged_ct <- NULL
    merged_any <- FALSE
    
    for (proj in names(proj_list)) {
      
      rds_path <- file.path(input_dir, paste0(proj, file_suffix))
      if (!file.exists(rds_path)) {
        if (verbose) message("  [skip] missing file: ", rds_path)
        next
      }
      
      cells_keep <- unique(proj_list[[proj]])
      
      obj <- read_obj(rds_path)
      if (is.null(obj)) next
      print(paste0("Loaded: ",rds_path ))
      
      cells_in_obj <- intersect(colnames(obj), cells_keep)
      if (length(cells_in_obj) == 0) {
        if (verbose) message("  [skip] no matching cells in ", proj)
        rm(obj); gc()
        next
      }
      
      obj_sub <- subset(obj, cells = cells_in_obj)
      obj_sub <- diet_obj(obj_sub)
      obj_sub$source_project <- proj
      
      # ensure unique cell names across projects
      obj_sub <- RenameCells(obj_sub, add.cell.id = proj)
      
      if (!merged_any) {
        merged_ct <- obj_sub
        merged_any <- TRUE
        if (verbose) message("  * init with ", proj, ": ", ncol(merged_ct), " cells")
      } else {
        merged_ct <- merge(x = merged_ct, y = obj_sub)
        if (verbose) message("  + merged ", proj, " | total now: ", ncol(merged_ct))
      }
      
      # free memory aggressively
      rm(obj, obj_sub)
      gc()
    }
    
    if (!merged_any) {
      if (verbose) message("  [skip] nothing merged for cell type: ", ct)
      next
    }
    
    out_file <- file.path(output_dir, paste0(gsub("[^A-Za-z0-9_\\-]+", "_", ct), "_merged.RDS"))
    saveRDS(merged_ct, out_file)
    print(paste0("saved in: ",out_file ))
    
    if (verbose) message("  -> saved: ", out_file, " | total cells: ", ncol(merged_ct))
    
    rm(merged_ct)
    gc()
  }
  
  invisible(TRUE)
}

#usage:
# all_cells <-
#   dir(file.path(atlas_data_dir),
#       pattern = "*.tsv",
#       full.names = T) %>%
#   # head() %>%
#   purrr::map(.f = function(x) {
#     print(x)
#     read_tsv(x,show_col_types = F) %>%
#       mutate(patient = as.character("patient")) %>%
#       dplyr::select(-any_of("predicted_doublet"))
#   })
# all_cells <-
#   all_cells %>%
#   list_rbind()
# saveRDS(
#   all_cells,
#   file = file.path(atlas_data_dir, "all_cells.RDS")
# )

atlas_data_dir <- Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects")

all_cells <- readRDS(file.path(atlas_data_dir, "all_cells.RDS"))
print ("Loaded all_cells")
merge_by_celltype_streaming(
  all_cells   = all_cells,
  input_dir   = atlas_data_dir,
  output_dir  = file.path(atlas_data_dir, "merged_by_cell_type"),
  celltype_col = "Consensus_Cell_Type",   # or "Consensus_Cell_Type_sep"
  project_col  = "project_name",
  file_suffix  = "_step4.RDS",
  diet         = TRUE,
  keep_layers  = c("counts"),
  verbose      = TRUE
)

