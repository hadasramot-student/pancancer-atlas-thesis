# functions.R
#
# Shared helper library used across Ch.1 (ingestion/QC) and Ch.2
# (annotation/consensus). Key functions:
#   creating_object()          - QC + normalize + cluster a raw Seurat object (Ch.1)
#   singleR_to_metadata()      - SingleR annotation (Ch.2)
#   sctype_cell_annotation()   - SCType annotation (Ch.2)
#   add_meta_from_scimilarity()- map Scimilarity labels to the meta-type ontology (Ch.2)
#   original_de_dup()          - cross-validate against original 3CA labels (Ch.2)
#   into_Meta()                - map SingleR/SCType labels to the meta-type ontology (Ch.2)
#   add_consensuses()          - maj2/maj3/maj6 consensus voting -> Consensus_Cell_Type (Ch.2)
#   merge_cancer_cell_types()  - streaming merge of per-object cells by consensus type (Ch.4)
#
# Requires a lab-internal Utils.R/Utils.MS.R (generic R helpers) on LAB_UTILS_DIR.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(tidyverse)
  library(SingleR)
  library(celldex)
  library(openxlsx)
  library(fastmap)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(readxl)
  library(SummarizedExperiment)
  library(SingleCellExperiment)
  library(RCurl)  # Allows downloading scripts without saving locally
  # library(usethis)   # disabled: requires fs >= 2.1.0 via devtools
  # library(devtools)  # disabled: requires fs >= 2.1.0, only 2.0.1 available in container
  library(logger)
  library(HGNChelper)
})

lab_utils_dir <- Sys.getenv("LAB_UTILS_DIR", unset = "config/lab_utils")
source(file.path(lab_utils_dir, "Utils.MS.R"))
source(file.path(lab_utils_dir, "Utils.R"))

#sources for sctype
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")     # often needed next


get_project_name <- function(path) {
  file_no_ext <- tools::file_path_sans_ext(basename(path))
  # replace commas with underscores
  clean_name <- gsub(",", "_", file_no_ext)
  return(clean_name)
}


extract_3CA_data_from_row <- function(projects_df , rownum) {
  
  sub_df <- projects_df %>%
    filter(technology == "10x", Cancer_type != "Hematologic") %>%
    group_by(project_name) %>%
    mutate(rn= row_number(), n= n()) %>%
    mutate(project_name_i = glue::glue("{project_name}_{rn}")) %>%
    mutate(project_name_i = if_else(n==1, project_name, project_name_i))
  
  project_name <- sub_df$project_name_i[rownum]
  mtx_path <- sub_df$mtx_file[rownum]
  cells_path <- sub_df$path[rownum]
  cancer_type <- sub_df$Cancer_type[rownum]
  genes_path <- sub_df$gene_file[rownum]
  return(list(
    project_name = project_name,
    mtx_path     = mtx_path,
    cells_path   = cells_path,
    cancer_type  = cancer_type,
    genes_path   = genes_path
  ))
  
  
}


save_metadata_tsv <- function(obj, out_path) {
  if (data_source == "breastscape"){
    cancer_tp <- "Breast"
  } else{
    cancer_tp <- unique(obj@meta.data$cancer_type)
  } 
  
  meta <- obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(cancer_type = cancer_tp)
  
  write_tsv(meta, out_path)
}


save_metadata_tsv_2 <- function(obj, out_path) {
  cancer_tp <- "Breast"
  meta <- obj@meta.data %>%
    rownames_to_column("cell_id") %>%
    mutate(cancer_type = cancer_tp)
  
  write_tsv(meta, out_path)
}

rownames_as_gene_symbols <- function(seurat_object,
                                     lookup_table,
                                     assay = "RNA",
                                     ensembl_col = "ensembl_gene_id",
                                     symbol_col  = "external_gene_name") {
  stopifnot(inherits(seurat_object, "Seurat"))
  stopifnot(assay %in% names(seurat_object@assays))
  stopifnot(all(c("ensembl_gene_id", "external_gene_name") %in% colnames(lookup_table)))
  
  feats <- rownames(seurat_object[[assay]])
  if (is.null(feats) || length(feats) == 0) stop("No feature names found in this assay.")
  
  ens0 <- sub("\\..*$", "", feats)  # strip version suffix
  
  # build a named vector: ensembl_gene_id -> symbol
  lt <- lookup_table[, c("ensembl_gene_id", "external_gene_name")]
  lt <- lt[!is.na(lt$external_gene_name) & lt$external_gene_name != "", ]
  lt <- lt[!duplicated(lt$ensembl_gene_id), ]
  m  <- setNames(as.character(lt$external_gene_name), lt$ensembl_gene_id)
  
  # map: only replace ENSG-like ids; otherwise keep as-is
  sym <- ens0
  is_ensg <- grepl("^ENSG[0-9]+$", ens0)
  hit <- match(ens0[is_ensg], names(m))
  sym[is_ensg][!is.na(hit)] <- m[ens0[is_ensg][!is.na(hit)]]
  
  # ensure uniqueness
  sym <- make.unique(sym)
  
  # feature metadata table (Seurat v5; can be 0 cols)
  feat_md <- seurat_object[[assay]][[]]
  if (nrow(feat_md) != length(feats)) stop("Feature metadata rows != number of features.")
  rownames(feat_md) <- feats
  
  feat_md[[ensembl_col]] <- ens0
  feat_md[[symbol_col]]  <- sym
  seurat_object[[assay]][[]] <- feat_md
  
  # rename features in the assay (updates layers)
  dn <- dimnames(seurat_object[[assay]])
  dimnames(seurat_object[[assay]]) <- list(sym, dn[[2]])
  
  return(seurat_object)
}

creating_object <- function(mSeurat.obj, project_name = NULL, mtx_path = NULL, cells_path = NULL, genes_path = NULL) {
  
  SC.3CA.folder = Sys.getenv("SC3CA_DIR", unset = "data/raw/3ca")
  SC.params = list()
  if(is.null(mSeurat.obj)){
    SC.params$Project.name <- project_name
    print(SC.params$Project.name)
  }
  
  SC.params$min.cells = 3
  SC.params$min.features = 200
  SC.params$scale.factor = 10000
  SC.params$n.var.genes = 10000
  SC.params$n.PC = 20
  SC.params$is.save.obj = T
  SC.params$is.load.obj = F
  SC.params$clust.resolutions = c(0.25,0.1, 0.5, 0.75, 1, 2)
  
  if (!is.null(SC.params$sub_project_name) && SC.params$sub_project_name != ""){ #there is a sub p
    SC.params$Project.name <- paste0(SC.params$Project.name,"_", SC.params$sub_project_name)
  }
  if (length(Layers(seurat_object[["RNA"]])) != 3){                                                                                                 
    seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
  }
  mSeurat.obj <- NormalizeData(mSeurat.obj, scale.factor = SC.params$scale.factor)
  mSeurat.obj <- FindVariableFeatures(mSeurat.obj,  nfeatures = SC.params$n.var.genes)
  mSeurat.obj <- ScaleData(mSeurat.obj)
  mSeurat.obj <- RunPCA(mSeurat.obj, features = VariableFeatures(mSeurat.obj))
  mSeurat.obj <- FindNeighbors(mSeurat.obj, dims = 1:SC.params$n.PC, k.param = 5)
  
  for (i in SC.params$clust.resolutions){
    mSeurat.obj <- FindClusters(mSeurat.obj, resolution = i)  
  }
  mSeurat.obj <- RunUMAP(mSeurat.obj, dims = 1:10)
  
  while(TRUE) {
    open_device <- tryCatch({
      dev.off()
      TRUE  
    }, error = function(e) {
      FALSE  
    })
    if (!open_device) {
      break
    }
  }
  return(mSeurat.obj)
}


singleR_to_metadata <- function(seurat_object) {
  
  hpca.se <- celldex::HumanPrimaryCellAtlasData()
  project_name <- seurat_object@project.name
  print(paste0("loaded: ",seurat_object@project.name ))
  
  #for a merged object from ScBaseCount
  if (length(Layers(seurat_object[["RNA"]])) != 3){                                                                                                 
    seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
  }
  seurat_matrix <- GetAssayData(seurat_object, assay = "RNA", layer = "data")
  common_genes <- intersect(rownames(seurat_matrix), rownames(hpca.se))
  seurat_matrix <- seurat_matrix[common_genes, ]
  hpca.se <- hpca.se[common_genes, ]
  singleR_results <- SingleR(test = seurat_matrix, ref = hpca.se, labels = hpca.se$label.main)
  
  seurat_object@meta.data$SingleR.labels <- singleR_results$labels
  
  return(seurat_object)
}

original_de_dup <- function(seurat_object, df_path) {
  df <- read_excel(df_path)
  colnames(df) <- c("original", "original_de_dup")  # Ensure column names are correct
  
  meta_data <- seurat_object@meta.data
  project_name <- seurat_object@project.name
  print(paste0("Loaded: ", project_name))
  
  # Skip if project name contains a space
  if (grepl(" ", project_name)) {
    print(paste0("Skipping project: ", project_name, " (contains space)"))
    return(NULL)  # Return NULL to avoid errors when combining results
  }
  if (!"cell_type" %in% colnames(seurat_object@meta.data)) {
    stop("cell_type column not found in metadata")
  }
  
  seurat_object@meta.data <- seurat_object@meta.data %>%
    left_join(df, by = c("cell_type" = "original"))
  
  print(paste0("Done: ", project_name))
  return(seurat_object)
}


sctype_cell_annotation <- function(seurat_object, db, tissue) {
  
  print(c("loaded: ", seurat_object@project.name))
  project_name <- seurat_object@project.name
  if (grepl(" ", project_name)) {
    print(paste0("Skipping project: ", project_name, " (contains space)"))
    return(NULL)  # Return NULL to avoid errors when combining results
  }
  colnames(seurat_object@meta.data)[which(colnames(seurat_object@meta.data) == "sctype_classification")] <- "sctype_w_cluster"
  
  seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seurat_object[["RNA"]])));
  scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(seurat_object[["RNA"]]$scale.data) else as.matrix(seurat_object[["RNA"]]@scale.data)
  
  # prepare gene sets
  gs_list <- gene_sets_prepare(db, tissue)
  
  if (length(intersect(rownames(scRNAseqData_scaled), unique(unlist(gs_list$gs_negative)))) == 0){
    mat <- scRNAseqData_scaled  # rename for clarity; should be a matrix
    
    # (Optional) strip version if you had ENSG000001234.5 — here it’s probably a no-op:
    ens_ids <- sub("\\..*$", "", rownames(mat))
    
    ###check if the rownames are ensembl-like
    rn <- rownames(scRNAseqData_scaled)
    is_ensembl <- any(grepl("^ENSG|^ENSMUSG", rn))
    # also handle versioned Ensembl like ENSG... .5
    is_ensembl <- is_ensembl || any(grepl("^ENSG|^ENSMUSG", sub("\\..*$", "", rn)))
    
    if (is_ensembl) {
      symbols <- mapIds(
        org.Hs.eg.db,
        keys     = ens_ids,
        keytype  = "ENSEMBL",
        column   = "SYMBOL",
        multiVals = "first"
      )
      
      # Keep only genes that successfully mapped to a symbol
      keep <- !is.na(symbols)
      
      mat2      <- mat[keep, ]
      symbols2  <- symbols[keep]
      
      # Replace Ensembl IDs by gene symbols
      rownames(mat2) <- symbols2
      
      # Collapse duplicates (several Ensembl IDs → same symbol)
      mat2 <- rowsum(mat2, group = rownames(mat2))
      scRNAseqData_scaled <- mat2
    }else{
      message("Row names look like gene symbols already; skipping ENSEMBL->SYMBOL mapping.")
    }
  }
  # run ScType
  es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # Prepare the ScType scores for each cell
  es.max.cl <- rowSums(es.max)  # Sum the scores across cells for each gene set
  
  # ScType annotation for each cell without using clusters
  sctype_scores <- data.frame(
    cell = colnames(es.max),  # assuming columns are cells
    type = apply(es.max, 2, function(x) names(sort(x, decreasing = TRUE)[1])),  # Assign cell type based on max score
    scores = apply(es.max, 2, max)  # The score for the chosen type
  )
  
  # classify low confidence cells to "Unknown" based on some threshold
  threshold <- 0.25  
  sctype_scores$type[sctype_scores$scores < threshold] <- "Unknown"
  # Add the cell type annotations to the Seurat metadata
  seurat_object@meta.data$sctype <- sctype_scores$type
  
  return(seurat_object)
}

into_Meta <- function(seurat_object, path) {
  
  ontology <- openxlsx::read.xlsx(path)
  names(ontology) <- trimws(names(ontology)) ########new
  ontology <- ontology %>%
    dplyr::mutate(dplyr::across(where(is.character), ~ trimws(.))) %>%
    dplyr::mutate(Cell_type = stringr::str_replace(Cell_type, "^ +", ""))
  
  #ontology_orig    <- ontology %>% dplyr::filter(Table == "Original") %>% dplyr::select(-Table)
  ontology_singler <- ontology %>% dplyr::filter(Table == "SingleR") %>% dplyr::select(-Table)
  ontology_sctype  <- ontology %>% dplyr::filter(Table == "Sctype_Immune_plus") %>% dplyr::select(-Table)
  
  meta_table <- seurat_object@meta.data
  project_name <- seurat_object@project.name
  print(paste0("Loaded: ", project_name))
  
  # Skip if project name contains a space
  if (grepl(" ", project_name)) {
    print(paste0("Skipping project: ", project_name, " (contains space)"))
    return(NULL)  # Return NULL to avoid errors when combining results
  }
  meta_table <- seurat_object@meta.data %>%
    rownames_to_column("cell_id")   ## new: keep cell names in a column
  
  meta_updated <- meta_table %>% 
    left_join(ontology_singler, by = join_by(SingleR.labels == Cell_type)) %>% 
    dplyr::rename(Meta.SingleR = Meta) %>% 
    left_join(ontology_sctype,  by = join_by(sctype == Cell_type)) %>% 
    dplyr::rename(Meta.Sctype = Meta) %>% 
    mutate(
      maj2_cell = case_when(
        Meta.Sctype == Meta.SingleR ~ Meta.Sctype,
        TRUE ~ "equivocal"
      )
    ) %>%
    column_to_rownames("cell_id")   # restore cell IDs as rownames
  
  seurat_object@meta.data <- meta_updated
  
  return(seurat_object)
  
}

add_meta_from_scimilarity <- function(
    obj,
    xlsx_path,
    sheet = 2,                       # your file has Sheet2 with the Scimilarity mapping
    seurat_scim_col = "scimilarity_celltype",
    out_col = "Meta.scimilarity",
    excel_meta_col = 1,              # "1st column" in the xlsx sheet
    excel_scim_col = "Scimilarity",  # the column name in the xlsx
    unknown_value = "unknown",
    ignore_case = TRUE
) {
  stopifnot(inherits(obj, "Seurat"))
  if (!seurat_scim_col %in% colnames(obj@meta.data)) {
    stop(sprintf("Column '%s' not found in obj@meta.data", seurat_scim_col))
  }
  
  # packages
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Please install readxl")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install tidyr")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Please install stringr")
  
  df <- readxl::read_excel(xlsx_path, sheet = 2)
  
  meta_name <- colnames(df)[excel_meta_col]
  if (!excel_scim_col %in% colnames(df)) {
    stop(sprintf("Column '%s' not found in Excel sheet. Columns are: %s",
                 excel_scim_col, paste(colnames(df), collapse = ", ")))
  }
  
  # Build a (scimilarity_label -> meta) mapping
  map_long <- df |>
    dplyr::select(meta = dplyr::all_of(meta_name), scim = dplyr::all_of(excel_scim_col)) |>
    dplyr::filter(!is.na(meta), !is.na(scim)) |>
    dplyr::mutate(scim = as.character(scim)) |>
    tidyr::separate_rows(scim, sep = ";") |>
    dplyr::mutate(
      scim = stringr::str_trim(scim),
      scim_key = if (ignore_case) tolower(scim) else scim
    ) |>
    dplyr::filter(scim_key != "")
  
  # keep first mapping if duplicates exist
  map_vec <- stats::setNames(map_long$meta, map_long$scim_key)
  
  # Apply mapping to Seurat meta.data
  scim_values <- as.character(obj@meta.data[[seurat_scim_col]])
  scim_key <- if (ignore_case) tolower(stringr::str_trim(scim_values)) else stringr::str_trim(scim_values)
  
  mapped <- unname(map_vec[scim_key])
  mapped[is.na(mapped)] <- unknown_value
  
  obj@meta.data[[out_col]] <- mapped
  return(obj)
}


add_precentage_and_project_name <- function(seurat_object) {
  
  print(paste0("loaded: ",seurat_object@project.name ))
  
  # Add project name to each cell name
  project_name <- seurat_object@project.name
  cell_names <- colnames(seurat_object)
  meta_data <- seurat_object@meta.data
  new_cell_names <- paste0(project_name, "_", cell_names)
  colnames(seurat_object) <- new_cell_names
  seurat_object$cell_name <- new_cell_names
  
  #adding the majority_cell_type and majority_percentage:
  
  majority_cell_cluster <- meta_data %>%
    group_by(seurat_clusters) %>%
    summarise(
      majority_cell_type = names(which.max(table(cell_type))),
      majority_percentage = max(table(cell_type)) / sum(table(cell_type)) * 100
    ) %>% 
    mutate(project_name = seurat_object@project.name)
  
  seurat_object@meta.data$majority_cell_type <- majority_cell_cluster$majority_cell_type[match(meta_data$seurat_clusters, majority_cell_cluster$seurat_clusters)]
  seurat_object@meta.data$majority_percentage <- majority_cell_cluster$majority_percentage[match(meta_data$seurat_clusters, majority_cell_cluster$seurat_clusters)]
  
  return (seurat_object)
}

annotate_majority_singler <- function(
    obj,
    singler_col = "Meta.SingleR",
    cluster_col = NULL,
    out_label_col = "majority_cell_type_SingleR",
    out_pct_col = "majority_percentage_SingleR"
    #digits = 1
) {
  # Basic checks
  stopifnot("meta.data" %in% slotNames(obj))
  md <- obj@meta.data
  if (!singler_col %in% colnames(md)) {
    stop(sprintf("Column '%s' not found in meta.data.", singler_col))
  }
  # Determine cluster vector
  if (is.null(cluster_col)) {
    clust <- as.character(Seurat::Idents(obj))
    if (is.null(clust) || length(clust) != nrow(md)) {
      stop("Idents(obj) is NULL or length mismatch. Provide cluster_col instead.")
    }
  } else {
    if (!cluster_col %in% colnames(md)) {
      stop(sprintf("Cluster column '%s' not found in meta.data.", cluster_col))
    }
    clust <- as.character(md[[cluster_col]])
  }
  
  labels <- as.character(md[[singler_col]])
  # Build per-cluster summaries (base R to avoid tidyverse dependency)
  df <- data.frame(cluster = clust, label = labels, stringsAsFactors = FALSE)
  # Remove NAs from label when computing majority (optional)
  df_nona <- df[!is.na(df$label), , drop = FALSE]
  
  # Tabulate counts per cluster x label
  # Split by cluster, then find majority label and pct for each cluster
  split_list <- split(df_nona$label, df_nona$cluster)
  majority_info <- lapply(names(split_list), function(k) {
    v <- split_list[[k]]
    tb <- sort(table(v), decreasing = TRUE)
    top_label <- names(tb)[1]
    top_count <- as.integer(tb[1])
    total <- length(v)
    pct <- if (total > 0) 100 * top_count / total else NA_real_
    data.frame(cluster = k,
               majority_label = top_label,
               majority_pct = round(pct, digits = 2),
               stringsAsFactors = FALSE)
  })
  majority_df <- do.call(rbind, majority_info)
  
  # Map back to all cells (including those with NA label originally)
  # Create lookup named by cluster
  lab_map <- setNames(majority_df$majority_label, majority_df$cluster)
  pct_map <- setNames(majority_df$majority_pct,   majority_df$cluster)
  
  md[[out_label_col]] <- unname(lab_map[clust])
  md[[out_pct_col]]   <- unname(pct_map[clust])
  
  obj@meta.data <- md
  return(obj)
}


add_precentage_sctype <- function(seurat_object) {
  
  meta_data <- seurat_object@meta.data
  project_name <- seurat_object@project.name
  print(paste0("Loaded: ", project_name))
  
  # Skip if project name contains a space
  if (grepl(" ", project_name)) {
    print(paste0("Skipping project: ", project_name, " (contains space)"))
    return(NULL)  # Return NULL to avoid errors when combining results
  }
  
  # Compute the majority classification and percentage in each cluster
  majority_cell_cluster <- meta_data %>%
    group_by(seurat_clusters) %>%
    summarise(
      majority_sctype = names(which.max(table(Meta.Sctype))),  # Most frequent cell type 
      total_cells = n(),  # Total cells in the cluster
      max_count = max(table(Meta.Sctype)),  # Count of the most common cell type
      majority_sctype_percentage = (max_count / total_cells) * 100  # Compute percentage
    ) %>%
    ungroup()
  
  # Merge back into metadata: Assign each cell the majority class of its cluster
  meta_data <- meta_data %>%
    left_join(majority_cell_cluster, by = "seurat_clusters")
  
  # Add to Seurat object metadata
  seurat_object@meta.data$majority_sctype <- meta_data$majority_sctype
  seurat_object@meta.data$majority_sctype_percentage <- meta_data$majority_sctype_percentage
  seurat_object@meta.data$project_name <- project_name
  meta_data <- seurat_object@meta.data
  
  return(seurat_object)
}


add_precentage_scimilarity <- function(seurat_object) {
  
  meta_data <- seurat_object@meta.data
  project_name <- seurat_object@project.name
  print(paste0("Loaded: ", project_name))
  
  # Skip if project name contains a space
  if (grepl(" ", project_name)) {
    print(paste0("Skipping project: ", project_name, " (contains space)"))
    return(NULL)  # Return NULL to avoid errors when combining results
  }
  
  # Compute the majority classification and percentage in each cluster
  majority_cell_cluster <- meta_data %>%
    group_by(seurat_clusters) %>%
    summarise(
      majority_scimilarity = names(which.max(table(Meta.scimilarity))),  # Most frequent cell type 
      total_cells = n(),  # Total cells in the cluster
      max_count = max(table(Meta.scimilarity)),  # Count of the most common cell type
      majority_scimilarity_percentage = (max_count / total_cells) * 100  # Compute percentage
    ) %>%
    ungroup()
  
  # Merge back into metadata: Assign each cell the majority class of its cluster
  meta_data <- meta_data %>%
    left_join(majority_cell_cluster, by = "seurat_clusters")
  
  # Add to Seurat object metadata
  seurat_object@meta.data$majority_scimilarity <- meta_data$majority_scimilarity
  seurat_object@meta.data$majority_scimilarity_percentage <- meta_data$majority_scimilarity_percentage
  seurat_object@meta.data$project_name <- project_name
  meta_data <- seurat_object@meta.data
  
  return(seurat_object)
}


add_sample_for_project <- function(
    obj,
    cells_path,
    meta_cell_col = "cell_name",
    file_cell_col = "cell_name",
    project_col = "project_name",
    sample_col = "sample",
    id_normalizer = NULL
){
  project <- obj@project.name
  cells_df <- read_csv(cells_path)
  
  # make sure we have the column in meta; if not, use rownames
  if (!(meta_cell_col %in% names(obj@meta.data))) {
    obj@meta.data[[meta_cell_col]] <- rownames(obj@meta.data)
  }
  
  # flatten list-cols and coerce to character
  to_chr <- function(x) {
    if (is.list(x)) {
      vapply(x, function(z) if (length(z)) as.character(z[[1]]) else NA_character_, character(1))
    } else as.character(x)
  }
  
  meta_ids_raw <- trimws(to_chr(obj@meta.data[[meta_cell_col]]))
  file_ids     <- trimws(to_chr(cells_df[[file_cell_col]]))
  file_samples <- trimws(to_chr(cells_df[[sample_col]]))
  
  # strip "<project_name>_" from meta cell_name if present
  prefix <- paste0(project, "_")
  meta_base <- ifelse(startsWith(meta_ids_raw, prefix),
                      substring(meta_ids_raw, nchar(prefix) + 1),
                      meta_ids_raw)
  
  # exact lookup
  lookup_exact <- setNames(file_samples, file_ids)
  out <- unname(lookup_exact[meta_base])
  
  obj[[sample_col]] <- out
  
  return(obj)
  
}

add_sample_for_project_h5ad <- function(seurat_object){
  if (seurat_object@meta.data$orig.ident %>% unique() %>% length() == 1 && seurat_object@meta.data$orig.ident %>% unique()!= "SeuratProject"){
    seurat_object@meta.data$sample <- seurat_object@meta.data$orig.ident
  }
  else{
    seurat_object@meta.data$sample <- seurat_object@meta.data$SRX_accession
  }
  return(seurat_object)
}

add_consensuses <- function(seurat_object, path) {
  
  ontology <- openxlsx::read.xlsx(path)
  names(ontology) <- trimws(names(ontology)) ########new
  ontology <- ontology %>%
    dplyr::mutate(dplyr::across(where(is.character), ~ trimws(.))) %>%
    dplyr::mutate(Cell_type = stringr::str_replace(Cell_type, "^ +", ""))
  
  ontology_scimilarity    <- ontology %>% dplyr::filter(Table == "Scimilarity") %>% dplyr::select(-Table)
  # ontology_orig    <- ontology %>% dplyr::filter(Table == "Original") %>% dplyr::select(-Table)
  ontology_singler <- ontology %>% dplyr::filter(Table == "SingleR") %>% dplyr::select(-Table)
  ontology_sctype  <- ontology %>% dplyr::filter(Table == "Sctype_Immune_plus") %>% dplyr::select(-Table)
  
  meta_table <- seurat_object@meta.data
  project_name <- seurat_object@project.name
  print(paste0("Loaded: ", project_name))
  
  # Skip if project name contains a space
  if (grepl(" ", project_name)) {
    print(paste0("Skipping project: ", project_name, " (contains space)"))
    return(NULL)  # Return NULL to avoid errors when combining results
  }
  
  meta_updated <- meta_table %>% 
    dplyr::mutate(
      maj3_cell = dplyr::if_else(Meta.scimilarity == Meta.Sctype & Meta.scimilarity == Meta.SingleR, Meta.scimilarity, "equivocal"),
      maj2_cell = dplyr::case_when(
        Meta.scimilarity == Meta.SingleR ~ Meta.scimilarity,
        Meta.scimilarity == Meta.Sctype  ~ Meta.scimilarity,
        Meta.Sctype == Meta.SingleR ~ Meta.Sctype,
        TRUE ~ "equivocal"
      ),
      maj3_cluster = dplyr::if_else(majority_scimilarity == majority_sctype & majority_scimilarity == majority_cell_type_SingleR, majority_scimilarity, "equivocal"),
      maj2_cluster = dplyr::case_when(
        majority_scimilarity == majority_cell_type_SingleR ~ majority_scimilarity,
        majority_scimilarity == majority_sctype  ~ majority_scimilarity,
        majority_sctype == majority_cell_type_SingleR ~ majority_sctype,
        TRUE ~ "equivocal"
      ),
      maj6 = dplyr::if_else(
        majority_scimilarity == majority_sctype & majority_scimilarity == majority_cell_type_SingleR &
          Meta.scimilarity == Meta.Sctype & Meta.scimilarity == Meta.SingleR,
        majority_scimilarity, "equivocal"
      )
    )
  meta_updated$Consensus_Cell_Type <- ifelse(meta_updated$maj2_cluster == meta_updated$maj2_cell, meta_updated$maj2_cell, "equivocal")
  meta_updated$Consensus_Cell_Type[meta_updated$Consensus_Cell_Type == "equivocal"] <- "other"
  seurat_object@meta.data <- meta_updated
  return(seurat_object)
}


#using at: merge_by_cancertype_celltype
# merge_cancer_cell_types <- function(proj_list, input_dir, file_suffix) {
#   
#   diet_obj <- function(obj) {
#     tryCatch(
#       obj <- DietSeurat(
#         object    = obj,
#         assays    = "RNA",
#         layers    = c("counts", "data"),
#         dimreducs = NULL,
#         graphs    = NULL,
#         misc      = FALSE
#       ),
#       error = function(e) obj
#     )
#     return(obj)
#   }
#   
#   seurat_list <- list()
#   
#   for (proj in names(proj_list)) {
#     rds_path <- file.path(input_dir, paste0(proj, file_suffix))
#     print(paste0("loaded: ", proj))
#     
#     if (!file.exists(rds_path)) next
#     
#     obj <- readRDS(rds_path)
#     cells_keep   <- unique(proj_list[[proj]])
#     cells_in_obj <- intersect(colnames(obj), cells_keep)
#     if (length(cells_in_obj) == 0) next
#     
#     obj_sub <- obj[, cells_in_obj]
#     obj_sub <- diet_obj(obj_sub)
#     
#     obj_sub$project_name <- proj
#     colnames(obj_sub) <- paste0(proj, "_", colnames(obj_sub))
#     
#     seurat_list[[proj]] <- obj_sub
#   }
#   
#   if (length(seurat_list) == 0) {
#     stop("No Seurat objects to merge; check meta_df and file paths.")
#   }
#   
#   merged <- Reduce(function(x, y) merge(x, y), seurat_list)
#   rm(seurat_list)
#   gc()
#   
#   return(merged)
# }

# merge_cancer_cell_types <- function(proj_list, input_dir, file_suffix) {
#   
#   diet_obj <- function(obj) {
#     tryCatch(
#       DietSeurat(
#         object    = obj,
#         assays    = "RNA",
#         layers    = c("counts", "data"),
#         dimreducs = NULL,
#         graphs    = NULL,
#         misc      = FALSE
#       ),
#       error = function(e) obj
#     )
#   }
#   
#   merged <- NULL
#   n_loaded <- 0L
#   
#   for (proj in names(proj_list)) {
#     rds_path <- file.path(input_dir, paste0(proj, file_suffix))
#     message("loaded: ", proj)
#     
#     if (!file.exists(rds_path)) next
#     
#     obj <- readRDS(rds_path)
#     
#     cells_keep <- unique(proj_list[[proj]])
#     cells_in_obj <- intersect(colnames(obj), cells_keep)
#     if (length(cells_in_obj) == 0) {
#       rm(obj)
#       gc()
#       next
#     }
#     
#     obj_sub <- obj[, cells_in_obj]
#     rm(obj)
#     gc()
#     
#     obj_sub <- diet_obj(obj_sub)
#     
#     obj_sub$project_name <- proj
#     colnames(obj_sub) <- paste0(proj, "_", colnames(obj_sub))
#     
#     if (is.null(merged)) {
#       merged <- obj_sub
#     } else {
#       merged <- merge(merged, y = obj_sub)
#     }
#     
#     rm(obj_sub)
#     gc()
#     n_loaded <- n_loaded + 1L
#   }
#   
#   if (is.null(merged)) {
#     stop("No Seurat objects to merge; check meta_df and file paths.")
#   }
#   
#   merged
# }

#this function with saving every rds after diet:
merge_cancer_cell_types <- function(proj_list, input_dir, file_suffix,
                                    max_cells_per_project = Inf,
                                    chunk_size = 10L,
                                    tmp_dir = tempdir()) {

  dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

  # Helper: memory usage in GB
  mem_gb <- function() {
    mem <- gc(verbose = FALSE)
    used_mb <- sum(mem[, 2])  # column 2 = used (Mb)
    sprintf("%.2f GB", used_mb / 1024)
  }
  
  diet_obj <- function(obj) {
    # keep only the RNA assay
    obj <- tryCatch(
      DietSeurat(object = obj, assays = "RNA",
                 dimreducs = NULL, graphs = NULL, misc = FALSE),
      error = function(e) obj
    )
    
    # explicitly drop every non-counts layer on RNA
    lyrs <- Layers(obj[["RNA"]])
    drop <- setdiff(lyrs, grep("^counts", lyrs, value = TRUE))
    for (ly in drop) {
      LayerData(obj[["RNA"]], layer = ly) <- NULL
    }
    obj
  }
  
  # ── Pass 1: read, subset, diet, subsample → save temp RDS ──────────────────
  tmp_files <- c()
  n_projects <- length(names(proj_list))
  
  log_msg(sprintf("=== Pass 1: processing %d projects === [mem: %s]",
                  n_projects, mem_gb()))
  
  for (idx in seq_along(names(proj_list))) {
    proj     <- names(proj_list)[idx]
    rds_path <- file.path(input_dir, paste0(proj, file_suffix))
    
    log_msg(sprintf("[%d/%d] Project: %s [mem: %s]",
                    idx, n_projects, proj, mem_gb()))
    
    if (!file.exists(rds_path)) {
      log_msg(sprintf("  SKIP (missing file): %s", rds_path))
      next
    }
    tmp_path <- file.path(tmp_dir, paste0("tmp_", proj, ".RDS"))
    
    # Skip if already processed
    if (file.exists(tmp_path)) {
      log_msg(sprintf("[%d/%d] SKIP (already exists): %s [mem: %s]",
                      idx, n_projects, proj, mem_gb()))
      tmp_files <- c(tmp_files, tmp_path)
      next
    }
    
    log_msg(sprintf("  Reading RDS: %s", rds_path))
    obj <- readRDS(rds_path)
    log_msg(sprintf("  Loaded: %d cells x %d genes [mem: %s]",
                    ncol(obj), nrow(obj), mem_gb()))
    
    cells_keep   <- unique(proj_list[[proj]])
    cells_in_obj <- intersect(colnames(obj), cells_keep)
    
    log_msg(sprintf("  Cells requested: %d | Found in object: %d",
                    length(cells_keep), length(cells_in_obj)))
    
    if (length(cells_in_obj) == 0) {
      log_msg(sprintf("  SKIP (no matching cells): %s", proj))
      rm(obj); gc()
      log_msg(sprintf("  After cleanup [mem: %s]", mem_gb()))
      next
    }
    
    obj_sub <- obj[, cells_in_obj]
    rm(obj); gc()
    log_msg(sprintf("  After subsetting to %d cells [mem: %s]",
                    ncol(obj_sub), mem_gb()))
    
    # # optional per-project cell cap
    # if (is.finite(max_cells_per_project) && ncol(obj_sub) > max_cells_per_project) {
    #   set.seed(42)
    #   keep    <- sample(colnames(obj_sub), max_cells_per_project)
    #   obj_sub <- obj_sub[, keep]
    #   log_msg(sprintf("  Downsampled to cap: %d cells [mem: %s]",
    #                   max_cells_per_project, mem_gb()))
    # }
    
    # # Hard cap at 100,000 cells per project
    # if (ncol(obj_sub) > 100000) {
    #   n_before <- ncol(obj_sub)
    #   set.seed(42)
    #   keep    <- sample(colnames(obj_sub), 100000)
    #   obj_sub <- obj_sub[, keep]
    #   log_msg(sprintf("  Capped %s: %d → 100,000 cells [mem: %s]",
    #                   proj, n_before, mem_gb()))
    # }
    
    log_msg(sprintf("  Running DietSeurat on %s", proj))
    obj_sub <- diet_obj(obj_sub)
    log_msg(sprintf("  After DietSeurat [mem: %s]", mem_gb()))
    
    obj_sub$project_name <- proj
    colnames(obj_sub)    <- paste0(proj, "_", colnames(obj_sub))
    
    tmp_path <- file.path(tmp_dir, paste0("tmp_", proj, ".RDS"))
    log_msg(sprintf("  Saving temp RDS: %s", tmp_path))
    saveRDS(obj_sub, tmp_path)
    
    tmp_files <- c(tmp_files, tmp_path)
    log_msg(sprintf("  Done: %s | %d cells saved [mem: %s]",
                    proj, ncol(obj_sub), mem_gb()))
    
    rm(obj_sub); gc()
    log_msg(sprintf("  After cleanup [mem: %s]", mem_gb()))
    log_msg(sprintf("  ----------------------------------------"))
  }
  
  if (length(tmp_files) == 0) {
    stop("No Seurat objects to merge; check proj_list and file paths.")
  }
  
  log_msg(sprintf("=== Pass 1 complete: %d/%d projects saved === [mem: %s]",
                  length(tmp_files), n_projects, mem_gb()))
  
  # ── Pass 2: chunked merge of temp files ────────────────────────────────────
  merge_file_list <- function(paths) {
    log_msg(sprintf("    Loading %d objects from disk", length(paths)))
    objs <- lapply(paths, readRDS)
    log_msg(sprintf("    Loaded %d objects [mem: %s]", length(objs), mem_gb()))
    
    if (length(objs) == 1) return(objs[[1]])
    
    merged <- merge(objs[[1]], y = objs[-1])
    log_msg(sprintf("    Merged → %d cells [mem: %s]", ncol(merged), mem_gb()))
    rm(objs); gc()
    log_msg(sprintf("    After cleanup [mem: %s]", mem_gb()))
    merged
  }
  
  chunks <- split(tmp_files, ceiling(seq_along(tmp_files) / chunk_size))
  log_msg(sprintf("=== Pass 2: merging %d projects in %d chunks of ≤%d === [mem: %s]",
                  length(tmp_files), length(chunks), chunk_size, mem_gb()))
  
  chunk_files <- c()
  for (i in seq_along(chunks)) {
    proj_names_in_chunk <- basename(chunks[[i]])
    proj_names_in_chunk <- gsub("^tmp_|.RDS$", "", proj_names_in_chunk)
    
    log_msg(sprintf("  Chunk %d/%d: projects = [%s] [mem: %s]",
                    i, length(chunks),
                    paste(proj_names_in_chunk, collapse=", "),
                    mem_gb()))
    
    chunk_merged <- merge_file_list(chunks[[i]])
    
    log_msg(sprintf("  Chunk %d result: %d cells x %d genes [mem: %s]",
                    i, ncol(chunk_merged), nrow(chunk_merged), mem_gb()))
    
    chunk_path <- file.path(tmp_dir, sprintf("chunk_%03d.RDS", i))
    log_msg(sprintf("  Saving chunk to: %s", chunk_path))
    saveRDS(chunk_merged, chunk_path)
    
    chunk_files <- c(chunk_files, chunk_path)
    rm(chunk_merged); gc()
    log_msg(sprintf("  Chunk %d saved, after cleanup [mem: %s]", i, mem_gb()))
    log_msg(sprintf("  ----------------------------------------"))
  }
  
  log_msg(sprintf("=== Final merge of %d chunk files === [mem: %s]",
                  length(chunk_files), mem_gb()))
  final <- merge_file_list(chunk_files)
  log_msg(sprintf("Final object: %d cells x %d genes [mem: %s]",
                  ncol(final), nrow(final), mem_gb()))
  
  # cleanup temp files
  # log_msg("Cleaning up temp files...")
  # unlink(c(tmp_files, chunk_files))
  # log_msg("Temp files deleted")
  
  log_msg(sprintf("=== merge_cancer_cell_types done: %d cells across %d projects === [mem: %s]",
                  ncol(final),
                  length(unique(final$project_name)),
                  mem_gb()))
  final
}

basic_normalizations <- function(merged, npcs = 50) {
  merged <- NormalizeData(merged)               # v5 will normalize each layer
  merged <- FindVariableFeatures(merged, nfeatures = 3000)
  merged <- ScaleData(merged)
  merged <- RunPCA(merged, npcs = npcs)
  return(merged)
  
}
