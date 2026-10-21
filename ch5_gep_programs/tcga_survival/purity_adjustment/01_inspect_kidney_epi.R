#!/usr/bin/env Rscript
# 01_inspect_kidney_epi.R
# ---------------------------------------------------------------------------
# Chapter 4 Kidney epithelial object -- inspection for the MP11 purity /
# "malignant-weighted MP11" analysis (thesis Section 5.3.6 robustness check).
#
# Goals:
#   (a) enumerate all metadata columns; find any tumour-vs-normal / tissue
#       annotation carried over from the source studies.
#   (b) derive a per-cell malignant-vs-normal-appearing proximal-tubule (PT)
#       call from canonical ccRCC (VHL/HIF) markers when no explicit label
#       exists.
#   (c) for every MP11 top-50 gene, tabulate mean normalised expression and
#       % expressing in malignant-PT vs normal-appearing-PT, so we can pick
#       the genes that do NOT separate tumour from adjacent normal at the
#       single-cell level (these get dropped from "malignant-weighted MP11").
#
# Container run:
#   ~/rstudio/rstudio.sbatch --script <CTR path to this file>
#   (paths below are /ofir_grp/... = container view of /groups/ofircohen-group)
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

OBJ  <- "/ofir_grp/users/hadasra/temp_objects/clustered_lineages/Epithelial_Kidney_clustered.RDS"
GENE <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical/gene_lists/KIRC_MP11.txt"
OUT  <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical/purity_adjustment"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n")

mp11 <- readLines(GENE)
mp11 <- trimws(mp11[nzchar(mp11)])
say("MP11 genes (%d): %s", length(mp11), paste(mp11, collapse = ", "))

say("[%s] loading %s", Sys.time(), OBJ)
obj <- readRDS(OBJ)
DefaultAssay(obj) <- "RNA"
say("[%s] loaded: %d cells x %d genes; assay class %s",
    Sys.time(), ncol(obj), nrow(obj), class(obj[["RNA"]])[1])

if (inherits(obj[["RNA"]], "Assay5") &&
    length(SeuratObject::Layers(obj[["RNA"]], search = "counts")) > 1) {
  say("[%s] JoinLayers", Sys.time())
  obj <- JoinLayers(obj)
}

md <- obj@meta.data
say("\n===== META.DATA COLUMNS (%d) =====", ncol(md))
for (cn in colnames(md)) {
  v <- md[[cn]]
  if (is.numeric(v)) {
    say("  %-32s numeric   range [%.3g, %.3g]  NA=%d", cn, min(v, na.rm = TRUE),
        max(v, na.rm = TRUE), sum(is.na(v)))
  } else {
    u <- unique(as.character(v))
    show <- head(u, 25)
    say("  %-32s %-8s  %d levels: %s%s", cn, class(v)[1], length(u),
        paste(show, collapse = " | "), if (length(u) > 25) " ..." else "")
  }
}

## ---- write a full per-cell metadata dump (no expression) ----------------
write.csv(data.frame(cell = rownames(md), md, check.names = FALSE),
          file.path(OUT, "kidney_epi_metadata.csv"), row.names = FALSE)
say("\nwrote kidney_epi_metadata.csv")

## ---- candidate tumour/normal columns ----------------------------------
cand <- grep(paste0("tumor|tumour|normal|malign|tissue|sample_?type|condition|",
                    "disease|status|region|site|origin|biopsy|adjacent|nat|",
                    "histolog|patholog|group|class"),
             colnames(md), ignore.case = TRUE, value = TRUE)
say("\n===== candidate tumour/normal/tissue columns: %s",
    if (length(cand)) paste(cand, collapse = ", ") else "(none)")
for (cn in cand) {
  say("--- %s ---", cn)
  print(sort(table(md[[cn]], useNA = "ifany"), decreasing = TRUE))
}

## also cross-tab project_name / study against those, to see which studies
## contribute normal tissue
proj_col <- intersect(c("project_name", "study", "Study", "dataset", "orig.ident",
                        "Project", "cancer_study"), colnames(md))
say("\nproject/study columns present: %s", paste(proj_col, collapse = ", "))
if (length(proj_col)) {
  say("--- cells per study (%s) ---", proj_col[1])
  print(sort(table(md[[proj_col[1]]]), decreasing = TRUE))
}

## ---- normalise & pull genes of interest ------------------------------
obj <- NormalizeData(obj, verbose = FALSE)

ccrcc_up   <- c("CA9", "NDUFA4L2", "VEGFA", "ANGPTL4", "NNMT", "EGLN3",
                "SPP1", "SLC17A3", "ENO2", "PFKP")               # VHL-loss / HIF, up in malignant ccRCC
pt_normal  <- c("SLC5A12", "SLC22A6", "SLC22A8", "SLC34A1", "ALDOB",
                "SLC13A3", "SLC5A2", "SLC7A13", "GPX3", "MT1G")  # canonical normal PT, NOT in MP11
markers    <- unique(c(mp11, ccrcc_up, pt_normal))
present    <- intersect(markers, rownames(obj))
say("\nmarker genes present in object: %d / %d  (missing: %s)",
    length(present), length(markers),
    paste(setdiff(markers, rownames(obj)), collapse = ", "))

expr <- FetchData(obj, vars = present, layer = "data")          # log-normalised
expr <- as.data.frame(expr)

## per-cell scores
sc <- function(df, g) {
  g <- intersect(g, colnames(df)); if (!length(g)) return(rep(NA_real_, nrow(df)))
  rowMeans(scale(df[, g, drop = FALSE]), na.rm = TRUE)
}
ccrcc_score  <- sc(expr, ccrcc_up)
ptnorm_score <- sc(expr, pt_normal)

## restrict to proximal-tubule cells (the compartment MP11 describes)
ct_col <- intersect(c("scimilarity_celltype", "scimilarity_pred", "cell_type",
                      "consensus_cell_type", "Consensus_Cell_Type",
                      "majority_cell_type"), colnames(md))
ct <- if (length(ct_col)) as.character(md[[ct_col[1]]]) else rep(NA, nrow(md))
is_pt <- grepl("proximal tubule|proximal convoluted", ct, ignore.case = TRUE)
say("\ncell-type column used: %s ; PT cells: %d / %d",
    if (length(ct_col)) ct_col[1] else "(none)", sum(is_pt), length(is_pt))
if (sum(is_pt) < 0.2 * nrow(md)) {           # fallback: all epithelial cells
  say("  -> too few PT cells; using ALL cells for the tumour/normal split")
  is_pt <- rep(TRUE, nrow(md))
}

## malignant call: within PT cells, high ccRCC-HIF score & (relative) low normal-PT score
q <- function(x) ecdf(x[is_pt])(x)
mal_call <- ifelse(!is_pt, NA,
             ifelse(q(ccrcc_score) >= 0.66 & q(ptnorm_score) <= 0.5, "malignant_like",
             ifelse(q(ccrcc_score) <= 0.33, "normal_like", "ambiguous")))
say("\nmalignant call (ccRCC-HIF score tertiles within PT cells):")
print(table(mal_call, useNA = "ifany"))

## if an explicit tumour/normal label exists, ALSO compute the split on it
expl <- NULL
for (cn in cand) {
  vals <- tolower(as.character(md[[cn]]))
  if (any(grepl("normal|adjacent|non[- ]?tumor|healthy|benign", vals)) &&
      any(grepl("tumor|tumour|malign|cancer|ccrcc|rcc", vals))) {
    expl <- cn; break
  }
}
say("\nexplicit tumour/normal label column: %s", if (is.null(expl)) "(none found)" else expl)

## ---- per-gene tumour vs normal table (both splits) -------------------
summ <- function(split_vec, label) {
  keep <- !is.na(split_vec)
  g1 <- split_vec == sort(unique(split_vec[keep]))[1]        # first alphabetically
  data <- expr[keep, intersect(mp11, colnames(expr)), drop = FALSE]
  sv <- split_vec[keep]
  grps <- sort(unique(sv))
  res <- lapply(intersect(mp11, colnames(expr)), function(g) {
    m <- tapply(data[[g]], sv, mean)
    p <- tapply(data[[g]] > 0, sv, mean)
    out <- data.frame(gene = g, split = label,
                      grpA = grps[1], grpB = grps[length(grps)],
                      mean_A = as.numeric(m[grps[1]]),
                      mean_B = as.numeric(m[grps[length(grps)]]),
                      pct_A  = as.numeric(p[grps[1]]),
                      pct_B  = as.numeric(p[grps[length(grps)]]))
    out$log2FC_BvsA <- log2((out$mean_B * (exp(1) - 1) + 1)) -   # crude, on log1p scale
                       log2((out$mean_A * (exp(1) - 1) + 1))
    out$abs_mean_diff <- out$mean_B - out$mean_A
    out
  })
  do.call(rbind, res)
}

tab_mal <- summ(mal_call, "ccRCC_HIF_score")
tabs <- tab_mal
if (!is.null(expl)) {
  lab <- tolower(as.character(md[[expl]]))
  lab2 <- ifelse(grepl("normal|adjacent|non[- ]?tumor|healthy|benign", lab), "normal",
          ifelse(grepl("tumor|tumour|malign|cancer|ccrcc|rcc", lab), "tumor", NA))
  lab2[!is_pt] <- NA
  tabs <- rbind(tabs, summ(lab2, paste0("label:", expl)))
}
write.csv(tabs, file.path(OUT, "mp11_tumor_vs_normal_PT.csv"), row.names = FALSE)
say("\nwrote mp11_tumor_vs_normal_PT.csv")
say("\n===== MP11 gene tumour(B) vs normal(A) -- ccRCC-HIF-score split =====")
o <- tab_mal[order(abs(tab_mal$abs_mean_diff)), ]
print(o[, c("gene","mean_A","mean_B","pct_A","pct_B","abs_mean_diff")], digits = 3)

## proposed drop list: genes with small tumour/normal separation
##   |mean diff| < 0.25 log-norm units AND |pct diff| < 0.10
drop <- o$gene[abs(o$abs_mean_diff) < 0.25 & abs(o$pct_A - o$pct_B) < 0.10]
say("\nPROPOSED non-discriminating genes to drop from malignant-weighted MP11 (%d):", length(drop))
say("  %s", paste(drop, collapse = ", "))
writeLines(drop, file.path(OUT, "mp11_drop_genes.txt"))
writeLines(setdiff(mp11, drop), file.path(OUT, "mp11_malignant_weighted_genes.txt"))

## also save compact per-cell frame for any re-analysis without the 8GB object
saveRDS(list(meta = md[, intersect(c(ct_col, proj_col, cand,
                                     grep("louvain", colnames(md), value = TRUE)),
                                   colnames(md)), drop = FALSE],
             expr = expr, ccrcc_score = ccrcc_score, ptnorm_score = ptnorm_score,
             mal_call = mal_call, is_pt = is_pt),
        file.path(OUT, "kidney_epi_compact.rds"))
say("\nwrote kidney_epi_compact.rds")
say("[%s] DONE", Sys.time())
