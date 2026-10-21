#!/usr/bin/env Rscript
# 01b_reanalyze_labels.R
# ---------------------------------------------------------------------------
# Redo the "which MP11 genes fail to separate tumour from adjacent-normal
# proximal tubule" step using the SOURCE-STUDY per-cell malignant calls that
# the Ch4 Kidney epithelial object carries (cell_type / majority_cell_type /
# original_de_dup: "Malignant" vs non-malignant "Epithelial"), instead of the
# HIF-score proxy used in 01_inspect_kidney_epi.R.
#
# Reads the artefacts 01 already wrote (no need to reload the 8 GB object):
#   kidney_epi_compact.rds   -- $expr : 99920 cells x 70 genes (log-norm)
#   kidney_epi_metadata.csv  -- all 53 meta columns, SAME cell order
# ---------------------------------------------------------------------------

PD <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical/purity_adjustment"
GENE <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical/gene_lists/KIRC_MP11.txt"
say <- function(...) cat(sprintf(...), "\n")

mp11 <- trimws(readLines(GENE)); mp11 <- mp11[nzchar(mp11)]

cmp  <- readRDS(file.path(PD, "kidney_epi_compact.rds"))
expr <- cmp$expr
md   <- read.csv(file.path(PD, "kidney_epi_metadata.csv"), stringsAsFactors = FALSE,
                 check.names = FALSE)
stopifnot(nrow(expr) == nrow(md))
say("cells: %d ; MP11 genes present: %d/%d", nrow(expr),
    length(intersect(mp11, colnames(expr))), length(mp11))
g_use <- intersect(mp11, colnames(expr))

## ---- tumour / normal groups from source-study labels -------------------
lab <- md$cell_type
lab[is.na(lab) | lab %in% c("", "NA")] <- NA
# also fold in majority_cell_type / original_de_dup where cell_type is blank
mj  <- md$majority_cell_type;  mj[mj %in% c("", "NA")] <- NA
od  <- md$original_de_dup
od2 <- ifelse(grepl("Malignant", od), "Malignant",
        ifelse(grepl("Epithelial", od), "Epithelial", NA))
lab_comb <- ifelse(!is.na(lab), lab, ifelse(!is.na(mj), mj, od2))

grp <- ifelse(lab_comb == "Malignant", "malignant",
        ifelse(lab_comb == "Epithelial", "normal", NA))
say("\nsource-study epithelial label (cell_type>majority>original_de_dup):")
print(table(grp, useNA = "ifany"))

is_pt <- grepl("proximal tubule|proximal convoluted", md$scimilarity_celltype, ignore.case = TRUE)
say("PT cells among labelled: malignant %d, normal %d",
    sum(is_pt & grp == "malignant", na.rm = TRUE),
    sum(is_pt & grp == "normal",    na.rm = TRUE))

## ---- per-gene tumour vs normal, pooled and per study ------------------
studies <- md$project_name
per_gene <- function(mask) {
  gg <- grp[mask]; E <- expr[mask, g_use, drop = FALSE]
  m_mal <- colMeans(E[gg == "malignant", , drop = FALSE])
  m_nor <- colMeans(E[gg == "normal",    , drop = FALSE])
  p_mal <- colMeans(E[gg == "malignant", , drop = FALSE] > 0)
  p_nor <- colMeans(E[gg == "normal",    , drop = FALSE] > 0)
  data.frame(gene = g_use,
             mean_malignant = m_mal, mean_normal = m_nor,
             pct_malignant = p_mal, pct_normal = p_nor,
             log2FC_MvsN = log2((expm1(m_mal) + 1e-9) / (expm1(m_nor) + 1e-9)),
             diff_MvsN = m_mal - m_nor,
             row.names = NULL)
}

pooled_pt  <- per_gene(is_pt & !is.na(grp))
pooled_all <- per_gene(!is.na(grp))

## per-study log2FC (studies with >=30 cells in BOTH groups, PT cells)
stud_ok <- names(which(sapply(unique(na.omit(studies)), function(s) {
  m <- is_pt & studies == s & !is.na(grp)
  sum(grp[m] == "malignant") >= 30 && sum(grp[m] == "normal") >= 30
})))
say("\nstudies with >=30 PT cells in both malignant & normal: %s",
    paste(stud_ok, collapse = ", "))
fc_by_study <- sapply(stud_ok, function(s) {
  per_gene(is_pt & studies == s & !is.na(grp))$log2FC_MvsN
})
rownames(fc_by_study) <- g_use
mean_fc_study <- rowMeans(fc_by_study)

tab <- pooled_pt
tab$log2FC_MvsN_meanAcrossStudies <- mean_fc_study[tab$gene]
tab <- tab[order(-tab$log2FC_MvsN_meanAcrossStudies), ]
say("\n===== MP11 genes: malignant vs normal PT (source-study labels) =====")
print(cbind(tab[, c("gene", "mean_malignant", "mean_normal", "pct_malignant",
                    "pct_normal", "log2FC_MvsN")],
            fc_study_mean = round(tab$log2FC_MvsN_meanAcrossStudies, 3)),
      digits = 3, row.names = FALSE)
print(round(fc_by_study[tab$gene, , drop = FALSE], 2))

write.csv(tab, file.path(PD, "mp11_tumor_vs_normal_PT_sourcelabels.csv"), row.names = FALSE)

## ---- decide drop list ------------------------------------------------
## KEEP a gene for "malignant-weighted MP11" only if it is genuinely
## tumour-enriched: mean per-study log2FC(malignant vs normal) >= 0.50
## (>= ~1.4x) AND pooled diff > 0. Everything else is dropped as
## "does not distinguish tumour from adjacent-normal PT".
keep <- tab$gene[tab$log2FC_MvsN_meanAcrossStudies >= 0.50 & tab$diff_MvsN > 0]
drop <- setdiff(mp11, keep)
say("\nKEEP (tumour-enriched, malignant-weighted MP11) [%d]: %s",
    length(keep), paste(keep, collapse = ", "))
say("DROP (non-discriminating / higher in normal) [%d]: %s",
    length(drop), paste(drop, collapse = ", "))

## also a laxer variant (threshold 0.25) for sensitivity
keep_lax <- tab$gene[tab$log2FC_MvsN_meanAcrossStudies >= 0.25 & tab$diff_MvsN > 0]
writeLines(setdiff(mp11, keep_lax), file.path(PD, "mp11_drop_genes_lax.txt"))
writeLines(keep_lax,               file.path(PD, "mp11_malignant_weighted_genes_lax.txt"))

writeLines(drop, file.path(PD, "mp11_drop_genes.txt"))
writeLines(keep, file.path(PD, "mp11_malignant_weighted_genes.txt"))
say("\nwrote mp11_drop_genes.txt (%d), mp11_drop_genes_lax.txt (%d)",
    length(drop), length(setdiff(mp11, keep_lax)))
say("[%s] DONE", Sys.time())
