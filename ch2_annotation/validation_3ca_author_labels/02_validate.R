# 02_validate.R -- Validate Chapter 2 consensus meta-type calls against the
# original 3CA author cell-type labels (held-out: not used in the SingleR/
# SCType/Scimilarity consensus vote). Produces confusion matrix, per-meta-type
# precision/recall/F1, overall accuracy + macro-F1, coverage, per-cancer-type
# breakdown, Other/Unknown representation check, and the mapping table used.

options(width = 220, stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(dplyr) })

WT  <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/3ca-author-label-validation-35ef05/work"
TO  <- "/ofir_grp/users/hadasra/temp_objects"
MAP_CSV <- "/ofir_grp/users/hadasra/create_sc_object_workflow/plots_and_outputs/Table2_ontology_mapping.csv"
ONT_XLSX <- "/ofir_grp/users/hadasra/create_sc_object_workflow/Ontology_4_methods.xlsx"
f <- file.path(TO, "all_cells.gapfix_20260727_164657.RDS")

say <- function(...) cat(..., "\n")
hr  <- function(t="") cat("\n", paste0(rep("=",78),collapse=""), "\n", t, "\n", sep="")

hr("LOAD")
say("file:", f, " mtime:", format(file.info(f)$mtime))
md <- as.data.frame(readRDS(f))
say("nrow (gapfix table):", nrow(md), "  ncol:", ncol(md))

# ---- canonical 9 meta-types -------------------------------------------------
canon   <- c("Epithelial","T-cells","Myeloid","B-cell","Fibroblasts",
             "Endothelial","NK-cells","Hematopoietic_Stem_and_Progenitor_cells","other")
disp    <- c("Epithelial/Cancer","T cells","Myeloid","B cells","Fibroblasts/Stroma",
             "Endothelial/Vascular","NK cells","HSC/Stem-like","Other/Unknown")
names(disp) <- canon

# ---- 5 cross-repo Breastscape duplicate studies to drop -------------------
dup5 <- c("WU_SZ_ET_AL_2021_P15_filtered","GAO_R_ET_AL_2021_P17_filtered",
          "GRIFFITHS_JI_ET_AL_2021_P18_filtered","GONZALEZ_H_ET_AL_2022_P31_filtered",
          "SAVAS_P_ET_AL_2018_P35_filtered")
hr("DEDUP ARITHMETIC")
say("rows in the 5 flagged Breastscape dup studies:")
print(sapply(dup5, function(x) sum(md$project_name == x, na.rm=TRUE)))
n_dup5 <- sum(md$project_name %in% dup5)
say("total in 5:", n_dup5)
say("nrow after removing 5:", nrow(md) - n_dup5, " (target corrected atlas: 8,107,665)")
n_atlas_corrected <- nrow(md) - n_dup5

# ---- derive repository of origin -----------------------------------------
hr("REPOSITORY OF ORIGIN")
pn <- md$project_name
is_bsc <- grepl("_ET_AL_.*_filtered$", pn)              # Breastscape: SURNAME_I_ET_AL_YEAR_Pn_filtered
is_scb <- grepl("^[SE]RP[0-9]+$", pn)                   # ScBaseCount: SRP/ERP accession
is_3ca <- !is_bsc & !is_scb                             # 3CA: Author<Year>_<CancerType>
repo <- ifelse(is_bsc, "Breastscape", ifelse(is_scb, "ScBaseCount", "3CA"))
say("cells per repository (project_name pattern):")
print(table(repo))
say("\ncross-check vs repo-specific columns (non-missing counts):")
chk <- data.frame(
  repo = repo,
  has_ORIGIN        = !is.na(md$ORIGIN),
  has_SRX_accession = !is.na(md$SRX_accession),
  has_cell_name     = !is.na(md$cell_name),
  has_cell_type     = !is.na(md$cell_type) & md$cell_type != "" ,
  has_original_dd   = !is.na(md$original_de_dup) & md$original_de_dup != ""
)
print(aggregate(. ~ repo, data = chk, FUN = sum))
say("\nn distinct project_name per repo:")
print(tapply(pn, repo, function(x) length(unique(x))))

# ---- 3CA author label column -------------------------------------------
hr("3CA AUTHOR LABEL COLUMN")
lab_raw  <- md$cell_type                       # raw author-supplied label
lab_dd   <- md$original_de_dup                 # lightly de-duplicated variant
blank <- function(x) is.na(x) | x %in% c("", "NA", "nan", "NaN", "unknown", "Unknown", "Unassigned")
say("cell_type        : non-missing =", sum(!blank(lab_raw)),
    " distinct =", length(unique(lab_raw[!blank(lab_raw)])))
say("original_de_dup  : non-missing =", sum(!blank(lab_dd)),
    " distinct =", length(unique(lab_dd[!blank(lab_dd)])))
say("\ncell_type non-missing by repo:")
print(tapply(!blank(lab_raw), repo, sum))
say("\n=> author label is populated ONLY for 3CA cells:",
    all(tapply(!blank(lab_raw), repo, sum)[c("Breastscape","ScBaseCount")] == 0))

# restrict to 3CA from here on
i3 <- which(repo == "3CA")
say("\n3CA cells total:", length(i3))
say("3CA cells with non-missing cell_type:", sum(!blank(lab_raw[i3])))

hr("3CA AUTHOR LABEL VOCABULARY (raw cell_type, 3CA cells only)")
tab_lab <- sort(table(lab_raw[i3][!blank(lab_raw[i3])]), decreasing = TRUE)
say("n distinct:", length(tab_lab))
print(tab_lab)

# ---- build author-label -> meta-type mapping --------------------------
hr("MAPPING TABLE")
m2 <- read.csv(MAP_CSV, check.names = FALSE)
names(m2) <- c("Method","Original_Label","Assigned_Meta_type")
m2$key <- tolower(trimws(m2$Original_Label))
m2$Assigned_Meta_type <- trimws(m2$Assigned_Meta_type)
orig_rows <- m2[m2$Method == "Original", ]
say("Table2_ontology_mapping.csv: ", nrow(m2), " rows total; Method=='Original' rows: ", nrow(orig_rows))

# also pull 'Original' rows from the ontology xlsx if present
xls_orig <- tryCatch({
  x <- openxlsx::read.xlsx(ONT_XLSX)
  names(x) <- trimws(names(x))
  if (all(c("Table","Cell_type","Meta") %in% names(x))) {
    x <- x[x$Table == "Original", c("Cell_type","Meta")]
    data.frame(key = tolower(trimws(x$Cell_type)), meta = trimws(x$Meta))
  } else data.frame(key=character(0), meta=character(0))
}, error = function(e) { say("xlsx read failed:", conditionMessage(e)); data.frame(key=character(0), meta=character(0)) })
say("ontology xlsx 'Original' rows: ", nrow(xls_orig))

# curated extension for residual 3CA raw labels NOT already in Table2 'Original'.
# Rationale follows the SAME rules already visible in the Table2 'Original' rows
# (singular/plural & underscore variants of an already-mapped label -> same meta-type).
extra_map <- c(
  "malignant"="Epithelial", "epithelial"="Epithelial", "epithelial cells"="Epithelial",
  "keratinocytes"="Epithelial", "hepatocyte"="other", "emt"="Epithelial",
  "acinar"="other", "ductal_1"="other", "ductal_2"="other", "alveolar"="other",
  "t_cell"="T-cells", "cd4+ t-cells"="T-cells", "cd8+ t-cells"="T-cells",
  "tgd"="T-cells", "til"="T-cells", "lymphoid"="T-cells",
  "t_nk_cell"="NK-cells", "nk_cell"="NK-cells", "ilc"="NK-cells",
  "b_cell"="B-cell", "plasma"="B-cell", "plasmablast"="B-cell",
  "myeloid"="Myeloid", "macrophage"="Myeloid", "monocyte"="Myeloid",
  "dendritic"="Myeloid", "dc"="Myeloid", "pdc"="Myeloid", "langerhans"="Myeloid",
  "mast"="Myeloid", "granulocyte"="Myeloid", "neutrophil"="Myeloid",
  "fibroblast"="Fibroblasts", "fibroblasts"="Fibroblasts", "myofibroblast"="Fibroblasts",
  "msc"="Fibroblasts", "mesenchymal"="Fibroblasts", "adipocytes"="Fibroblasts",
  "stellate"="Fibroblasts", "schwann"="Fibroblasts", "mesangial cells"="Fibroblasts",
  "endothelial"="Endothelial", "pericyte"="Endothelial",
  "vascular_smooth_muscle"="Endothelial",
  "hsc"="Hematopoietic_Stem_and_Progenitor_cells",
  "hpc-like"="Hematopoietic_Stem_and_Progenitor_cells",
  "erythroblast"="Hematopoietic_Stem_and_Progenitor_cells",
  "erythrocytes"="Hematopoietic_Stem_and_Progenitor_cells",
  "megakaryocyte"="Hematopoietic_Stem_and_Progenitor_cells",
  "chondrocytes"="Hematopoietic_Stem_and_Progenitor_cells",
  "osteoblast"="other","osteoclast"="other","osteoblast_proli"="other",
  "chondrocyte"="other","myocyte"="other","myocytes"="other","myoblast"="other",
  "smooth_muscle"="other","skeletal muscle"="other","endocrine"="other",
  "neuron"="other","astrocyte"="other","oligodendrocyte"="other","glia"="other",
  "mesothelial"="other","melanocytes"="other","ovary_stromal"="other",
  "lymphovascular"="other","etc"="other","tzbtb16"="other",
  "multiplet"="other","doublet"="other"
)

lut <- c(setNames(orig_rows$Assigned_Meta_type, orig_rows$key),
         setNames(xls_orig$meta, xls_orig$key))
lut <- lut[!duplicated(names(lut))]
# add extras only where not already present from the thesis tables
for (k in names(extra_map)) if (is.na(lut[k]) || is.null(lut[k])) lut[k] <- extra_map[k]

map_label <- function(v) {
  out <- unname(lut[tolower(trimws(v))])
  out[!out %in% canon] <- NA
  out
}
mapped3 <- map_label(lab_raw[i3])
mapped3[blank(lab_raw[i3])] <- NA

# report the mapping actually used, with provenance + cell counts
lab_names <- names(tab_lab)
prov <- ifelse(tolower(trimws(lab_names)) %in% orig_rows$key, "Table2 'Original' row",
        ifelse(tolower(trimws(lab_names)) %in% xls_orig$key, "ontology xlsx 'Original' row",
        ifelse(tolower(trimws(lab_names)) %in% names(extra_map), "added (this analysis)", "UNMAPPED")))
map_used <- data.frame(
  author_label   = lab_names,
  n_cells_3CA    = as.integer(tab_lab),
  assigned_meta  = unname(lut[tolower(trimws(lab_names))]),
  assigned_disp  = disp[unname(lut[tolower(trimws(lab_names))])],
  provenance     = prov
)
map_used$assigned_meta[!map_used$assigned_meta %in% canon] <- NA
map_used$provenance[is.na(map_used$assigned_meta)] <- "UNMAPPED"
say("\n--- author-label -> meta-type mapping used ---")
print(map_used, row.names = FALSE)
write.csv(map_used, file.path(WT, "author_label_mapping_used.csv"), row.names = FALSE)

hr("UNMAPPED AUTHOR LABELS")
um <- map_used[is.na(map_used$assigned_meta), c("author_label","n_cells_3CA")]
um <- um[order(-um$n_cells_3CA), ]
print(um, row.names = FALSE)
say("total unmapped 3CA cells (have a label but no clean meta-type):", sum(um$n_cells_3CA))
write.csv(um, file.path(WT, "unmapped_labels.csv"), row.names = FALSE)

# ---- confusion matrix + metrics --------------------------------------
hr("CONFUSION MATRIX  (rows = consensus meta_type / Consensus_Cell_Type ; cols = mapped 3CA author label = ground truth)")
cons3 <- md$Consensus_Cell_Type[i3]
keep  <- !is.na(mapped3) & !is.na(cons3)
say("3CA cells in confusion matrix (both consensus + mapped author label):", sum(keep))
pred <- factor(cons3[keep],    levels = canon)
truth<- factor(mapped3[keep],  levels = canon)
cm <- table(consensus = pred, author_3CA = truth)
cm_disp <- cm; dimnames(cm_disp) <- list(consensus = disp[canon], author_3CA = disp[canon])
print(cm_disp)
write.csv(as.data.frame.matrix(cm_disp), file.path(WT, "confusion_matrix_9x9.csv"))

N <- sum(cm)
acc <- sum(diag(cm)) / N
prec <- diag(cm) / rowSums(cm)
rec  <- diag(cm) / colSums(cm)
f1   <- 2*prec*rec / (prec+rec)
metrics <- data.frame(
  meta_type = disp[canon],
  n_consensus = as.integer(rowSums(cm)),
  n_author_3CA = as.integer(colSums(cm)),
  precision = round(prec,4), recall = round(rec,4), f1 = round(f1,4)
)
say("\n--- per-meta-type precision / recall / F1 (3CA author label = ground truth) ---")
print(metrics, row.names = FALSE)
macroF1 <- mean(f1, na.rm = TRUE)
say(sprintf("\nOverall accuracy (9-class): %.4f  (%d / %d)", acc, sum(diag(cm)), N))
say(sprintf("Macro-averaged F1 (9-class): %.4f", macroF1))

# 8-class variant (drop Other/Unknown from both sides)
c8 <- setdiff(canon, "other")
k8 <- keep & mapped3 %in% c8 & cons3 %in% c8
p8 <- factor(cons3[k8], levels=c8); t8 <- factor(mapped3[k8], levels=c8)
cm8 <- table(p8,t8)
acc8 <- sum(diag(cm8))/sum(cm8)
pr8 <- diag(cm8)/rowSums(cm8); rc8 <- diag(cm8)/colSums(cm8); f18 <- 2*pr8*rc8/(pr8+rc8)
say(sprintf("\n[sensitivity] 8-class (excl. Other/Unknown both sides): n=%d  accuracy=%.4f  macroF1=%.4f",
            sum(cm8), acc8, mean(f18, na.rm=TRUE)))
write.csv(metrics, file.path(WT, "per_metatype_metrics.csv"), row.names = FALSE)

# ---- coverage --------------------------------------------------------
hr("COVERAGE")
n_3ca        <- length(i3)
n_3ca_lab    <- sum(!blank(lab_raw[i3]))
n_3ca_mapped <- sum(keep)
cov <- data.frame(
  quantity = c("3CA-origin cells (total)",
               "3CA cells with non-missing author label",
               "3CA cells whose author label mapped to one of 9 meta-types",
               "validated subset as % of corrected atlas (8,107,665)",
               "validated subset as % of gapfix table (8,270,370)",
               "validated subset as % of thesis atlas (8,113,115)"),
  value = c(n_3ca, n_3ca_lab, n_3ca_mapped,
            sprintf("%.2f%%", 100*n_3ca_mapped/n_atlas_corrected),
            sprintf("%.2f%%", 100*n_3ca_mapped/nrow(md)),
            sprintf("%.2f%%", 100*n_3ca_mapped/8113115)))
print(cov, row.names = FALSE)
write.csv(cov, file.path(WT, "coverage.csv"), row.names = FALSE)

# ---- per-cancer-type breakdown -------------------------------------
hr("PER-CANCER-TYPE BREAKDOWN (validated 3CA subset)")
ct <- md$cancer_type[i3][keep]
per_ct <- do.call(rbind, lapply(sort(unique(ct)), function(cc) {
  ii <- ct == cc
  cmi <- table(factor(cons3[keep][ii], levels=canon), factor(mapped3[keep][ii], levels=canon))
  ni <- sum(cmi)
  pr <- diag(cmi)/rowSums(cmi); rc <- diag(cmi)/colSums(cmi); f1i <- 2*pr*rc/(pr+rc)
  data.frame(cancer_type = cc, n_cells = ni,
             accuracy = round(sum(diag(cmi))/ni, 4),
             macroF1  = round(mean(f1i, na.rm=TRUE), 4),
             n_metatypes_present = sum(colSums(cmi) > 0))
}))
per_ct <- per_ct[order(per_ct$accuracy), ]
print(per_ct, row.names = FALSE)
write.csv(per_ct, file.path(WT, "per_cancertype_metrics.csv"), row.names = FALSE)

# ---- Other/Unknown representation --------------------------------
hr("OTHER / UNKNOWN REPRESENTATION")
cons_all_corrected <- md$Consensus_Cell_Type[!(md$project_name %in% dup5)]
p_other_atlas <- mean(cons_all_corrected == "other")
p_other_3ca_all <- mean(md$Consensus_Cell_Type[i3] == "other")
p_other_3ca_val <- mean(cons3[keep] == "other")
p_other_authormap <- mean(mapped3[keep] == "other")
orep <- data.frame(
  subset = c("full corrected atlas (8,107,665)",
             "all 3CA cells (by consensus)",
             "validated 3CA subset (by consensus)",
             "validated 3CA subset (by mapped author label)"),
  pct_Other_Unknown = round(100*c(p_other_atlas, p_other_3ca_all, p_other_3ca_val, p_other_authormap), 2))
print(orep, row.names = FALSE)
write.csv(orep, file.path(WT, "other_representation.csv"), row.names = FALSE)
say(sprintf("\nOther/Unknown is %s in the 3CA subset vs the full atlas (%.1f%% vs %.1f%% by consensus).",
            ifelse(p_other_3ca_all > p_other_atlas, "OVER-represented", "UNDER-represented"),
            100*p_other_3ca_all, 100*p_other_atlas))

hr("DONE")
say("outputs written to:", WT)
