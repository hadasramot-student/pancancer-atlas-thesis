# 01_exclusions.R -- For the validated 3CA population (consensus meta-type call
# present AND mapped author label present, same 9-meta-type ontology as the
# prior 3ca-author-label-validation analysis), report cell count + accuracy
# under 4 exclusion definitions:
#   (1) no further exclusion
#   (2) excluding cells where CONSENSUS call is Other/Unknown
#   (3) excluding cells where MAPPED AUTHOR label is Other/Unknown
#   (4) excluding cells where EITHER is Other/Unknown
#
# Reuses the exact population-construction logic from
# 3ca-author-label-validation-35ef05/work/02_validate.R (same RDS, same
# Table2_ontology_mapping.csv + curated extension) so the base population
# should reproduce ~3,262,688 cells.

options(width = 220, stringsAsFactors = FALSE)

WT  <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/atlas-metadata-consensus-validation-e53a85/work"
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

canon <- c("Epithelial","T-cells","Myeloid","B-cell","Fibroblasts",
           "Endothelial","NK-cells","Hematopoietic_Stem_and_Progenitor_cells","other")
disp  <- c("Epithelial/Cancer","T cells","Myeloid","B cells","Fibroblasts/Stroma",
           "Endothelial/Vascular","NK cells","HSC/Stem-like","Other/Unknown")
names(disp) <- canon

dup5 <- c("WU_SZ_ET_AL_2021_P15_filtered","GAO_R_ET_AL_2021_P17_filtered",
          "GRIFFITHS_JI_ET_AL_2021_P18_filtered","GONZALEZ_H_ET_AL_2022_P31_filtered",
          "SAVAS_P_ET_AL_2018_P35_filtered")
n_dup5 <- sum(md$project_name %in% dup5)
n_atlas_corrected <- nrow(md) - n_dup5
say("corrected atlas cell count:", n_atlas_corrected, " (target 8,107,665)")

pn <- md$project_name
is_bsc <- grepl("_ET_AL_.*_filtered$", pn)
is_scb <- grepl("^[SE]RP[0-9]+$", pn)
is_3ca <- !is_bsc & !is_scb
repo <- ifelse(is_bsc, "Breastscape", ifelse(is_scb, "ScBaseCount", "3CA"))

lab_raw <- md$cell_type
blank <- function(x) is.na(x) | x %in% c("", "NA", "nan", "NaN", "unknown", "Unknown", "Unassigned")

i3 <- which(repo == "3CA" & !(md$project_name %in% dup5))
say("3CA cells (post-QC, dedup'd) total:", length(i3))

hr("MAPPING TABLE")
m2 <- read.csv(MAP_CSV, check.names = FALSE)
names(m2) <- c("Method","Original_Label","Assigned_Meta_type")
m2$key <- tolower(trimws(m2$Original_Label))
m2$Assigned_Meta_type <- trimws(m2$Assigned_Meta_type)
orig_rows <- m2[m2$Method == "Original", ]

xls_orig <- tryCatch({
  x <- openxlsx::read.xlsx(ONT_XLSX)
  names(x) <- trimws(names(x))
  if (all(c("Table","Cell_type","Meta") %in% names(x))) {
    x <- x[x$Table == "Original", c("Cell_type","Meta")]
    data.frame(key = tolower(trimws(x$Cell_type)), meta = trimws(x$Meta))
  } else data.frame(key=character(0), meta=character(0))
}, error = function(e) { say("xlsx read failed:", conditionMessage(e)); data.frame(key=character(0), meta=character(0)) })

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
for (k in names(extra_map)) if (is.na(lut[k]) || is.null(lut[k])) lut[k] <- extra_map[k]

map_label <- function(v) {
  out <- unname(lut[tolower(trimws(v))])
  out[!out %in% canon] <- NA
  out
}
mapped3 <- map_label(lab_raw[i3])
mapped3[blank(lab_raw[i3])] <- NA

cons3 <- md$Consensus_Cell_Type[i3]

hr("BASE POPULATION")
keep <- !is.na(mapped3) & !is.na(cons3)
n_base <- sum(keep)
say("Base population (consensus present AND mapped author label present):", n_base)
say("Reference check (prior analysis): 3,262,688  -- diff:", n_base - 3262688)

pred  <- cons3[keep]
truth <- mapped3[keep]

is_other <- function(x) x == "other"

defs <- list(
  "1_no_exclusion"        = rep(TRUE, n_base),
  "2_excl_consensus_other"= !is_other(pred),
  "3_excl_author_other"   = !is_other(truth),
  "4_excl_either_other"   = !is_other(pred) & !is_other(truth)
)

hr("RESULTS")
res <- do.call(rbind, lapply(names(defs), function(nm) {
  sel <- defs[[nm]]
  n <- sum(sel)
  n_match <- sum(pred[sel] == truth[sel])
  data.frame(definition = nm, n_cells = n, n_match = n_match,
             pct_match = round(100 * n_match / n, 4))
}))
print(res, row.names = FALSE)
write.csv(res, file.path(WT, "exclusion_definitions_results.csv"), row.names = FALSE)

hr("DONE")
say("outputs written to:", WT)
