# 01_subsample.R -- draw the fixed 200,000-cell Figure 1/2 subsample.
#
# Method (from Prompt 15 scoping, job 20767794):
#   * N = 200,000, seed = 42
#   * plain PROPORTIONAL stratified sampling on the joint
#     (cancer_type x Consensus_Cell_Type) cross-tab of the DEDUPED atlas
#     (8,107,665 cells; 5 cross-repo Breastscape duplicate studies removed).
#     -> preserves BOTH marginals simultaneously; no floor / no oversample.
#   * per joint-cell target = round(n_ij * 200000 / N_atlas), largest-remainder
#     rounding to land on exactly 200,000; within a joint cell: uniform random.
#
# Outputs (to temp_objects/figure1_2_scanvi/):
#   drawn_cells.csv            one row / drawn cell: cell_id, project_name,
#                              cancer_type, Consensus_Cell_Type, maj3_cell,
#                              maj2_cell, maj2_cluster, source_repo, agreement
#   per_project_barcodes/<project_name>.txt   barcodes to pull from each step4 RDS
#   composition_check.csv      realized vs expected per cancer_type & per meta-type
#   subsample_params.txt       exact parameters / R version / drift summary
options(width = 200, stringsAsFactors = FALSE)

N_SUB <- 200000L
SEED  <- 42L
TO    <- "/ofir_grp/users/hadasra/temp_objects"
ROOT  <- file.path(TO, "figure1_2_scanvi")
BCDIR <- file.path(ROOT, "per_project_barcodes")
dir.create(BCDIR, recursive = TRUE, showWarnings = FALSE)
lg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))

DUPS <- c("WU_SZ_ET_AL_2021_P15_filtered","GRIFFITHS_JI_ET_AL_2021_P18_filtered",
          "GONZALEZ_H_ET_AL_2022_P31_filtered","SAVAS_P_ET_AL_2018_P35_filtered",
          "GAO_R_ET_AL_2021_P17_filtered")

## ------------------------------------------------------------------ load ------
f <- file.path(TO, "all_cells.gapfix_20260727_164657.RDS")
lg("reading ", f)
md <- as.data.frame(readRDS(f))
lg("raw: ", paste(dim(md), collapse = " x "))
n_raw <- nrow(md)
md <- md[!(md$project_name %in% DUPS), , drop = FALSE]
N_ATLAS <- nrow(md)
lg(sprintf("dropped %d dup rows -> atlas N = %d", n_raw - N_ATLAS, N_ATLAS))
stopifnot(N_ATLAS == 8107665L)
## NB: cell_id is NOT unique across projects in the deduped atlas -> we key on
## the composite (project_name, cell_id) == row_uid, enforced below.
n_dup_cellid <- sum(duplicated(md$cell_id))
lg("non-unique cell_id rows (cross-project collisions): ", n_dup_cellid)
stopifnot(!anyDuplicated(paste(md$project_name, md$cell_id, sep = "|")))

## ---------------------------------------------- derived: source repository ---
## ScBaseCount: SRX_accession non-NA (exact, per Ch.3 audit).
## Breastscape : project_name ends _filtered (all Breastscape carry the suffix
##               in this table) and is not a bare SRP/ERP accession.
## 3CA         : everything else (Author<Year>_<Type>).
repo <- rep(NA_character_, N_ATLAS)
repo[!is.na(md$SRX_accession)] <- "ScBaseCount"
repo[is.na(repo) & grepl("_filtered$", md$project_name)] <- "Breastscape"
repo[is.na(repo) & grepl("^(SRP|ERP|DRP|SRR|ERR)[0-9]+$", md$project_name)] <- "ScBaseCount"
repo[is.na(repo)] <- "3CA"
md$source_repo <- repo
lg("source_repo split (deduped atlas):")
print(table(md$source_repo))
## cross-checks
lg("  cross-check ScBaseCount vs accession regex: ",
   sum(md$source_repo == "ScBaseCount"), " vs ",
   sum(grepl("^(SRP|ERP|DRP|SRR|ERR)[0-9]+$", md$project_name)))
lg("  3CA rows with non-NA cell_type (author label): ",
   sum(md$source_repo == "3CA" & !is.na(md$cell_type)), " / ",
   sum(md$source_repo == "3CA"))

## ---------------------------------------------- derived: agreement level -----
## Prompt 15 (4d): maj3_cell -> "Full"; maj2_cell / maj2_cluster -> "Partial";
##                 Consensus_Cell_Type == "other" -> "Other/Unknown".
is_real <- function(x) !is.na(x) & !(x %in% c("", "NA", "equivocal", "Unknown", "unknown"))
agreement <- ifelse(is_real(md$maj3_cell), "Full",
             ifelse(is_real(md$maj2_cell) | is_real(md$maj2_cluster), "Partial", NA_character_))
agreement[md$Consensus_Cell_Type == "other"] <- "Other/Unknown"
agreement[is.na(agreement)] <- "Other/Unknown"
md$agreement <- agreement
lg("agreement level split (deduped atlas):")
print(table(md$agreement))

## --------------------------------------------- proportional joint-cell draw --
can <- as.character(md$cancer_type)
lin <- as.character(md$Consensus_Cell_Type)
key <- paste(can, lin, sep = "\r")
tab <- table(key)
frac <- N_SUB / N_ATLAS
raw  <- as.numeric(tab) * frac
base <- floor(raw)
rem  <- N_SUB - sum(base)
if (rem > 0) {
  ord <- order(raw - base, decreasing = TRUE)
  headroom <- as.numeric(tab) - base
  ord <- ord[headroom[ord] > 0]
  take <- ord[seq_len(min(rem, length(ord)))]
  base[take] <- base[take] + 1
}
alloc <- setNames(as.integer(base), names(tab))
stopifnot(sum(alloc) == N_SUB)
lg("allocated exactly ", sum(alloc), " across ", sum(alloc > 0), " non-empty joint cells")

set.seed(SEED)
idx_by_key <- split(seq_len(N_ATLAS), key)
sel <- unlist(lapply(names(alloc), function(k) {
  pool <- idx_by_key[[k]]
  a <- alloc[[k]]
  if (is.null(pool) || a <= 0) return(integer(0))
  if (a >= length(pool)) return(pool)
  sample(pool, a)
}), use.names = FALSE)
sel <- sort(unique(sel))
stopifnot(length(sel) == N_SUB)
lg("selected ", length(sel), " unique cells")

drawn <- md[sel, c("cell_id","project_name","cancer_type","Consensus_Cell_Type",
                   "maj3_cell","maj2_cell","maj2_cluster","source_repo","agreement")]
## cell_id is NOT unique across projects in the deduped atlas -> composite key
drawn$row_uid <- paste(drawn$project_name, drawn$cell_id, sep = "|")
stopifnot(!anyDuplicated(drawn$row_uid))
drawn <- drawn[, c("row_uid", setdiff(colnames(drawn), "row_uid"))]
write.csv(drawn, file.path(ROOT, "drawn_cells.csv"), row.names = FALSE)
lg("wrote drawn_cells.csv")

## per-project barcode lists (for the step4 counts pull)
bp <- split(drawn$cell_id, drawn$project_name)
for (p in names(bp)) writeLines(bp[[p]], file.path(BCDIR, paste0(p, ".txt")))
pj <- data.frame(project_name = names(bp), n_drawn = lengths(bp))
pj <- pj[order(-pj$n_drawn), ]
write.csv(pj, file.path(ROOT, "per_project_n_drawn.csv"), row.names = FALSE)
lg("wrote ", nrow(pj), " per-project barcode files; total ", sum(pj$n_drawn))

## ----------------------------------------- composition-preservation check ----
comp <- function(col) {
  full <- table(factor(md[[col]]))
  sub  <- table(factor(drawn[[col]], levels = names(full)))
  data.frame(
    variable       = col,
    level          = names(full),
    n_atlas        = as.integer(full),
    pct_atlas      = round(100 * as.integer(full) / N_ATLAS, 4),
    expected_200k  = round(as.integer(full) * frac, 1),
    n_drawn        = as.integer(sub),
    pct_drawn      = round(100 * as.integer(sub) / N_SUB, 4),
    pct_point_diff = round(100 * as.integer(sub) / N_SUB - 100 * as.integer(full) / N_ATLAS, 4),
    row.names = NULL)
}
cc <- rbind(comp("cancer_type"), comp("Consensus_Cell_Type"),
            comp("source_repo"), comp("agreement"))
write.csv(cc, file.path(ROOT, "composition_check.csv"), row.names = FALSE)
cat("\n=== COMPOSITION CHECK ===\n"); print(cc, row.names = FALSE)
lg("max |pct_point_diff| cancer_type      : ",
   sprintf("%.4f", max(abs(cc$pct_point_diff[cc$variable == "cancer_type"]))))
lg("max |pct_point_diff| Consensus_Cell_Type: ",
   sprintf("%.4f", max(abs(cc$pct_point_diff[cc$variable == "Consensus_Cell_Type"]))))

## joint-cell realized vs expected (full detail, scoping-style)
jc <- as.data.frame(table(cancer_type = can, lineage = lin), responseName = "n_atlas")
jc$expected_200k <- round(jc$n_atlas * frac, 1)
dk <- as.data.frame(table(cancer_type = drawn$cancer_type,
                          lineage = drawn$Consensus_Cell_Type), responseName = "n_drawn")
jc <- merge(jc, dk, all.x = TRUE); jc$n_drawn[is.na(jc$n_drawn)] <- 0
jc$diff <- jc$n_drawn - jc$expected_200k
jc <- jc[order(-jc$n_atlas), ]
write.csv(jc, file.path(ROOT, "composition_check_jointcells.csv"), row.names = FALSE)
lg("max |n_drawn - expected| over joint cells: ", sprintf("%.1f", max(abs(jc$diff))))

## ------------------------------------------------------------- params dump ---
pf <- file.path(ROOT, "subsample_params.txt")
writeLines(c(
  "Figure 1 / Figure 2 -- 200k-cell atlas scANVI subsample",
  sprintf("generated              : %s", format(Sys.time())),
  sprintf("R version              : %s", paste(R.version$major, R.version$minor, sep = ".")),
  sprintf("source metadata        : %s", f),
  sprintf("atlas N (deduped)      : %d  (raw %d - %d rows from 5 cross-repo dup studies)",
          N_ATLAS, n_raw, n_raw - N_ATLAS),
  sprintf("excluded dup studies   : %s", paste(DUPS, collapse = ", ")),
  sprintf("subsample N            : %d", N_SUB),
  sprintf("random seed            : %d  (set.seed, base R Mersenne-Twister)", SEED),
  sprintf("sampling fraction      : %.6f", frac),
  "stratification         : joint (cancer_type x Consensus_Cell_Type) cross-tab",
  "per-stratum target     : round(n_ij * frac), largest-remainder to exactly N_SUB",
  "within stratum         : uniform random, no replacement",
  sprintf("non-empty strata       : %d", sum(alloc > 0)),
  sprintf("max cancer_type drift  : %.4f pct-points",
          max(abs(cc$pct_point_diff[cc$variable == "cancer_type"]))),
  sprintf("max meta-type drift    : %.4f pct-points",
          max(abs(cc$pct_point_diff[cc$variable == "Consensus_Cell_Type"]))),
  sprintf("projects contributing  : %d", nrow(pj))
), pf)
lg("wrote ", pf)
lg("DONE 01_subsample")
