#!/usr/bin/env Rscript
# 02_purity_survival.R
# ---------------------------------------------------------------------------
# KIRC MP11 overall-survival result (thesis Section 5.3.6) -- tumour-purity
# robustness analysis.
#
#   1. Reproduce the original unadjusted continuous Cox (z-scored MP11).
#   2. Multivariable Cox: MP11 + age + sex + AJCC stage + grade + CPE purity.
#      Schoenfeld / cox.zph PH check.
#   3. Re-fit the adjusted model on CPE >= 0.6 and CPE >= 0.7 subsets.
#   4. "Malignant-weighted MP11": drop MP11 genes that do not separate
#      tumour from adjacent-normal proximal tubule at single-cell level
#      (list from 01_inspect_kidney_epi.R) and re-run unadjusted + adjusted.
#   5. Multiplicity count for the Section 5.3.6 statement.
#   6. Exact cohort-subsetting numbers (KIRC vs BRCA clarification).
#
# Base R + survival only. Container paths (/ofir_grp = /groups/ofircohen-group).
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(survival))

DD  <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical"
PD  <- file.path(DD, "purity_adjustment")
dir.create(PD, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")
hr_line <- function(m, term = 1) {
  s <- summary(m)
  ci <- s$conf.int[term, ]
  cf <- s$coefficients[term, ]
  sprintf("HR=%.3f  95%% CI %.3f-%.3f  p=%.3g   (n=%d, events=%d)",
          ci["exp(coef)"], ci["lower .95"], ci["upper .95"],
          cf[ncol(s$coefficients)], s$n, s$nevent)
}

## ========================================================================
## 1. CLINICAL  (replicate run_tcga_survival.R filtering exactly)
## ========================================================================
clin <- read.csv(file.path(DD, "TCGA-CDR.csv"), stringsAsFactors = FALSE,
                 na.strings = c("#N/A", "[Not Available]", "[Not Applicable]", ""))
clin$OS      <- suppressWarnings(as.numeric(clin$OS))
clin$OS.time <- suppressWarnings(as.numeric(clin$OS.time))
clin_f <- clin[(is.na(clin$Redaction) | clin$Redaction == "") &
               !is.na(clin$OS) & !is.na(clin$OS.time), ]

say("=== TCGA-CDR cohort accounting ===")
for (ty in c("KIRC", "BRCA")) {
  a <- sum(clin$type == ty, na.rm = TRUE)
  b <- sum(clin_f$type == ty)
  say("  %s : %d rows in TCGA-CDR ; %d after (Redaction blank & OS & OS.time)", ty, a, b)
}

kirc <- clin_f[clin_f$type == "KIRC", ]

## covariates
kirc$age <- suppressWarnings(as.numeric(kirc$age_at_initial_pathologic_diagnosis))
kirc$sex <- factor(kirc$gender)

st <- gsub("Stage\\s+", "", kirc$ajcc_pathologic_tumor_stage)
st <- gsub("[ABC]$", "", trimws(st))
st[grepl("Discrepancy", st, ignore.case = TRUE)] <- NA
st[st %in% c("", "NA")] <- NA
kirc$stage  <- factor(st, levels = c("I", "II", "III", "IV"))
kirc$stage34 <- factor(ifelse(is.na(st), NA, ifelse(st %in% c("III", "IV"), "III-IV", "I-II")),
                       levels = c("I-II", "III-IV"))

gr <- trimws(kirc$histological_grade)
gr[gr %in% c("GX", "", "NA")] <- NA
kirc$grade   <- factor(gr, levels = c("G1", "G2", "G3", "G4"))
kirc$grade3  <- factor(ifelse(is.na(gr), NA,
                       ifelse(gr %in% c("G1", "G2"), "G1-2", gr)),
                       levels = c("G1-2", "G3", "G4"))

say("\nKIRC covariate availability (of %d):", nrow(kirc))
for (v in c("age", "sex", "stage", "grade", "grade3"))
  say("  %-7s missing %d", v, sum(is.na(kirc[[v]])))

## ========================================================================
## 2. EXPRESSION MATRIX  (replicate build_matrix from run_tcga_survival.R)
## ========================================================================
build_matrix <- function(expr_long) {
  patients <- unique(expr_long$patient_barcode)
  genes    <- unique(expr_long$gene)
  mat <- matrix(NA_real_, nrow = length(patients), ncol = length(genes),
                dimnames = list(patients, genes))
  mat[cbind(match(expr_long$patient_barcode, patients),
            match(expr_long$gene, genes))] <- log2(expr_long$tpm + 1)
  scale(mat)                      # per-gene z across patients
}
score_program <- function(matz, genes) {
  present <- intersect(genes, colnames(matz))
  rowMeans(matz[, present, drop = FALSE], na.rm = TRUE)
}

expr_long <- read.csv(file.path(DD, "KIRC_expr_long.csv"), stringsAsFactors = FALSE)
matz <- build_matrix(expr_long)
say("\nKIRC expression matrix: %d patients x %d genes", nrow(matz), ncol(matz))

mp11 <- trimws(readLines(file.path(DD, "gene_lists", "KIRC_MP11.txt")))
mp11 <- mp11[nzchar(mp11)]
drop_file <- file.path(PD, "mp11_drop_genes.txt")
drop_genes <- if (file.exists(drop_file)) trimws(readLines(drop_file)) else character(0)
drop_genes <- drop_genes[nzchar(drop_genes)]
mp11_mw <- setdiff(mp11, drop_genes)
say("MP11 genes: %d ; dropped for malignant-weighted: %d (%s) ; kept: %d",
    length(mp11), length(drop_genes),
    if (length(drop_genes)) paste(drop_genes, collapse = ",") else "-", length(mp11_mw))

score_orig <- score_program(matz, mp11)
score_mw   <- score_program(matz, mp11_mw)

## ========================================================================
## 3. TUMOUR PURITY  (Aran et al. 2015 consensus, via TCGAbiolinks Tumor.purity)
## ========================================================================
nm <- load(file.path(PD, "Tumor.purity.rda"))
tp <- get(nm[1])
say("\nTumor.purity: %d rows, cols: %s", nrow(tp), paste(colnames(tp), collapse = ", "))
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
for (cc in c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"))
  if (cc %in% colnames(tp)) tp[[cc]] <- num(tp[[cc]])

tk <- tp[tp$Cancer.type == "KIRC", ]
tk$patient <- substr(tk$Sample.ID, 1, 12)
tk$stype   <- substr(tk$Sample.ID, 14, 15)
say("KIRC purity rows: %d ; sample-type codes: %s", nrow(tk),
    paste(names(table(tk$stype)), table(tk$stype), sep = ":", collapse = " "))
tk <- tk[tk$stype == "01", ]                      # primary solid tumour only
cpe <- aggregate(cbind(CPE, ESTIMATE, ABSOLUTE) ~ patient, data = tk,
                 FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
say("KIRC primary-tumour patients with a purity row: %d ; with non-NA CPE: %d",
    nrow(cpe), sum(!is.na(cpe$CPE)))

## ========================================================================
## 4. ASSEMBLE ANALYSIS FRAME
## ========================================================================
df <- data.frame(bcr_patient_barcode = rownames(matz),
                 score_orig = score_orig, score_mw = score_mw,
                 stringsAsFactors = FALSE)
df <- merge(df, kirc[, c("bcr_patient_barcode", "OS", "OS.time",
                         "age", "sex", "stage", "stage34", "grade", "grade3")],
            by = "bcr_patient_barcode")
df <- merge(df, cpe, by.x = "bcr_patient_barcode", by.y = "patient", all.x = TRUE)
df$z_mp11    <- as.numeric(scale(df$score_orig))
df$z_mp11_mw <- as.numeric(scale(df$score_mw))

say("\n=== analysis frame ===")
say("  patients with expression + KIRC clinical : %d  (events=%d)", nrow(df), sum(df$OS))
say("  ... also with non-NA CPE                  : %d  (events=%d)",
    sum(!is.na(df$CPE)), sum(df$OS[!is.na(df$CPE)]))
say("  Spearman r(MP11 score, CPE)   = %.3f", cor(df$score_orig, df$CPE, method = "spearman", use = "complete.obs"))
say("  Spearman r(MP11 score, ESTIMATE purity) = %.3f", cor(df$score_orig, df$ESTIMATE, method = "spearman", use = "complete.obs"))
say("  Pearson  r(MP11 orig, MP11 malignant-weighted) = %.3f", cor(df$score_orig, df$score_mw, use = "complete.obs"))

fmt_grp <- function(x) sprintf("%.1f%%", 100 * mean(x))

## ========================================================================
## 5. MODELS
## ========================================================================
res <- list()
add <- function(tag, m, term = "z_mp11") {
  s <- summary(m); ci <- s$conf.int[term, ]; cf <- s$coefficients[term, ]
  res[[length(res) + 1]] <<- data.frame(
    model = tag, term = term, n = s$n, events = s$nevent,
    HR = ci["exp(coef)"], lo95 = ci["lower .95"], hi95 = ci["upper .95"],
    p = cf[ncol(s$coefficients)], row.names = NULL)
  say("[%-42s] %s", tag, hr_line(m, term))
}

say("\n================  UNADJUSTED  ================")
m_raw <- coxph(Surv(OS.time, OS) ~ score_orig, data = df)          # exact original parametrisation
say("original parametrisation (raw mean-z score, matches thesis pipeline):")
say("   %s", hr_line(m_raw, "score_orig"))
m0 <- coxph(Surv(OS.time, OS) ~ z_mp11, data = df); add("unadj_z_MP11_full", m0)
m0c <- coxph(Surv(OS.time, OS) ~ z_mp11, data = df[!is.na(df$CPE), ]); add("unadj_z_MP11_CPEcohort", m0c)

say("\n================  MULTIVARIABLE ADJUSTED (full)  ================")
m_adj <- coxph(Surv(OS.time, OS) ~ z_mp11 + age + sex + stage + grade + CPE, data = df)
print(summary(m_adj))
add("adj_MP11+age+sex+stage+grade+CPE", m_adj)
say("\n-- Schoenfeld / cox.zph --")
zph <- cox.zph(m_adj); print(zph)
saveRDS(zph, file.path(PD, "coxzph_adjusted_full.rds"))

say("\n-- sensitivity: collapsed stage (I-II/III-IV) + grade (G1-2/G3/G4) --")
m_adj_c <- coxph(Surv(OS.time, OS) ~ z_mp11 + age + sex + stage34 + grade3 + CPE, data = df)
add("adj_collapsed_stage_grade", m_adj_c)
print(cox.zph(m_adj_c))

say("\n-- covariate-only vs +MP11 likelihood-ratio test (complete cases) --")
cc <- df[stats::complete.cases(df[, c("OS.time","OS","z_mp11","age","sex","stage","grade","CPE")]), ]
m_cov  <- coxph(Surv(OS.time, OS) ~ age + sex + stage + grade + CPE, data = cc)
m_full <- coxph(Surv(OS.time, OS) ~ z_mp11 + age + sex + stage + grade + CPE, data = cc)
print(anova(m_cov, m_full))

say("\n-- purity-only Cox (continuous CPE) --")
say("   %s", hr_line(coxph(Surv(OS.time, OS) ~ CPE, data = df), "CPE"))

say("\n================  HIGH-PURITY SUBSETS  ================")
for (thr in c(0.6, 0.7)) {
  sub <- df[!is.na(df$CPE) & df$CPE >= thr, ]
  say("--- CPE >= %.1f : n=%d, events=%d ---", thr, nrow(sub), sum(sub$OS))
  ms_u <- coxph(Surv(OS.time, OS) ~ z_mp11, data = sub)
  add(sprintf("subset_CPE>=%.1f_unadj", thr), ms_u)
  ms_a <- tryCatch(coxph(Surv(OS.time, OS) ~ z_mp11 + age + sex + stage + grade + CPE, data = sub),
                   error = function(e) {say("   full adj failed (%s); using collapsed", conditionMessage(e)); NULL})
  if (!is.null(ms_a)) add(sprintf("subset_CPE>=%.1f_adj", thr), ms_a)
  ms_ac <- coxph(Surv(OS.time, OS) ~ z_mp11 + age + sex + stage34 + grade3, data = sub)
  add(sprintf("subset_CPE>=%.1f_adj_collapsed_noCPE", thr), ms_ac)
}

say("\n================  MALIGNANT-WEIGHTED MP11  ================")
if (length(drop_genes) == 0)
  say("(no drop list yet -- run 01_inspect_kidney_epi.R first; z_mp11_mw == z_mp11 here)")
mw0 <- coxph(Surv(OS.time, OS) ~ z_mp11_mw, data = df); add("MW_unadj", mw0, "z_mp11_mw")
mwa <- coxph(Surv(OS.time, OS) ~ z_mp11_mw + age + sex + stage + grade + CPE, data = df)
add("MW_adj_full", mwa, "z_mp11_mw")
print(cox.zph(mwa))
mwac <- coxph(Surv(OS.time, OS) ~ z_mp11_mw + age + sex + stage34 + grade3 + CPE, data = df)
add("MW_adj_collapsed", mwac, "z_mp11_mw")

## ========================================================================
## 6. MULTIPLICITY COUNT  (from the original results table)
## ========================================================================
say("\n================  MULTIPLICITY  ================")
orig <- read.csv(file.path(DD, "tcga_survival_results.csv"), stringsAsFactors = FALSE)
say("original results table rows (program x cohort combinations tested): %d", nrow(orig))
print(orig[, c("cohort", "program", "n", "n_events")])
framings <- c("median-split log-rank", "tertile-split log-rank", "continuous Cox")
say("statistical framings per program: %d  (%s)", length(framings), paste(framings, collapse = "; "))
say("TOTAL nominal tests behind Section 5.3.6 = %d program-cohort x %d framings = %d",
    nrow(orig), length(framings), nrow(orig) * length(framings))
say("MP11 continuous-Cox p = %.3g ; Bonferroni x%d = %.3g",
    orig$p_cox[orig$program == "MP11"], nrow(orig) * length(framings),
    min(1, orig$p_cox[orig$program == "MP11"] * nrow(orig) * length(framings)))

## ========================================================================
out <- do.call(rbind, res)
write.csv(out, file.path(PD, "mp11_purity_survival_results.csv"), row.names = FALSE)
say("\nwrote mp11_purity_survival_results.csv")
saveRDS(df, file.path(PD, "kirc_mp11_analysis_frame.rds"))
say("[%s] DONE", Sys.time())
