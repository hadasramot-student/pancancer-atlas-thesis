#!/usr/bin/env Rscript
# run_tcga_survival.R
#
# Ch.5 (GEP programs) — outcome linkage. Scores TCGA-BRCA and TCGA-KIRC
# patients on selected high-confidence cross-method (cNMF/GeneNMF) programs
# from the atlas, via mean z-scored log2(TPM+1) expression of each program's
# top-50 genes, then tests for a survival association (median-split log-rank
# + continuous Cox) against TCGA-CDR overall survival (OS, OS.time).
#
# Inputs (under DATA_DIR = ATLAS_DATA_DIR/tcga_clinical):
#   TCGA-CDR.csv               - Liu et al. 2018 pan-cancer clinical resource
#   <Cohort>_expr_long.csv     - long-format patient_barcode,gene,tpm
#   gene_lists/<Cohort>_<Program>.txt - one gene symbol per line, top-50
#
# Programs tested (chosen for high cross-method jaccard + clean, non-MT-dominated
# biology; see Breast/Kidney _highconfidence_programs.csv):
#   BRCA: MP7  (immediate-early/stress, FOS/EGR1/ATF3; matches cNMF TNFA-NFkB)
#         MP10 (hormone-receptor/luminal identity: ESR1/AR/PRLR/MLPH)
#   KIRC: MP11 (proximal-tubule differentiation/metabolism: PCK1/GATM/AZGP1)
#         MP13 (immediate-early/stress: HSPA1A/FOSB/ATF3/EGR1)
#         MP15 (epithelial/luminal differentiation: KRT19/TACSTD2/CLDN4/MUC1)
#
# Output: tcga_survival_results.csv (all programs, both split methods, Cox HR)
#         km_plots/<Cohort>_<Program>_km.png (only for p < 0.1 on either split)
#
# Container note: run inside the apptainer Rstudio image via
# ~/rstudio/rstudio.sbatch --script <this file, CTR path>
# (submitted by submit_tcga_survival.sh, which handles HOST/CTR path translation).

userlib <- Sys.getenv("R_USER_LIBRARY", unset = "")
if (nzchar(userlib)) .libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(survival)
})

DATA_DIR <- file.path(Sys.getenv("ATLAS_DATA_DIR", unset = "data/temp_objects"), "tcga_clinical")
GENE_DIR <- file.path(DATA_DIR, "gene_lists")
OUT_DIR  <- DATA_DIR
KM_DIR   <- file.path(OUT_DIR, "km_plots")
dir.create(KM_DIR, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) cat(sprintf(...), "\n")

# ── Clinical data ──────────────────────────────────────────────────────────
clin <- read.csv(file.path(DATA_DIR, "TCGA-CDR.csv"), stringsAsFactors = FALSE,
                  na.strings = c("#N/A", "[Not Available]", "[Not Applicable]", ""))
clin$OS      <- suppressWarnings(as.numeric(clin$OS))
clin$OS.time <- suppressWarnings(as.numeric(clin$OS.time))
# Drop patients GDC/TCGA flagged for redaction, and rows missing OS or OS.time
clin <- clin[(is.na(clin$Redaction) | clin$Redaction == "") &
             !is.na(clin$OS) & !is.na(clin$OS.time), ]

# ── Program definitions ───────────────────────────────────────────────────
cohorts <- list(
  BRCA = c("MP7", "MP10"),
  KIRC = c("MP11", "MP13", "MP15")
)

score_program <- function(matz, genes) {
  present <- intersect(genes, colnames(matz))
  if (length(present) < 0.5 * length(genes)) {
    warning(sprintf("only %d/%d program genes found in expression matrix", length(present), length(genes)))
  }
  rowMeans(matz[, present, drop = FALSE], na.rm = TRUE)
}

build_matrix <- function(expr_long) {
  patients <- unique(expr_long$patient_barcode)
  genes    <- unique(expr_long$gene)
  mat <- matrix(NA_real_, nrow = length(patients), ncol = length(genes),
                dimnames = list(patients, genes))
  mat[cbind(match(expr_long$patient_barcode, patients),
            match(expr_long$gene, genes))] <- log2(expr_long$tpm + 1)
  scale(mat)  # per-gene z-score across patients
}

results <- list()

for (cohort in names(cohorts)) {

  expr_long <- read.csv(file.path(DATA_DIR, sprintf("%s_expr_long.csv", cohort)),
                         stringsAsFactors = FALSE)
  matz <- build_matrix(expr_long)

  cohort_clin <- clin[clin$type == cohort, c("bcr_patient_barcode", "OS", "OS.time")]
  common <- intersect(rownames(matz), cohort_clin$bcr_patient_barcode)
  log_msg("%s: %d patients with expression, %d with matched clinical data",
          cohort, nrow(matz), length(common))

  df <- merge(cohort_clin, data.frame(bcr_patient_barcode = rownames(matz)), by = "bcr_patient_barcode")
  df <- df[match(common, df$bcr_patient_barcode), ]

  for (program in cohorts[[cohort]]) {
    genes <- readLines(file.path(GENE_DIR, sprintf("%s_%s.txt", cohort, program)))
    score <- score_program(matz, genes)
    df$score <- score[match(df$bcr_patient_barcode, names(score))]

    # Median split
    med <- median(df$score, na.rm = TRUE)
    df$grp_median <- factor(ifelse(df$score > med, "High", "Low"), levels = c("Low", "High"))

    # Tertile split (drop middle third) as a robustness check
    terts <- quantile(df$score, c(1/3, 2/3), na.rm = TRUE)
    df$grp_tertile <- ifelse(df$score <= terts[1], "Low",
                       ifelse(df$score >= terts[2], "High", NA))
    df$grp_tertile <- factor(df$grp_tertile, levels = c("Low", "High"))

    sd_med <- survdiff(Surv(OS.time, OS) ~ grp_median, data = df)
    p_med  <- 1 - pchisq(sd_med$chisq, df = 1)

    df_tert <- df[!is.na(df$grp_tertile), ]
    sd_tert <- survdiff(Surv(OS.time, OS) ~ grp_tertile, data = df_tert)
    p_tert  <- 1 - pchisq(sd_tert$chisq, df = 1)

    cox_cont <- coxph(Surv(OS.time, OS) ~ score, data = df)
    cox_sum  <- summary(cox_cont)
    hr       <- cox_sum$coefficients[1, "exp(coef)"]
    hr_lo    <- cox_sum$conf.int[1, "lower .95"]
    hr_hi    <- cox_sum$conf.int[1, "upper .95"]
    p_cox    <- cox_sum$coefficients[1, "Pr(>|z|)"]

    n_events <- sum(df$OS == 1, na.rm = TRUE)

    results[[length(results) + 1]] <- data.frame(
      cohort = cohort, program = program, n = nrow(df), n_events = n_events,
      n_genes_used = sum(genes %in% colnames(matz)), n_genes_total = length(genes),
      p_logrank_median = p_med, p_logrank_tertile = p_tert,
      cox_HR_per_1SD = hr, cox_HR_lo95 = hr_lo, cox_HR_hi95 = hr_hi, p_cox = p_cox
    )

    log_msg("  %s/%s: n=%d events=%d | median-split logrank p=%.4f | tertile-split p=%.4f | Cox HR=%.2f [%.2f-%.2f] p=%.4f",
            cohort, program, nrow(df), n_events, p_med, p_tert, hr, hr_lo, hr_hi, p_cox)

    # Only plot if there's a real signal on either split (p < 0.1)
    if (min(p_med, p_tert, na.rm = TRUE) < 0.1) {
      png(file.path(KM_DIR, sprintf("%s_%s_km.png", cohort, program)),
          width = 1600, height = 1400, res = 220)
      fit <- survfit(Surv(OS.time, OS) ~ grp_median, data = df)
      plot(fit, col = c("steelblue", "firebrick"), lwd = 2.5,
           xlab = "Days", ylab = "Overall survival probability",
           main = sprintf("%s — %s (median split)\nlog-rank p = %.4f, Cox HR = %.2f [%.2f-%.2f]",
                           cohort, program, p_med, hr, hr_lo, hr_hi))
      legend("bottomleft", legend = c(sprintf("Low (n=%d)", sum(df$grp_median == "Low")),
                                       sprintf("High (n=%d)", sum(df$grp_median == "High"))),
             col = c("steelblue", "firebrick"), lwd = 2.5, bty = "n")
      dev.off()
      log_msg("  -> KM plot saved: %s", file.path(KM_DIR, sprintf("%s_%s_km.png", cohort, program)))
    }
  }
}

out <- do.call(rbind, results)
write.csv(out, file.path(OUT_DIR, "tcga_survival_results.csv"), row.names = FALSE)
log_msg("Done. Results table: %s", file.path(OUT_DIR, "tcga_survival_results.csv"))
