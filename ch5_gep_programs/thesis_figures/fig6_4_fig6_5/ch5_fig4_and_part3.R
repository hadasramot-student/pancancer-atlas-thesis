#!/usr/bin/env Rscript
# ch5_fig4_and_part3.R  (comment 15/16, 2026-10-21)
# Part 1: redraw Fig 6.4 (KIRC MP11 median-split KM) -> ch5_fig4.png/.pdf
# Part 3: purity-adjusted Cox for the MP11 Jaccard>=0.25 family {MP2,MP6,MP8,MP11}
#         (and, for comparison only, the {MP6,MP11,MP12,MP14} set named in the thesis text)
# Data loading reused from:
#   pancancer-atlas-thesis/ch5_gep_programs/tcga_survival/run_tcga_survival.R (original Fig)
#   temp_objects/tcga_clinical/purity_adjustment/03_purity_survival_v2.R (analysis frame)
#   temp_objects/tcga_clinical/consolidated_pt/01_consolidated_pt_survival.R (family scoring)
# Nothing existing is overwritten; all output goes to OUT.
suppressPackageStartupMessages(library(survival))

DD   <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical"
PD   <- file.path(DD, "purity_adjustment")
STAR <- "/ofir_grp/datasets/RNA-SEQ-TCGA-PANCANCER/RNA-SEQ"
GENNMF <- "/ofir_grp/users/hadasra/temp_objects/genenmf_output/Kidney/Kidney_genenmf_MP_top50_genes.csv"
CACHE0 <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/kirc-survival-consolidated-pt-da69f5/consolidated_pt/kirc_star_tpm_matrix.rds"
OUT  <- "/ofir_grp/users/hadasra/thesis_figures/comment15_16_2026-10-21"
say  <- function(...) cat(sprintf(...), "\n")

## ---------------------------------------------------------------- Part 1 ----
df <- readRDS(file.path(PD, "kirc_mp11_analysis_frame.rds"))
df$grade3 <- factor(df$grade3, levels = c("G1-2","G3","G4"))
df$stage  <- factor(df$stage,  levels = c("I","II","III","IV"))
med <- median(df$score_orig)
df$grp <- factor(ifelse(df$score_orig > med, "High", "Low"), levels = c("Low","High"))

sdf   <- survdiff(Surv(OS.time, OS) ~ grp, df)
p_lr  <- 1 - pchisq(sdf$chisq, 1)
cx    <- summary(coxph(Surv(OS.time, OS) ~ score_orig, df))
hr    <- cx$conf.int[1, c("exp(coef)","lower .95","upper .95")]
p_cox <- cx$coefficients[1, "Pr(>|z|)"]
sd_score <- sd(df$score_orig)
cx1   <- summary(coxph(Surv(OS.time, OS) ~ I(score_orig / sd_score), df))
hr1   <- cx1$conf.int[1, c("exp(coef)","lower .95","upper .95")]
ev    <- tapply(df$OS, df$grp, sum)
nn    <- table(df$grp)

say("PART1 n=%d events=%d", nrow(df), sum(df$OS))
say("PART1 Low: n=%d events=%d | High: n=%d events=%d", nn["Low"], ev["Low"], nn["High"], ev["High"])
say("PART1 log-rank chisq=%.3f p=%.2g", sdf$chisq, p_lr)
say("PART1 Cox (score as in original script: mean of per-gene z, 1 unit) HR=%.3f [%.3f-%.3f] p=%.3g",
    hr[1], hr[2], hr[3], p_cox)
say("PART1 SD of score_orig = %.3f ; Cox per 1 SD of score HR=%.3f [%.3f-%.3f]", sd_score, hr1[1], hr1[2], hr1[3])

fmt_p <- function(p) if (p < 1e-4) "p < 0.0001" else sprintf("p = %s", signif(p, 2))
fit <- survfit(Surv(OS.time / 365.25, OS) ~ grp, df)
COL <- c(Low = "#0072B2", High = "#C0392B")
brk <- seq(0, 12, 2)
risk <- summary(fit, times = brk, extend = TRUE)
rk   <- split(risk$n.risk, sub("grp=", "", as.character(risk$strata)))

draw <- function() {
  par(family = "sans", ps = 9, mar = c(7.2, 6.2, 0.6, 0.8), mgp = c(2.1, 0.55, 0),
      tcl = -0.3, las = 1, cex.axis = 8/9, cex.lab = 1)
  plot(fit, col = COL[c("Low","High")], lwd = 1.8, xlim = c(0, 12), xaxs = "i",
       ylim = c(0, 1), xaxt = "n", conf.int = FALSE, mark.time = TRUE, cex = 0.5,
       xlab = "", ylab = "Overall survival probability")
  axis(1, at = brk)
  mtext("Years", side = 1, line = 1.7, cex = 1)
  legend("topright", bty = "n", lwd = 1.8, col = COL[c("High","Low")], cex = 8/9,
         legend = c(sprintf("High MP11 score (n = %d)", nn["High"]),
                    sprintf("Low MP11 score (n = %d)",  nn["Low"])))
  text(0.2, 0.13, adj = c(0, 0.5), cex = 8/9,
       sprintf("Log-rank %s", fmt_p(p_lr)))
  text(0.2, 0.05, adj = c(0, 0.5), cex = 8/9,
       sprintf("Cox HR = %.2f [%.2f–%.2f] per unit score (unadjusted)", hr[1], hr[2], hr[3]))
  # numbers at risk
  mtext("Number at risk", side = 1, line = 3.1, adj = 0, at = -2.7, cex = 8/9, font = 2)
  mtext("High", side = 1, line = 4.2, at = -0.5, adj = 1, cex = 8/9, col = COL["High"], font = 2)
  mtext("Low",  side = 1, line = 5.2, at = -0.5, adj = 1, cex = 8/9, col = COL["Low"],  font = 2)
  for (i in seq_along(brk)) {
    mtext(rk$High[i], side = 1, line = 4.2, at = brk[i], cex = 8/9)
    mtext(rk$Low[i],  side = 1, line = 5.2, at = brk[i], cex = 8/9)
  }
}
W <- 5.4; H <- 3.9
png(file.path(OUT, "ch5_fig4.png"), width = W, height = H, units = "in", res = 400); draw(); dev.off()
cairo_pdf(file.path(OUT, "ch5_fig4.pdf"), width = W, height = H); draw(); dev.off()
say("PART1 at-risk years %s", paste(brk, collapse = ","))
say("PART1 at-risk High %s", paste(rk$High, collapse = ","))
say("PART1 at-risk Low  %s", paste(rk$Low,  collapse = ","))
say("PART1 max follow-up years = %.2f", max(df$OS.time) / 365.25)
say("PART1 figures written")

## ---------------------------------------------------------------- Part 3 ----
g   <- read.csv(GENNMF, stringsAsFactors = FALSE)
ord <- order(g$meta_program, g$rank)
mps <- lapply(split(as.character(g$gene[ord]), g$meta_program[ord]), unique)
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
GROUPS <- list(MP11_alone = "MP11",
               family_MP2_MP6_MP8_MP11  = c("MP2","MP6","MP8","MP11"),
               thesis_text_MP6_MP11_MP12_MP14 = c("MP6","MP11","MP12","MP14"))
for (nm in names(GROUPS)[-1]) {
  G <- GROUPS[[nm]]
  for (i in 1:(length(G)-1)) for (j in (i+1):length(G))
    say("PART3 %s J(%s,%s)=%.3f", nm, G[i], G[j], jacc(mps[[G[i]]], mps[[G[j]]]))
}
SIGS <- lapply(GROUPS, function(G) sort(unique(unlist(mps[G]))))   # UNION of member top-50 lists
want <- sort(unique(unlist(SIGS)))

matched <- read.csv(file.path(DD, "KIRC_matched.csv"), header = FALSE, stringsAsFactors = FALSE)
colnames(matched) <- c("patient","file")
mat <- readRDS(CACHE0)
need <- setdiff(want, colnames(mat))
say("PART3 genes needed beyond cached matrix: %d", length(need))
if (length(need)) {
  add <- matrix(NA_real_, nrow(matched), length(need), dimnames = list(matched$patient, need))
  for (k in seq_len(nrow(matched))) {
    x <- read.delim(file.path(STAR, matched$file[k]), comment.char = "#", stringsAsFactors = FALSE)
    x <- x[!startsWith(x$gene_id, "N_") & x$gene_name %in% need, c("gene_name","tpm_unstranded")]
    v <- tapply(as.numeric(x$tpm_unstranded), x$gene_name, sum)
    add[k, names(v)] <- v
  }
  mat <- cbind(mat[matched$patient, , drop = FALSE], add)
}
matz <- scale(log2(mat + 1))

cox_row <- function(m, term) {
  s <- summary(m); ci <- s$conf.int[term, ]
  sprintf("HR=%.3f [%.3f-%.3f] p=%.3g n=%d events=%d", ci[1], ci[3], ci[4],
          s$coefficients[term, ncol(s$coefficients)], s$n, s$nevent)
}
for (nm in names(SIGS)) {
  pres <- intersect(SIGS[[nm]], colnames(matz))
  pres <- pres[colSums(!is.na(matz[, pres, drop = FALSE])) > 0]
  raw  <- rowMeans(matz[, pres, drop = FALSE], na.rm = TRUE)
  d <- df; d$raw <- raw[d$bcr_patient_barcode]; d$z <- as.numeric(scale(d$raw))
  say("PART3 %s: genes %d/%d ; r_spearman(score, original MP11 score)=%.3f", nm, length(pres),
      length(SIGS[[nm]]), cor(d$raw, d$score_orig, method = "spearman"))
  say("PART3 %s unadjusted per 1 SD : %s", nm, cox_row(coxph(Surv(OS.time,OS) ~ z, d), "z"))
  m <- coxph(Surv(OS.time,OS) ~ z + age + sex + strata(stage) + grade3 + CPE, d)
  say("PART3 %s ADJUSTED per 1 SD (age+sex+strata(stage)+grade3+CPE): %s", nm, cox_row(m, "z"))
  say("PART3 %s cox.zph p(z)=%.3f", nm, cox.zph(m)$table["z", "p"])
  say("PART3 %s adjusted, stage as covariate: %s", nm,
      cox_row(coxph(Surv(OS.time,OS) ~ z + age + sex + stage + grade3 + CPE, d), "z"))
}
say("DONE")
