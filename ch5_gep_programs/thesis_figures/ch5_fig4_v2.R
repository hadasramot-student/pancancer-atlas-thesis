#!/usr/bin/env Rscript
# ch5_fig4_v2.R  (comment 15/16 follow-up, 2026-10-21)
# Identical to Part 1 of ch5_fig4_and_part3.R except the Cox text line in the panel:
# unadjusted Cox HR per 1 SD of the MP11 score (same scale as the adjusted HR 0.71).
# Writes ch5_fig4_v2.png/.pdf only; nothing existing is overwritten.
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
say("PART1 SD of score_orig = %.4f ; Cox per 1 SD of score HR=%.4f [%.4f-%.4f] p=%.3g", sd_score, hr1[1], hr1[2], hr1[3], cx1$coefficients[1, "Pr(>|z|)"])
say("PART1 printed: Cox HR = %.2f [%.2f-%.2f] per 1 SD (unadjusted)", hr1[1], hr1[2], hr1[3])

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
       sprintf("Cox HR = %.2f [%.2f–%.2f] per 1 SD (unadjusted)", hr1[1], hr1[2], hr1[3]))
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
png(file.path(OUT, "ch5_fig4_v2.png"), width = W, height = H, units = "in", res = 400); draw(); dev.off()
cairo_pdf(file.path(OUT, "ch5_fig4_v2.pdf"), width = W, height = H); draw(); dev.off()
say("PART1 at-risk years %s", paste(brk, collapse = ","))
say("PART1 at-risk High %s", paste(rk$High, collapse = ","))
say("PART1 at-risk Low  %s", paste(rk$Low,  collapse = ","))
say("PART1 max follow-up years = %.2f", max(df$OS.time) / 365.25)
say("PART1 v2 figures written")

