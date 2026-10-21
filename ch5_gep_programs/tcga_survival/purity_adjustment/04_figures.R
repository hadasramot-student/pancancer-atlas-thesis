#!/usr/bin/env Rscript
# 04_figures.R -- forest plot of the MP11 HRs + KM (CPE>=0.7 median split)
suppressPackageStartupMessages(library(survival))
PD <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical/purity_adjustment"
df <- readRDS(file.path(PD, "kirc_mp11_analysis_frame.rds"))

fp <- read.csv(file.path(PD, "mp11_purity_survival_results_v2.csv"), stringsAsFactors = FALSE)
keep <- c("unadj_zMP11","adj_stageCovar","adj_strataStage_PRIMARY","adj_strataStage_ESTIMATE",
          "CPE>=0.6_adj_strata","CPE>=0.7_adj_strata","MW_unadj","MW_adj_strataStage")
lab <- c("MP11 unadjusted (n=510)","MP11 + age/sex/stage/grade/CPE (n=498)",
         "MP11 + age/sex/strata(stage)/grade/CPE  [primary] (n=498)",
         "  \" , ESTIMATE purity instead of CPE (n=499)",
         "MP11 adjusted, CPE>=0.6 (n=331)","MP11 adjusted, CPE>=0.7 (n=204)",
         "Malignant-weighted (C19orf33+CRYAB) unadjusted (n=510)",
         "Malignant-weighted adjusted (n=498)")
fp <- fp[match(keep, fp$model), ]

png(file.path(PD, "mp11_forest.png"), width = 1900, height = 1050, res = 200)
par(mar = c(4.5, 0.5, 3, 1))
n <- nrow(fp); ys <- n:1
xlim <- c(0.4, 1.15)
plot(NA, xlim = xlim, ylim = c(0.5, n + 0.5), axes = FALSE, xlab = "Hazard ratio (per 1 SD)  —  log scale", ylab = "", log = "x")
abline(v = 1, col = "grey60", lty = 2)
axis(1, at = c(0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1))
segments(fp$lo95, ys, fp$hi95, ys, lwd = 2)
pch_ <- ifelse(grepl("^MW|Malignant", lab), 22, 15)
points(fp$HR, ys, pch = pch_, cex = 1.6, bg = "white")
text(rep(xlim[1], n), ys, lab, pos = 4, xpd = NA, cex = 0.8, offset = 0)
text(rep(xlim[2], n), ys,
     sprintf("%.2f (%.2f-%.2f)  p=%.1e", fp$HR, fp$lo95, fp$hi95, fp$p),
     pos = 2, xpd = NA, cex = 0.75)
mtext("TCGA-KIRC overall survival — MP11 score, purity adjustment", 3, line = 1, cex = 1.05, font = 2)
dev.off()

# KM: CPE>=0.7, MP11 median split
s <- df[!is.na(df$CPE) & df$CPE >= 0.7, ]
s$grp <- factor(ifelse(s$score_orig > median(s$score_orig), "High MP11", "Low MP11"), c("Low MP11","High MP11"))
sd <- survdiff(Surv(OS.time, OS) ~ grp, s); p <- 1 - pchisq(sd$chisq, 1)
png(file.path(PD, "mp11_km_CPE0.7.png"), width = 1300, height = 1150, res = 200)
plot(survfit(Surv(OS.time, OS) ~ grp, s), col = c("steelblue","firebrick"), lwd = 3,
     xlab = "Days", ylab = "Overall survival",
     main = sprintf("TCGA-KIRC, CPE >= 0.7 (n=%d, %d events)\nMP11 median split — log-rank p = %.4f", nrow(s), sum(s$OS), p))
legend("bottomleft", c(sprintf("Low MP11 (n=%d)", sum(s$grp=="Low MP11")),
                        sprintf("High MP11 (n=%d)", sum(s$grp=="High MP11"))),
       col = c("steelblue","firebrick"), lwd = 3, bty = "n")
dev.off()
cat("wrote mp11_forest.png, mp11_km_CPE0.7.png\n")
