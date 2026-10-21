#!/usr/bin/env Rscript
# 03_purity_survival_v2.R  -- clean re-fit
#   * grade collapsed to G1-2 / G3 / G4  (4-level grade -> monotone-likelihood
#     separation because KIRC has only 14 G1 patients, ~0 events)
#   * stage: AJCC I/II/III/IV violates PH (cox.zph p=2e-4) -> PRIMARY adjusted
#     model stratifies on stage; a covariate-in version is reported alongside
#   * adds: KM median-split log-rank inside CPE>=0.7 (thesis' own framing),
#     ESTIMATE-purity sensitivity, per-stratum checks.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(survival))
DD <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical"
PD <- file.path(DD, "purity_adjustment")
say <- function(...) cat(sprintf(...), "\n")
line <- function(m, term) {
  s <- summary(m); ci <- s$conf.int[term, ]; p <- s$coefficients[term, ncol(s$coefficients)]
  sprintf("HR=%.3f  95%% CI %.3f-%.3f  p=%.3g   (n=%d, ev=%d)",
          ci["exp(coef)"], ci["lower .95"], ci["upper .95"], p, s$n, s$nevent)
}

df <- readRDS(file.path(PD, "kirc_mp11_analysis_frame.rds"))
say("frame: n=%d ev=%d ; CPE non-NA=%d", nrow(df), sum(df$OS), sum(!is.na(df$CPE)))
# grade3 / stage already on the frame from script 02; recompute defensively
df$grade3  <- factor(df$grade3,  levels = c("G1-2","G3","G4"))
df$stage   <- factor(df$stage,   levels = c("I","II","III","IV"))
df$stage34 <- factor(df$stage34, levels = c("I-II","III-IV"))

res <- list()
add <- function(tag, m, term = "z_mp11") {
  s <- summary(m); ci <- s$conf.int[term, ]
  res[[length(res)+1]] <<- data.frame(model = tag, term = term, n = s$n, events = s$nevent,
    HR = ci["exp(coef)"], lo95 = ci["lower .95"], hi95 = ci["upper .95"],
    p = s$coefficients[term, ncol(s$coefficients)], row.names = NULL)
  say("[%-40s] %s", tag, line(m, term))
}

say("\n===== UNADJUSTED (reference) =====")
add("unadj_zMP11", coxph(Surv(OS.time,OS) ~ z_mp11, df))

say("\n===== ADJUSTED, stage as COVARIATE (grade3) =====")
m_cov <- coxph(Surv(OS.time,OS) ~ z_mp11 + age + sex + stage + grade3 + CPE, df)
print(summary(m_cov)); add("adj_stageCovar", m_cov)
say("-- cox.zph --"); print(cox.zph(m_cov))

say("\n===== ADJUSTED, PRIMARY: strata(stage) + grade3 + CPE =====")
m_str <- coxph(Surv(OS.time,OS) ~ z_mp11 + age + sex + strata(stage) + grade3 + CPE, df)
print(summary(m_str)); add("adj_strataStage_PRIMARY", m_str)
zph <- cox.zph(m_str); say("-- cox.zph (primary) --"); print(zph)
saveRDS(zph, file.path(PD, "coxzph_primary_strataStage.rds"))

say("\n-- LR test: covariate-only vs +MP11 (strata model, complete cases) --")
cc <- df[complete.cases(df[,c("OS.time","OS","z_mp11","age","sex","stage","grade3","CPE")]),]
print(anova(coxph(Surv(OS.time,OS) ~ age+sex+strata(stage)+grade3+CPE, cc),
            coxph(Surv(OS.time,OS) ~ z_mp11+age+sex+strata(stage)+grade3+CPE, cc)))

say("\n-- purity-only & purity-vs-MP11 --")
say("   CPE  alone : %s", line(coxph(Surv(OS.time,OS) ~ CPE, df), "CPE"))
say("   ESTIMATE   : %s", line(coxph(Surv(OS.time,OS) ~ ESTIMATE, df), "ESTIMATE"))
say("   Spearman r(MP11,CPE)=%.3f  r(MP11,ESTIMATE)=%.3f  r(MP11,ABSOLUTE)=%.3f",
    cor(df$score_orig,df$CPE,method="spearman",use="complete.obs"),
    cor(df$score_orig,df$ESTIMATE,method="spearman",use="complete.obs"),
    cor(df$score_orig,df$ABSOLUTE,method="spearman",use="complete.obs"))

say("\n-- ESTIMATE-purity adjustment (sensitivity) --")
add("adj_strataStage_ESTIMATE", coxph(Surv(OS.time,OS) ~ z_mp11+age+sex+strata(stage)+grade3+ESTIMATE, df))

say("\n===== HIGH-PURITY SUBSETS =====")
for (thr in c(0.6, 0.7)) {
  s <- df[!is.na(df$CPE) & df$CPE >= thr, ]
  say("--- CPE>=%.1f : n=%d ev=%d ---", thr, nrow(s), sum(s$OS))
  add(sprintf("CPE>=%.1f_unadj", thr),        coxph(Surv(OS.time,OS) ~ z_mp11, s))
  add(sprintf("CPE>=%.1f_adj_strata", thr),   coxph(Surv(OS.time,OS) ~ z_mp11+age+sex+strata(stage)+grade3, s))
  add(sprintf("CPE>=%.1f_adj_strata+CPE",thr),coxph(Surv(OS.time,OS) ~ z_mp11+age+sex+strata(stage)+grade3+CPE, s))
  # thesis' own framing: median split log-rank
  s$grp <- factor(ifelse(s$score_orig > median(s$score_orig), "High","Low"), c("Low","High"))
  sd <- survdiff(Surv(OS.time,OS) ~ grp, s)
  say("   median-split log-rank p = %.4g", 1 - pchisq(sd$chisq, 1))
}

say("\n===== MALIGNANT-WEIGHTED MP11 (2 genes kept: C19orf33, CRYAB) =====")
add("MW_unadj",          coxph(Surv(OS.time,OS) ~ z_mp11_mw, df), "z_mp11_mw")
add("MW_adj_strataStage", coxph(Surv(OS.time,OS) ~ z_mp11_mw+age+sex+strata(stage)+grade3+CPE, df), "z_mp11_mw")
for (thr in c(0.6,0.7)) {
  s <- df[!is.na(df$CPE) & df$CPE >= thr, ]
  add(sprintf("MW_CPE>=%.1f_adj", thr), coxph(Surv(OS.time,OS) ~ z_mp11_mw+age+sex+strata(stage)+grade3, s), "z_mp11_mw")
}

out <- do.call(rbind, res)
write.csv(out, file.path(PD, "mp11_purity_survival_results_v2.csv"), row.names = FALSE)
say("\nwrote mp11_purity_survival_results_v2.csv"); say("[%s] DONE", Sys.time())
