#!/usr/bin/env Rscript
# 01_consolidated_pt_survival.R
# Ofir Point 15, Fix 3 (rerun branch): KIRC MP11 OS analysis re-run on the
# CONSOLIDATED proximal-tubule program containing MP11 = {MP2,MP6,MP8,MP11}
# (Point 11f union-find, top-50 Jaccard >= 0.25). Consolidated signature =
# UNION of member top-50 lists (primary); intersection & ">=2 of 4" as
# sensitivity checks. Scoring + framings identical to run_tcga_survival.R and
# the Fix-1 purity-adjusted model. Expression re-assembled from raw GDC STAR
# gene-count TSVs (tpm_unstranded). /ofir_grp == /groups/ofircohen-group
suppressPackageStartupMessages(library(survival))
set.seed(1)

DD    <- "/ofir_grp/users/hadasra/temp_objects/tcga_clinical"
PD    <- file.path(DD, "purity_adjustment")
OD    <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/kirc-survival-consolidated-pt-da69f5/consolidated_pt"
STAR  <- "/ofir_grp/datasets/RNA-SEQ-TCGA-PANCANCER/RNA-SEQ"
GENNMF<- "/ofir_grp/users/hadasra/temp_objects/genenmf_output/Kidney/Kidney_genenmf_MP_top50_genes.csv"
CACHE <- file.path(OD, "kirc_star_tpm_matrix.rds")
dir.create(OD, showWarnings = FALSE, recursive = TRUE)
say <- function(...) cat(sprintf(...), "\n")

GROUP    <- c("MP2","MP6","MP8","MP11")
JACC_THR <- 0.25

## 1. Signatures --------------------------------------------------------------
g   <- read.csv(GENNMF, stringsAsFactors = FALSE)
ord <- order(g$meta_program, g$rank)
mps <- split(as.character(g$gene[ord]), g$meta_program[ord])
mps <- lapply(mps, unique)
jacc <- function(a,b){a<-unique(a);b<-unique(b);length(intersect(a,b))/length(union(a,b))}
say("=== Point 11f check: pairwise top-50 Jaccard within {%s} ===", paste(GROUP,collapse=","))
for (i in 1:(length(GROUP)-1)) for (j in (i+1):length(GROUP)) {
  a<-GROUP[i]; b<-GROUP[j]; jv<-jacc(mps[[a]],mps[[b]])
  say("  J(%s,%s) = %.3f%s", a, b, jv, ifelse(jv>=JACC_THR,"  (>= thr, linked)",""))
}
memb <- mps[GROUP]; tab <- table(unlist(memb))
sig_union <- names(tab)
sig_freq2 <- names(tab)[tab>=2]
sig_inter <- names(tab)[tab==length(GROUP)]
mp11_only <- mps[["MP11"]]
say("\nsignature sizes:  MP11-alone=%d  union=%d  freq>=2=%d  intersection=%d",
    length(mp11_only), length(sig_union), length(sig_freq2), length(sig_inter))
say("intersection genes: %s", paste(sig_inter, collapse=", "))
writeLines(sig_union, file.path(OD,"sig_consolidated_union.txt"))
writeLines(sig_freq2, file.path(OD,"sig_consolidated_freq2.txt"))
writeLines(sig_inter, file.path(OD,"sig_consolidated_intersection.txt"))
SIGS <- list(MP11_alone=mp11_only, consolidated_union=sig_union,
             consolidated_freq2=sig_freq2, consolidated_inter=sig_inter)
want_genes <- sort(unique(unlist(SIGS)))

## 2. KIRC bulk expression (cache-aware) -----------------------------------
matched <- read.csv(file.path(DD,"KIRC_matched.csv"), header=FALSE, stringsAsFactors=FALSE)
colnames(matched) <- c("patient","file")
say("\nKIRC matched patients (bulk STAR files): %d", nrow(matched))

if (file.exists(CACHE)) {
  mat <- readRDS(CACHE)
  if (!all(want_genes %in% colnames(mat)) || !all(matched$patient %in% rownames(mat))) {
    say("cache missing genes/patients -> rebuilding"); mat <- NULL
  } else { say("loaded cached TPM matrix: %d x %d", nrow(mat), ncol(mat)) }
} else mat <- NULL

if (is.null(mat)) {
  say("distinct signature genes to pull from STAR TSVs: %d", length(want_genes))
  read_one <- function(fp){
    x <- tryCatch(read.delim(fp, comment.char="#", stringsAsFactors=FALSE), error=function(e) NULL)
    if (is.null(x)) return(NULL)
    x <- x[!startsWith(x$gene_id,"N_"),]
    x <- x[x$gene_name %in% want_genes, c("gene_name","tpm_unstranded")]
    tapply(as.numeric(x$tpm_unstranded), x$gene_name, sum)
  }
  mat <- matrix(NA_real_, nrow=nrow(matched), ncol=length(want_genes),
                dimnames=list(matched$patient, want_genes))
  n_ok <- 0
  for (k in seq_len(nrow(matched))) {
    v <- read_one(file.path(STAR, matched$file[k]))
    if (is.null(v)) { say("  !! could not read %s", matched$file[k]); next }
    mat[k, names(v)] <- v; n_ok <- n_ok + 1
  }
  say("STAR TSVs read OK: %d / %d", n_ok, nrow(matched))
  saveRDS(mat, CACHE)
}
miss_g <- want_genes[colSums(!is.na(mat))==0]
if (length(miss_g)) say("genes absent from STAR annotation: %s", paste(miss_g, collapse=", "))

tryCatch({
  el <- read.csv(file.path(DD,"KIRC_expr_long.csv"), stringsAsFactors=FALSE)
  shared <- intersect(unique(el$gene), colnames(mat))
  el <- el[el$gene %in% shared,]
  w <- reshape(el, idvar="patient_barcode", timevar="gene", direction="wide")
  rownames(w) <- w$patient_barcode; w$patient_barcode <- NULL
  colnames(w) <- sub("^tpm\\.","",colnames(w))
  cp <- intersect(rownames(w), rownames(mat))
  d  <- sapply(shared, function(gn) max(abs(w[cp,gn]-mat[cp,gn]), na.rm=TRUE))
  rel<- sapply(shared, function(gn){a<-w[cp,gn];b<-mat[cp,gn];max(abs(a-b)/pmax(a,b,1e-6),na.rm=TRUE)})
  say("\nSANITY vs cached KIRC_expr_long.csv: %d shared genes, %d patients", length(shared), length(cp))
  say("  max |dTPM| = %.5g   max relative diff = %.5g", max(d), max(rel))
}, error=function(e) say("sanity check skipped: %s", conditionMessage(e)))

logmat <- log2(mat + 1)
matz   <- scale(logmat)
score_program <- function(genes){
  present <- intersect(genes, colnames(matz))
  present <- present[colSums(!is.na(matz[,present,drop=FALSE]))>0]
  list(score=rowMeans(matz[,present,drop=FALSE], na.rm=TRUE),
       n_used=length(present), n_tot=length(genes))
}

## 3. Analysis frame ------------------------------------------------------
fr <- readRDS(file.path(PD,"kirc_mp11_analysis_frame.rds"))
fr$grade3 <- factor(fr$grade3, levels=c("G1-2","G3","G4"))
fr$stage  <- factor(fr$stage,  levels=c("I","II","III","IV"))
say("\nFix-1 analysis frame: n=%d  events=%d  CPE non-NA=%d",
    nrow(fr), sum(fr$OS), sum(!is.na(fr$CPE)))
sc <- data.frame(bcr_patient_barcode=rownames(matz), stringsAsFactors=FALSE)
meta <- list()
for (nm in names(SIGS)) {
  s <- score_program(SIGS[[nm]])
  sc[[paste0("raw_",nm)]] <- s$score[sc$bcr_patient_barcode]
  sc[[paste0("z_",nm)]]   <- as.numeric(scale(sc[[paste0("raw_",nm)]]))
  meta[[nm]] <- c(n_used=s$n_used, n_tot=s$n_tot)
}
df <- merge(fr, sc, by="bcr_patient_barcode", all.x=TRUE)
say("frame after merge: n=%d", nrow(df))
say("\nscore diagnostics (Spearman):")
for (nm in names(SIGS)) {
  r_mp11 <- cor(df[[paste0("raw_",nm)]], df$raw_MP11_alone, method="spearman", use="complete.obs")
  r_orig <- cor(df[[paste0("raw_",nm)]], df$score_orig,     method="spearman", use="complete.obs")
  r_cpe  <- cor(df[[paste0("raw_",nm)]], df$CPE,            method="spearman", use="complete.obs")
  say("  %-20s genes %2d/%2d   r(.,MP11new)=%+.3f  r(.,score_orig)=%+.3f  r(.,CPE)=%+.3f",
      nm, meta[[nm]]["n_used"], meta[[nm]]["n_tot"], r_mp11, r_orig, r_cpe)
}

## 4. Survival ----------------------------------------------------------
res <- list()
f0  <- function(x) if (is.null(x) || length(x)==0 || all(is.na(x))) NA_real_ else unname(x)[1]
push <- function(sig,framing,n,ev,HR,lo,hi,p,extra=""){
  res[[length(res)+1]] <<- data.frame(
    signature=sig, framing=framing, n=f0(n), events=f0(ev),
    HR=f0(HR), lo95=f0(lo), hi95=f0(hi), p=f0(p),
    note=if (length(extra)==0) "" else as.character(extra)[1],
    stringsAsFactors=FALSE)
}
cox_row <- function(m,term){
  s<-summary(m); ci<-s$conf.int[term,]
  c(HR=unname(ci["exp(coef)"]),lo=unname(ci["lower .95"]),hi=unname(ci["upper .95"]),
    p=unname(s$coefficients[term,ncol(s$coefficients)]),n=s$n,ev=s$nevent)
}
for (nm in names(SIGS)) {
  d <- df; d$raw <- d[[paste0("raw_",nm)]]; d$z <- d[[paste0("z_",nm)]]
  d <- d[!is.na(d$raw),]
  r <- cox_row(coxph(Surv(OS.time,OS)~raw,d),"raw")
  push(nm,"continuous Cox (raw mean-z, orig parametrisation)",r["n"],r["ev"],r["HR"],r["lo"],r["hi"],r["p"])
  r <- cox_row(coxph(Surv(OS.time,OS)~z,d),"z")
  push(nm,"continuous Cox (per 1 SD)",r["n"],r["ev"],r["HR"],r["lo"],r["hi"],r["p"])
  med <- median(d$raw,na.rm=TRUE)
  d$grp_med <- factor(ifelse(d$raw>med,"High","Low"),c("Low","High"))
  p_m <- 1-pchisq(survdiff(Surv(OS.time,OS)~grp_med,d)$chisq,1)
  r <- cox_row(coxph(Surv(OS.time,OS)~grp_med,d),"grp_medHigh")
  push(nm,"median-split log-rank",sum(!is.na(d$grp_med)),sum(d$OS[!is.na(d$grp_med)]),r["HR"],r["lo"],r["hi"],p_m,"HR High vs Low")
  q <- quantile(d$raw,c(1/3,2/3),na.rm=TRUE)
  d$grp_t <- factor(ifelse(d$raw<=q[1],"Low",ifelse(d$raw>=q[2],"High",NA)),c("Low","High"))
  dt <- d[!is.na(d$grp_t),]
  p_t <- 1-pchisq(survdiff(Surv(OS.time,OS)~grp_t,dt)$chisq,1)
  r <- cox_row(coxph(Surv(OS.time,OS)~grp_t,dt),"grp_tHigh")
  push(nm,"tertile-split log-rank (drop middle)",nrow(dt),sum(dt$OS),r["HR"],r["lo"],r["hi"],p_t,"HR High vs Low")
  m_adj <- coxph(Surv(OS.time,OS)~z+age+sex+strata(stage)+grade3+CPE,d)
  r <- cox_row(m_adj,"z")
  push(nm,"purity-adj Cox: z + age + sex + strata(stage) + grade3 + CPE  [PRIMARY]",r["n"],r["ev"],r["HR"],r["lo"],r["hi"],r["p"],"per 1 SD")
  zt <- cox.zph(m_adj)$table
  push(nm,"  cox.zph signature-term p",r["n"],r["ev"],NA,NA,NA,zt[grep("^z$",rownames(zt)),"p"],"PH diagnostic")
  r <- cox_row(coxph(Surv(OS.time,OS)~z+age+sex+stage+grade3+CPE,d),"z")
  push(nm,"purity-adj Cox: z + age + sex + stage(covar) + grade3 + CPE",r["n"],r["ev"],r["HR"],r["lo"],r["hi"],r["p"],"per 1 SD")
  cc <- d[complete.cases(d[,c("OS.time","OS","z","age","sex","stage","grade3","CPE")]),]
  lr <- anova(coxph(Surv(OS.time,OS)~age+sex+strata(stage)+grade3+CPE,cc),
              coxph(Surv(OS.time,OS)~z+age+sex+strata(stage)+grade3+CPE,cc))
  pcol  <- grep("^Pr\\(", names(lr), value=TRUE)
  p_lr  <- if (length(pcol)) lr[[pcol]][2] else NA
  chi_lr<- if (!is.null(lr$Chisq)) lr$Chisq[2] else NA
  push(nm,"  LR test add signature to covariates p",nrow(cc),sum(cc$OS),NA,NA,NA,p_lr,
       sprintf("chisq=%.2f df=1", ifelse(is.null(chi_lr),NA,chi_lr)))
}
out <- do.call(rbind,res)
out[,c("HR","lo95","hi95")] <- lapply(out[,c("HR","lo95","hi95")],round,3)
out$p <- signif(out$p,3)
write.csv(out, file.path(OD,"consolidated_pt_survival_results.csv"), row.names=FALSE)
say("\n================  FULL RESULTS  ================")
print(out, row.names=FALSE)
mk <- function(sig,fr_key){
  x <- out[out$signature==sig & out$framing==fr_key,]
  if (!nrow(x)) return(c(n=NA,HR=NA,CI=NA,p=NA))
  c(n=x$n[1],HR=sprintf("%.2f",x$HR[1]),CI=sprintf("%.2f-%.2f",x$lo95[1],x$hi95[1]),p=format(signif(x$p[1],3),scientific=TRUE))
}
say("\n======  MP11-ALONE  vs  CONSOLIDATED {MP2,MP6,MP8,MP11}  ======")
say("(original MP11-alone thesis result: n=510  HR 0.58  CI 0.45-0.74  p<0.0001)")
frs <- c("continuous Cox (raw mean-z, orig parametrisation)","continuous Cox (per 1 SD)",
         "median-split log-rank","tertile-split log-rank (drop middle)",
         "purity-adj Cox: z + age + sex + strata(stage) + grade3 + CPE  [PRIMARY]")
for (f in frs) {
  say("\n-- %s", f)
  for (s in names(SIGS)) { v <- mk(s,f)
    say("   %-20s n=%-4s HR=%-6s CI %-12s p=%s", s, v["n"], v["HR"], v["CI"], v["p"]) }
}
saveRDS(df, file.path(OD,"kirc_consolidated_pt_analysis_frame.rds"))
d <- df[!is.na(df$raw_consolidated_union),]
d$grp <- factor(ifelse(d$raw_consolidated_union>median(d$raw_consolidated_union),"High","Low"),c("Low","High"))
pm <- 1-pchisq(survdiff(Surv(OS.time,OS)~grp,d)$chisq,1)
cx <- summary(coxph(Surv(OS.time,OS)~grp,d))$conf.int[1,]
png(file.path(OD,"km_consolidated_union_median.png"), width=1600, height=1400, res=220)
plot(survfit(Surv(OS.time,OS)~grp,d), col=c("steelblue","firebrick"), lwd=2.5,
     xlab="Days", ylab="Overall survival probability",
     main=sprintf("TCGA-KIRC consolidated PT program {MP2,MP6,MP8,MP11}\nmedian split: log-rank p = %.4f, HR = %.2f [%.2f-%.2f]",
                  pm, cx["exp(coef)"], cx["lower .95"], cx["upper .95"]))
legend("bottomleft", c(sprintf("Low (n=%d)",sum(d$grp=="Low")), sprintf("High (n=%d)",sum(d$grp=="High"))),
       col=c("steelblue","firebrick"), lwd=2.5, bty="n")
dev.off()
say("\n[%s] DONE -> %s", Sys.time(), OD)
