suppressPackageStartupMessages({ library(Seurat); library(data.table) })
rds <- "/ofir_grp/users/hadasra/temp_objects/clustered_lineages/Epithelial_Kidney_clustered.RDS"
out <- "/ofir_grp/users/hadasra/temp_objects/clustered_lineages/kidney_c5_c9_check"
dir.create(out, showWarnings = FALSE)
obj <- readRDS(rds)
obj <- JoinLayers(obj)
md <- obj@meta.data
cat("meta cols:\n"); print(colnames(md))
cl <- as.character(md$louvain_res0.50)

genes <- c("IGHM","JCHAIN","FCN1","IGKC","IGHA1","IGHG1","MZB1","CD79A","PTPRC","LYZ","CD14","S100A8","S100A9",
           "HBA1","HBA2","HBB","ALAS2","HBD",
           "LRP2","CUBN","SLC34A1","ALDOB","GPX3","MIOX","CA9","NDUFA4L2","EPCAM","KRT8","KRT18","PAX8")
genes <- intersect(genes, rownames(obj))
# counts layer = SCimilarity log1p-normalized (single normalization); data = NormalizeData on top
L <- list(counts = LayerData(obj, "counts")[genes, ], data = LayerData(obj, "data")[genes, ])
for (ln in names(L)) {
  m <- as.matrix(L[[ln]])
  rows <- list()
  for (g in genes) for (grp in c("5", "9", "rest")) {
    idx <- if (grp == "rest") !(cl %in% c("5","9")) else cl == grp
    v <- m[g, idx]
    q <- quantile(v, c(.25,.5,.75,.9,.95,.99))
    rows[[length(rows)+1]] <- data.table(layer=ln, gene=g, group=grp, n=length(v),
      frac_gt0=mean(v>0), frac_gt0.5=mean(v>0.5), frac_gt1=mean(v>1), frac_gt2=mean(v>2), frac_gt3=mean(v>3),
      mean=mean(v), mean_if_pos=if (any(v>0)) mean(v[v>0]) else NA,
      q25=q[1], q50=q[2], q75=q[3], q90=q[4], q95=q[5], q99=q[6],
      share_signal_top10pct = { s <- sort(v, decreasing=TRUE); if (sum(s)>0) sum(s[seq_len(ceiling(.1*length(s)))])/sum(s) else NA })
  }
  fwrite(rbindlist(rows), file.path(out, paste0("gene_dist_", ln, ".tsv")), sep="\t")
}
m <- as.matrix(L$counts)
anyset <- function(gs, thr) colSums(m[intersect(gs, genes), , drop=FALSE] > thr) > 0
ig <- c("IGHM","JCHAIN","FCN1"); hb <- c("HBA1","HBA2","HBB")
pt <- c("LRP2","CUBN","SLC34A1","ALDOB","GPX3","MIOX")
comb <- rbindlist(lapply(c("5","9","rest"), function(grp) {
  idx <- if (grp=="rest") !(cl %in% c("5","9")) else cl==grp
  data.table(group=grp, n=sum(idx),
    ig_any_gt0=mean(anyset(ig,0)[idx]), ig_any_gt1=mean(anyset(ig,1)[idx]),
    ig_all3_gt0=mean((colSums(m[ig,,drop=FALSE]>0)==3)[idx]),
    hb_any_gt0=mean(anyset(hb,0)[idx]), hb_any_gt1=mean(anyset(hb,1)[idx]),
    hb_all3_gt0=mean((colSums(m[hb,,drop=FALSE]>0)==3)[idx]),
    pt_any_gt0=mean(anyset(pt,0)[idx]), PTPRC_gt0=mean((m["PTPRC",]>0)[idx]),
    EPCAM_or_KRT8_gt0=mean(anyset(c("EPCAM","KRT8","KRT18"),0)[idx]),
    pt_pos_among_ig_pos = mean(anyset(pt,0)[idx & anyset(ig,0)]),
    ptprc_pos_among_ig_pos = mean((m["PTPRC",]>0)[idx & anyset(ig,0)]),
    pt_pos_among_hb_pos = mean(anyset(pt,0)[idx & anyset(hb,0)]),
    alas2_pos_among_hb_pos = if ("ALAS2" %in% genes) mean((m["ALAS2",]>0)[idx & anyset(hb,0)]) else NA)
}))
fwrite(comb, file.path(out, "combined_sets.tsv"), sep="\t")

sc <- data.table(cluster=cl, ig=colSums(m[ig,]), hb=colSums(m[hb,]), pt=colMeans(m[pt,]),
                 ig_pos=anyset(ig,0), hb_pos=anyset(hb,0))
fwrite(sc[, .(n=.N, ig_score=mean(ig), ig_frac_pos=mean(ig_pos), hb_score=mean(hb), hb_frac_pos=mean(hb_pos), pt_score=mean(pt)),
          by=cluster][order(as.integer(cluster))], file.path(out, "per_cluster_scores.tsv"), sep="\t")

for (col in intersect(c("project_name","sample","orig.ident","sample_id","patient","dataset","repository","Consensus_Cell_Type","scimilarity_celltype","meta_type"), colnames(md))) {
  d <- data.table(cluster=cl, v=as.character(md[[col]]), ig_pos=anyset(ig,0), hb_pos=anyset(hb,0))
  tab <- d[cluster %in% c("5","9"), .(n=.N, ig_frac_pos=mean(ig_pos), hb_frac_pos=mean(hb_pos)), by=.(cluster, v)]
  tot <- d[, .(n_in_object=.N, obj_ig_frac_pos=mean(ig_pos), obj_hb_frac_pos=mean(hb_pos)), by=v]
  fwrite(merge(tab, tot, by="v")[order(cluster,-n)], file.path(out, paste0("by_", col, ".tsv")), sep="\t")
}
num <- names(md)[sapply(md, is.numeric)]
q <- as.data.table(md[, num, drop=FALSE]); q$cluster <- cl
fwrite(q[, lapply(.SD, median, na.rm=TRUE), by=cluster][order(as.integer(cluster))], file.path(out, "per_cluster_numeric_medians.tsv"), sep="\t")
cat("done\n")
