args <- commandArgs(trailingOnly=TRUE)
obj <- readRDS("/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS")
cat("class:", class(obj), "\n")
cat("dim:", paste(dim(obj), collapse=" x "), "\n")
cn <- colnames(obj)
cat("colnames:\n")
print(cn)

# find likely project/study id column and repo column
if ("project_name" %in% cn) {
  pn <- sort(unique(obj$project_name))
  cat("\nnum unique project_name:", length(pn), "\n")
  write.table(pn, "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/cross-repo-duplicate-check-e81250/work/final_atlas_project_names.txt",
              row.names=FALSE, col.names=FALSE, quote=FALSE)
}

# check for any obvious repo/source column
for (cand in c("repo","repository","source","dataset","source_repo","cohort")) {
  if (cand %in% cn) {
    cat("\nfound candidate col:", cand, "\n")
    print(table(obj[[cand]]))
  }
}

targets <- c("Karaayvas2018_Breast","KARAAYVAZ_M_ET_AL_2018_P20","Wu2021_Breast","WU_SZ_ET_AL_2021_P15",
             "Griffiths2021_Breast","GRIFFITHS_JI_ET_AL_2021_P18","Gao2021_Breast","GAO_R_ET_AL_2021_P17",
             "Savas2018_Breast","SAVAS_P_ET_AL_2018_P35")
if ("project_name" %in% cn) {
  cat("\n--- grep-like match of target names against project_name unique values ---\n")
  for (t in targets) {
    hits <- pn[grepl(t, pn, ignore.case=TRUE)]
    cat(t, " -> ", paste(hits, collapse="; "), "\n")
  }
  cat("\n--- broader substring search (author name only) ---\n")
  authors <- c("Karaayva","Wu2021","WU_SZ","Griffiths","GAO_R","Gao2021","Savas")
  for (a in authors) {
    hits <- pn[grepl(a, pn, ignore.case=TRUE)]
    cat(a, " -> ", paste(hits, collapse="; "), "\n")
  }
}
