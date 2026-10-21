# Read-only audit: per-project and per-(project, cancer_type) counts on the deduplicated atlas.
md <- readRDS("/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS")
dups <- c("WU_SZ_ET_AL_2021_P15_filtered", "GRIFFITHS_JI_ET_AL_2021_P18_filtered",
          "GONZALEZ_H_ET_AL_2022_P31_filtered", "SAVAS_P_ET_AL_2018_P35_filtered",
          "GAO_R_ET_AL_2021_P17_filtered")
cat("all projects before dedup:", length(unique(md$project_name)), "\n")
md <- md[!md$project_name %in% dups, ]
cat("cells:", nrow(md), " projects:", length(unique(md$project_name)), "\n")
cat("columns:\n"); print(names(md))
keep <- intersect(c("project_name", "cancer_type", "SRX_accession", "sample", "orig.ident",
                    "sample_name", "study", "dataset", "source", "repository", "origin", "GSE",
                    "SRP", "study_id", "cell_name"), names(md))
pc <- as.data.frame(table(project_name = md$project_name, cancer_type = md$cancer_type))
pc <- pc[pc$Freq > 0, ]
names(pc)[3] <- "n_cells"
pc$has_SRX <- tapply(!is.na(md$SRX_accession), md$project_name, any)[as.character(pc$project_name)]
pc$has_cell_name <- tapply(!is.na(md$cell_name), md$project_name, any)[as.character(pc$project_name)]
write.csv(pc, "/ofir_grp/users/hadasra/create_sc_object_workflow/plots_and_outputs/project_count_audit/project_cancer_pairs_raw.csv", row.names = FALSE)
