# Aggregate Consensus_Cell_Type x cancer_type counts for Ch.2 Fig 1 (no change to cell selection:
# full gapfix atlas minus the 5 cross-repo duplicate Breastscape projects = 8,107,665 cells).
suppressPackageStartupMessages(library(dplyr))
md <- readRDS("/ofir_grp/users/hadasra/temp_objects/all_cells.gapfix_20260727_164657.RDS")
dups <- c("WU_SZ_ET_AL_2021_P15_filtered", "GRIFFITHS_JI_ET_AL_2021_P18_filtered",
          "GONZALEZ_H_ET_AL_2022_P31_filtered", "SAVAS_P_ET_AL_2018_P35_filtered",
          "GAO_R_ET_AL_2021_P17_filtered")
md <- md[!md$project_name %in% dups, c("project_name", "cancer_type", "Consensus_Cell_Type")]
cat("cells:", nrow(md), " projects:", length(unique(md$project_name)), "\n")
cat("NA cancer_type:", sum(is.na(md$cancer_type)), " NA Consensus:", sum(is.na(md$Consensus_Cell_Type)), "\n")
print(table(md$Consensus_Cell_Type, useNA = "ifany"))
out <- md %>% count(cancer_type, Consensus_Cell_Type, name = "n")
write.csv(out, "/ofir_grp/users/hadasra/create_sc_object_workflow/plots_and_outputs/ch2_fig1_print/counts_long.csv", row.names = FALSE)
