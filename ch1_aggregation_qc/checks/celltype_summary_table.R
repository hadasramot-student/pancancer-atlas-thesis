library(dplyr)
library(tidyr)
library(purrr)

# assuming all_cells is already loaded as a data.frame / tibble

# 1. total number of Breast projects
total_breast_projects <- all_cells %>%
  filter(cancer_type == "Breast") %>%
  distinct(project_name) %>%
  nrow()

# 2. function to calculate stats for each cell type
get_celltype_stats <- function(cell_type_name) {
  
  breast_df <- all_cells %>%
    filter(cancer_type == "Breast")
  
  tibble(
    Cell_Type = cell_type_name,
    
    Total_Breast_Projects = total_breast_projects,
    
    Consensus_projects = breast_df %>%
      filter(Consensus_Cell_Type == cell_type_name) %>%
      distinct(project_name) %>%
      nrow(),
    
    Scimilarity_projects = breast_df %>%
      filter(Meta.scimilarity == cell_type_name) %>%
      distinct(project_name) %>%
      nrow(),
    
    SingleR_projects = breast_df %>%
      filter(Meta.SingleR == cell_type_name) %>%
      distinct(project_name) %>%
      nrow(),
    
    Sctype_projects = breast_df %>%
      filter(Meta.Sctype == cell_type_name) %>%
      distinct(project_name) %>%
      nrow()
  )
}

# 3. get all unique cell types from all annotation methods
all_cell_types <- unique(c(
  all_cells$Consensus_Cell_Type,
  all_cells$Meta.scimilarity,
  all_cells$Meta.SingleR,
  all_cells$Meta.Sctype
))

all_cell_types <- all_cell_types[!is.na(all_cell_types)]

# 4. create final summary table
celltype_summary_table <- map_dfr(all_cell_types, get_celltype_stats)

# 5. view result
View(celltype_summary_table)

# optional: save to csv
write.csv(
  celltype_summary_table,
  "breast_celltype_project_summary.csv",
  row.names = FALSE
)