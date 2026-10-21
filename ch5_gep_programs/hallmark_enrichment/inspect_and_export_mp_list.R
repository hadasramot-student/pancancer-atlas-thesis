userlib <- "/home/hadasra/R/x86_64-pc-linux-gnu-library/4.5"
.libPaths(c(userlib, .libPaths()))

mp <- readRDS("/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/gavish_2023/tiroshlab-3ca-b124fc7/ITH_hallmarks/MPs_distribution/MP_list.RDS")

cat("class:", class(mp), "\n")
cat("length:", length(mp), "\n")
cat("names:", paste(names(mp), collapse=" | "), "\n\n")

str(mp, max.level = 2, list.len = 60)

# Try to flatten into a long-format data frame: cell_type, meta_program, gene
rows <- list()
for (ct in names(mp)) {
  inner <- mp[[ct]]
  cat("\n--", ct, "--\n")
  cat("  class:", class(inner), "length:", length(inner), "\n")
  if (is.list(inner)) {
    cat("  names:", paste(names(inner), collapse=" | "), "\n")
    for (mpname in names(inner)) {
      genes <- inner[[mpname]]
      rows[[length(rows)+1]] <- data.frame(cell_type=ct, meta_program=mpname, gene=genes, stringsAsFactors = FALSE)
    }
  }
}

if (length(rows) > 0) {
  df <- do.call(rbind, rows)
  write.csv(df, "/mnt/new_groups/ofircohen-group/users/hadasra/temp_objects/gavish_2023/gavish_MP_gene_lists_long.csv", row.names = FALSE)
  cat("\nWrote", nrow(df), "rows ->", "gavish_MP_gene_lists_long.csv\n")
  cat("Unique cell_types:", paste(unique(df$cell_type), collapse=", "), "\n")
  cat("Unique meta_programs:", length(unique(paste(df$cell_type, df$meta_program))), "\n")
}
