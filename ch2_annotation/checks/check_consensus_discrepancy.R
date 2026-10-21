df <- readRDS("/ofir_grp/users/hadasra/thesis/sources/all_cells.RDS")

cat("=== class ===\n")
print(class(df))
cat("\n=== dim ===\n")
print(dim(df))
cat("\n=== colnames ===\n")
print(colnames(df))

needed <- c("Consensus_Cell_Type", "Meta.SingleR", "Meta.Sctype", "Meta.scimilarity",
            "maj2_cell", "maj3_cell", "maj2_cluster", "maj3_cluster")
cat("\n=== which needed cols exist ===\n")
print(needed[needed %in% colnames(df)])
print(needed[!needed %in% colnames(df)])

if (!"Consensus_Cell_Type" %in% colnames(df)) {
  cat("\nConsensus_Cell_Type column not found -- stopping structural report only.\n")
  quit(save = "no")
}

total_n <- nrow(df)
cat("\n=== TOTAL ROWS ===\n")
print(total_n)

cons <- df$Consensus_Cell_Type
is_na_cons <- is.na(cons)
cat("\n=== Consensus_Cell_Type NA vs non-NA ===\n")
cat("non-NA:", sum(!is_na_cons), "\n")
cat("NA:", sum(is_na_cons), "\n")

cat("\n=== Consensus_Cell_Type value counts (top 20) ===\n")
print(head(sort(table(cons, useNA = "ifany"), decreasing = TRUE), 20))

is_other <- !is_na_cons & cons == "other"
cat("\n=== 'other' (Other/Unknown) count ===\n")
cat("n other:", sum(is_other), "\n")

m1 <- df$Meta.SingleR
m2 <- df$Meta.Sctype
m3 <- df$Meta.scimilarity

na_count_per_cell <- is.na(m1) + is.na(m2) + is.na(m3)

cat("\n=== Among 'other' cells: how many of the 3 method columns (Meta.SingleR/Meta.Sctype/Meta.scimilarity) are NA ===\n")
print(table(na_count_per_cell[is_other]))

cat("\n=== Among NA-consensus cells: how many of the 3 method columns are NA ===\n")
print(table(na_count_per_cell[is_na_cons]))

cat("\n=== Overlap check: is_na_cons vs is_other ===\n")
cat("n is_na_cons:", sum(is_na_cons), "\n")
cat("n is_other:", sum(is_other), "\n")
cat("n both TRUE (overlap):", sum(is_na_cons & is_other), "\n")
cat("n neither:", sum(!is_na_cons & !is_other), "\n")
cat("n exactly one of the two:", sum(xor(is_na_cons, is_other)), "\n")

cat("\n=== For 'other' cells, breakdown of pairwise agreement (0,1,2,3 populated methods that agree) ===\n")
other_df <- data.frame(m1 = m1[is_other], m2 = m2[is_other], m3 = m3[is_other])
n_na <- rowSums(is.na(other_df))
agree_count <- mapply(function(a,b,c) {
  vals <- c(a,b,c)
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(NA_integer_)
  tab <- table(vals)
  max(tab)
}, other_df$m1, other_df$m2, other_df$m3)

cat("n_na distribution among 'other' cells (0/1/2/3 methods missing):\n")
print(table(n_na))

cat("\nmax-agreement-count distribution among 'other' cells (i.e. size of largest matching subgroup among populated labels):\n")
print(table(agree_count, useNA = "ifany"))

cat("\n=== Cross classification: n_na (missing methods) x agree_count, for 'other' cells ===\n")
print(table(n_na, agree_count, useNA = "ifany"))

cat("\nDONE\n")
