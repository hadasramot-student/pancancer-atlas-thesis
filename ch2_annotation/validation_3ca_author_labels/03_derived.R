# 03_derived.R -- derived views from the 9x9 confusion matrix (no re-load of atlas)
options(width = 220)
WT <- "/ofir_grp/users/hadasra/create_sc_object_workflow/.claude/worktrees/3ca-author-label-validation-35ef05/work"
cm <- as.matrix(read.csv(file.path(WT,"confusion_matrix_9x9.csv"), row.names=1, check.names=FALSE))
lv <- colnames(cm)
N  <- sum(cm)
say <- function(...) cat(..., "\n")

say("total validated cells:", N)
say("overall 9-class accuracy:", round(sum(diag(cm))/N, 4))

# row-normalised (what does the consensus call become, per author-truth? -> col-normalised)
col_norm <- sweep(cm, 2, colSums(cm), "/")
row_norm <- sweep(cm, 1, rowSums(cm), "/")
say("\n=== COLUMN-normalised (each col = one 3CA author meta-type; entries = fraction of those cells given each consensus call) ===")
print(round(col_norm, 3))
say("\n=== ROW-normalised (each row = one consensus call; entries = fraction that were each 3CA author meta-type) ===")
print(round(row_norm, 3))

# treat consensus 'Other/Unknown' as an abstention, not a contradiction
other <- "Other/Unknown"
diag_sum <- sum(diag(cm))
abstain  <- sum(cm[other, ]) - cm[other, other]        # consensus said Other, author said something specific
say("\ndiagonal (exact agreement):", diag_sum, sprintf("(%.1f%%)", 100*diag_sum/N))
say("consensus abstained (=Other/Unknown) while author gave a specific lineage:", abstain,
    sprintf("(%.1f%%)", 100*abstain/N))
say("=> exact-agreement OR consensus-abstention ('not contradicted'):", diag_sum+abstain,
    sprintf("(%.1f%%)", 100*(diag_sum+abstain)/N))

# restrict to cells where consensus made a positive (non-Other) call
pos <- setdiff(lv, other)
cm_pos <- cm[pos, , drop=FALSE]
npos <- sum(cm_pos)
agree_pos <- sum(diag(cm_pos[, pos, drop=FALSE]))
say("\ncells where consensus made a positive lineage call:", npos)
say("  ...of those, author agreed (same lineage):", agree_pos, sprintf("(%.1f%%)", 100*agree_pos/npos))

# 8-class, both sides positive
cm8 <- cm[pos, pos, drop=FALSE]
say("\n8-class (both sides give a specific lineage): n =", sum(cm8),
    " accuracy =", round(sum(diag(cm8))/sum(cm8), 4))
p8 <- diag(cm8)/rowSums(cm8); r8 <- diag(cm8)/colSums(cm8); f8 <- 2*p8*r8/(p8+r8)
print(data.frame(meta=pos, precision=round(p8,3), recall=round(r8,3), f1=round(f8,3)))
say("8-class macro-F1:", round(mean(f8, na.rm=TRUE), 4))
say("8-class weighted-F1 (by author n):", round(sum(f8*colSums(cm8), na.rm=TRUE)/sum(colSums(cm8)), 4))

# biggest off-diagonal flows
od <- cm; diag(od) <- 0
idx <- which(od > 0, arr.ind=TRUE)
fl <- data.frame(consensus=rownames(od)[idx[,1]], author_3CA=colnames(od)[idx[,2]], n=od[idx])
fl <- fl[order(-fl$n), ][1:15, ]
say("\n=== 15 largest disagreement flows (consensus -> author truth) ===")
print(fl, row.names=FALSE)
