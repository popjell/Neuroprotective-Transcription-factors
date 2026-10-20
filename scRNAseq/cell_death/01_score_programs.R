# UCell score of every gene set in GENESETS in every analysed cell, and which sets are testable in which subtype
source("common.R")
library(UCell)

out <- res("program_scores.rds")
no_overwrite(out)

d <- readRDS("data/cortex_neurons.rds")
sets <- lapply(read_gmt(genesets), function(g) setdiff(intersect(g, rownames(d$counts)), gene))
sets <- sets[lengths(sets) > 0]

testable <- do.call(rbind, lapply(names(sets), function(p) {
  do.call(rbind, lapply(d$subtypes, function(s) {
    n <- sum(d$frac[sets[[p]], s] >= detect_frac)
    data.frame(program = p, subtype = s, n_genes = length(sets[[p]]), n_detected = n, testable = n >= detect_genes)
  }))
}))
write.csv(testable, res("testability.csv"), row.names = FALSE)
cat("sets:", length(sets), " testable set x subtype pairs:", sum(testable$testable), "of", nrow(testable), "\n")

ids <- setNames(sets, paste0("S", seq_along(sets)))
scores <- ScoreSignatures_UCell(d$counts, features = ids, maxRank = max_rank, ncores = 1, name = "")
scores <- scores[, names(ids), drop = FALSE]
colnames(scores) <- names(sets)
saveRDS(scores, out)
