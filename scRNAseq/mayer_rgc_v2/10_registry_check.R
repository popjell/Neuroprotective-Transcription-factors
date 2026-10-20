# Registry check: every test has a tier, Tier 2 padj is BH across exactly the five Tier 2 tests, counts by tier. Writes tests_registry_with_padj.csv (the registry itself is append-only).
source("common_v2.R")

reg <- read.csv(file.path(out_root, "tests_registry.csv"), stringsAsFactors = FALSE)
stopifnot(!anyNA(reg$tier), all(reg$tier %in% 1:3), !anyDuplicated(reg$test_id))
t2 <- reg$tier == 2
cat("Tier 2 tests:", paste(reg$test_id[t2], collapse = ", "), "\n")
stopifnot(sum(t2) == 5, setequal(reg$test_id[t2], c("B_ITM2C", "D_PSAP_interaction", "C_PSAP", "D_PSAP_resilient", "B_roast_set31")))
reg$padj_BH_tier2 <- NA_real_
reg$padj_BH_tier2[t2] <- p.adjust(reg$p[t2], method = "BH")
write.csv(reg, file.path(out_root, "tests_registry_with_padj.csv"), row.names = FALSE)
print(reg[reg$tier < 3, c("test_id", "tier", "estimate", "se", "p", "padj_BH_tier2")], digits = 4)
cnt <- table(reg$tier)
log_line("registry: ", nrow(reg), " tests (Tier 1: ", cnt["1"], ", Tier 2: ", cnt["2"], ", Tier 3: ", cnt["3"], "); all rows have a tier; Tier 2 padj is BH over exactly ", sum(t2), " tests")
