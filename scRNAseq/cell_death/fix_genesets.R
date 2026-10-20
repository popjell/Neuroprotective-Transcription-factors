# Writes data/yu_tier2_primary_fixed.gmt and data/yu_tier2_primary_fixes.csv. The original gmt is left as it is.
# Symbols follow the Trobisch gene annotation (older HGNC names, e.g. SEPT4, MARCH7, HIST1H2AA), which is what the analysis matches against.
# Only unambiguous fixes are made: OCR typos, current -> annotation symbol, one-gene aliases. Everything else that is not in the data is dropped.
suppressMessages(library(Matrix))
source("common.R")

d <- readRDS("data/cortex_neurons.rds")
genes_in_data <- rownames(d$counts)
sets <- read_gmt("data/yu_tier2_primary.gmt")

fix <- c(# OCR typos and transposed letters
         PTPMTI = "PTPMT1", PTTGIIP = "PTTG1IP", RACKI = "RACK1", RBICC1 = "RB1CC1", RHOTI = "RHOT1", RRPS = "RRP8", SENPI = "SENP1",
         SIAHl = "SIAH1", SIRTI = "SIRT1", SNAII = "SNAI1", SNWI = "SNW1", SYVNI = "SYVN1", TFDPi = "TFDP1", UBOLN1 = "UBQLN1",
         SsT = "SST", NFKB1A = "NFKBIA", IFR1 = "IRF1", IFR2 = "IRF2",
         # IFNA genes with the "NA" lost
         setNames(paste0("IFNA", c(1, 2, 4, 5, 6, 7, 8, 10, 13, 14, 16, 17, 21)), paste0("IF", c(1, 2, 4, 5, 6, 7, 8, 10, 13, 14, 16, 17, 21))),
         # current HGNC symbol -> symbol in the Trobisch annotation
         SEPTIN4 = "SEPT4", MARCHF7 = "MARCH7", MARCHF1 = "MARCH1", QARS1 = "QARS", GBA1 = "GBA", STING1 = "TMEM173", ELAPOR1 = "KIAA1324",
         COA8 = "APOPT1", SHFL = "C19orf66", H2AX = "H2AFX", H2AZ1 = "H2AFZ", H2AZ2 = "H2AFV", MACROH2A1 = "H2AFY", MACROH2A2 = "H2AFY2",
         H2AJ = "H2AFJ", H2AB1 = "H2AFB1", H2AB2 = "H2AFB2", H2AB3 = "H2AFB3",
         setNames(c(paste0("HIST1H2A", c("A", "B", "C", "D", "E", "G", "H", "I", "J", "K", "L", "M")), "HIST2H2AA3", "HIST2H2AA4", "HIST2H2AC", "HIST2H2AB", "HIST3H2A"),
                  paste0("H2AC", c(1, 4, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25))),
         # one-gene aliases and protein names
         TNFSF6 = "FASLG", TNFRSF6 = "FAS", NOX2 = "CYBB", TRIF = "TICAM1", TIRP = "TICAM2", PARP = "PARP1", HSP90A = "HSP90AA1", GLNA = "GLUL",
         NRF2 = "NFE2L2", PI3KC3 = "PIK3C3", VPS24 = "CHMP3", NCX1 = "SLC8A1")

log <- list()
fixed <- lapply(names(sets), function(s) {
  out <- character()
  for (g in unique(sets[[s]])) {
    if (g %in% genes_in_data) { out <- c(out, g); next }
    new <- if (g %in% names(fix)) fix[[g]] else NA_character_
    action <- if (is.na(new)) "dropped" else if (!new %in% genes_in_data) "dropped" else if (new %in% sets[[s]] || new %in% out) "dropped (replacement already in set)" else "renamed"
    reason <- if (is.na(new)) {
      if (grepl("^MIR|^TRC|^TRD|^TRG|^ENSG|^NONHSAG|-AS1$|-DT$|^SNHG|^LINC", g)) "non-coding RNA, not in the data"
      else "not in the data (family or protein name, pseudogene or unresolved symbol)"
    } else if (!new %in% genes_in_data) paste("replacement", new, "is not in the data") else ""
    log[[length(log) + 1]] <<- data.frame(set = s, original = g, action = action, new = if (action == "renamed") new else "", reason = reason)
    if (action == "renamed") out <- c(out, new)
  }
  out
})
names(fixed) <- names(sets)
writeLines(vapply(names(fixed), function(s) paste(c(s, "Yu tier2 primary, symbols fixed to the Trobisch annotation", fixed[[s]]), collapse = "\t"), ""), "data/yu_tier2_primary_fixed.gmt")
write.csv(do.call(rbind, log), "data/yu_tier2_primary_fixes.csv", row.names = FALSE)

# Before and after: genes in the data and subtypes where the set is testable (at least detect_genes genes detected in detect_frac of cells)
testable_n <- function(x) sum(sapply(d$subtypes, function(s) sum(d$frac[intersect(x, rownames(d$frac)), s] >= detect_frac) >= detect_genes))
summary <- data.frame(set = names(sets), in_gmt = lengths(lapply(sets, unique)),
                      data_before = sapply(sets, function(x) sum(unique(x) %in% genes_in_data)), data_after = lengths(fixed),
                      testable_before = sapply(sets, function(x) testable_n(unique(x))), testable_after = sapply(fixed, testable_n),
                      row.names = NULL)
summary$pct_in_data_after <- round(100 * summary$data_after / summary$in_gmt, 1)
print(summary)
cat("\nactions:\n"); print(table(do.call(rbind, log)$action))
stopifnot(all(unlist(fixed) %in% genes_in_data), !anyDuplicated(unlist(lapply(fixed, function(x) x[duplicated(x)]))))
