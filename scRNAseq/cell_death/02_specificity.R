# Is GENE associated with each program more strongly than expression-matched genes are? One test per program, plus a donor bootstrap.
# Per donor x subtype: partial Spearman r of a gene with the program score (adjusting rank nUMI and rank genes detected), combined across donors by random effects.
# Statistic per program = median across subtypes. GENE and the matched genes go through exactly the same steps.
source("common.R")

out <- res("program_specificity.csv")
no_overwrite(out)

d <- readRDS("data/cortex_neurons.rds")
x <- d$counts
cells <- d$cells
scores <- readRDS(res("program_scores.rds"))
testable <- read.csv(res("testability.csv"))
programs <- colnames(scores)
tsub <- sort(unique(testable$subtype[testable$testable]))
donor_ids <- sort(unique(cells$donor))
stopifnot(gene %in% rownames(x))
row_medians <- function(m) matrixStats::rowMedians(m, na.rm = TRUE)

# Log-normalised expression, mean expression and fraction detected, per subtype
lognorm <- list()
gene_mean <- gene_frac <- matrix(NA_real_, nrow(x), length(tsub), dimnames = list(rownames(x), tsub))
for (s in tsub) {
  i <- cells$subtype == s
  l <- x[, i]
  l@x <- log1p(l@x / rep(cells$nUMI[i], diff(l@p)) * norm_scale)
  lognorm[[s]] <- l
  gene_mean[, s] <- Matrix::rowSums(l) / ncol(l)
  gene_frac[, s] <- Matrix::rowSums(l > 0) / ncol(l)
}

# One set of matched genes for all subtypes: closest to GENE in standardised (mean expression, fraction detected), averaged over subtypes
exclude <- unique(c(gene, readLines("data/exclude_genes.txt"), unlist(read_gmt(genesets))))
pool <- setdiff(rownames(x), exclude)
pool <- pool[apply(gene_frac[pool, , drop = FALSE] > 0, 1, all)]
distance <- sapply(tsub, function(s) {
  sqrt(((gene_mean[pool, s] - gene_mean[gene, s]) / sd(gene_mean[pool, s]))^2 + ((gene_frac[pool, s] - gene_frac[gene, s]) / sd(gene_frac[pool, s]))^2)
})
null_genes <- pool[order(rowMeans(distance))[seq_len(n_null)]]
writeLines(null_genes, res("null_genes.txt"))
cat("pool:", length(pool), " matched genes:", n_null, " max mean distance:", round(max(rowMeans(distance)[null_genes]), 3), "\n")

# Donor-level Fisher z of every gene (GENE first) with every testable program score, per subtype
genes <- c(gene, null_genes)
fisher <- var_z <- list()
for (s in tsub) {
  i <- which(cells$subtype == s)
  progs <- testable$program[testable$subtype == s & testable$testable]
  m <- as.matrix(t(lognorm[[s]][genes, ]))
  sc <- scores[i, progs, drop = FALSE]
  don <- cells$donor[i]
  donors_s <- sort(unique(don))
  z <- array(NA_real_, c(length(donors_s), length(genes), length(progs)), dimnames = list(donors_s, genes, progs))
  v <- setNames(numeric(length(donors_s)), donors_s)
  for (k in seq_along(donors_s)) {
    j <- which(don == donors_s[k])
    v[k] <- 1 / (length(j) - 3)
    covariates <- qr(cbind(1, rank(cells$nUMI[i][j]), rank(cells$ngen[i][j])))
    z[k, , ] <- clipz(suppressWarnings(cor(qr.resid(covariates, apply(m[j, ], 2, rank)), qr.resid(covariates, apply(sc[j, , drop = FALSE], 2, rank)))))
  }
  fisher[[s]] <- z
  var_z[[s]] <- v
  cat(s, ":", length(donors_s), "donors,", length(progs), "programs\n")
}

# Median across subtypes of the combined r (rows = genes, columns = programs). donor_count: how often each donor is drawn (NULL = every donor once)
program_stat <- function(donor_count = NULL) {
  r <- array(NA_real_, c(length(genes), length(tsub), length(programs)), dimnames = list(genes, tsub, programs))
  for (s in tsub) {
    pick <- seq_along(var_z[[s]])
    if (!is.null(donor_count)) pick <- rep(pick, as.integer(donor_count[names(var_z[[s]])]))
    if (length(unique(pick)) < 2) next
    for (p in dimnames(fisher[[s]])[[3]]) r[, s, p] <- dl_r(fisher[[s]][pick, , p], var_z[[s]][pick])
  }
  list(stat = vapply(programs, function(p) row_medians(r[, , p]), numeric(length(genes))), r = r)
}
compare <- function(stat, p) {
  r_gene <- stat[1, p]
  r_null <- stat[-1, p][is.finite(stat[-1, p])]
  if (!is.finite(r_gene) || length(r_null) < 100) return(NULL)
  list(r = r_gene, n_null = length(r_null), null_mean = mean(r_null), null_sd = sd(r_null), null_median = median(r_null),
       z = (r_gene - mean(r_null)) / sd(r_null), p = (1 + sum(abs(r_null) >= abs(r_gene))) / (length(r_null) + 1))
}

full <- program_stat()
tab <- do.call(rbind, lapply(programs, function(p) {
  a <- compare(full$stat, p)
  if (is.null(a)) return(NULL)
  data.frame(program = p, n_subtypes = sum(is.finite(full$r[1, , p])), n_subtypes_positive = sum(full$r[1, , p] > 0, na.rm = TRUE),
             r_median = a$r, n_null = a$n_null, null_mean = a$null_mean, null_sd = a$null_sd, null_median = a$null_median,
             z = a$z, excess_r = a$r - a$null_median, p_two_sided = a$p)
}))
tab$q_BH <- p.adjust(tab$p_two_sided, "BH")
cat("programs:", nrow(tab), " BH q < 0.05:", sum(tab$q_BH < 0.05), "\n")

# Donor bootstrap: donors resampled with replacement once per replicate, shared by all subtypes (a subtype needs at least 2 distinct donors)
boot_z <- boot_r <- matrix(NA_real_, n_boot, length(programs), dimnames = list(NULL, programs))
for (b in seq_len(n_boot)) {
  drawn <- table(factor(sample(donor_ids, length(donor_ids), replace = TRUE), levels = donor_ids))
  stat <- program_stat(drawn)$stat
  for (p in programs) {
    a <- compare(stat, p)
    if (!is.null(a)) {
      boot_z[b, p] <- a$z
      boot_r[b, p] <- a$r
    }
  }
  if (b %% 100 == 0) cat("bootstrap", b, "of", n_boot, "\n")
}
write.csv(data.frame(replicate = seq_len(n_boot), boot_z), res("bootstrap_z.csv"), row.names = FALSE)
boot_q <- function(m, p, q) unname(quantile(m[, p], q, na.rm = TRUE))
tab$z_boot_lo <- sapply(tab$program, boot_q, m = boot_z, q = 0.025)
tab$z_boot_hi <- sapply(tab$program, boot_q, m = boot_z, q = 0.975)
tab$frac_boot_z_ge_1.96 <- sapply(tab$program, function(p) mean(boot_z[, p] >= 1.96, na.rm = TRUE))
tab$r_boot_lo <- sapply(tab$program, boot_q, m = boot_r, q = 0.025)
tab$r_boot_hi <- sapply(tab$program, boot_q, m = boot_r, q = 0.975)
tab <- tab[order(-tab$z), ]
write.csv(tab, out, row.names = FALSE)
writeLines(capture.output(sessionInfo()), res("sessionInfo.txt"))
