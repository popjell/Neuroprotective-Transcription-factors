library(tidyr)
library(dplyr)
suppressMessages({
    library(ggplot2)
    library(gridExtra)
    library(RColorBrewer)
    library(ggrepel)
})

#ask if we analyse TFs or clusters
use_clusters <- tolower(trimws(readline("Cluster the TFs? (y/n): "))) %in% c("y", "yes")

base_dir <- "filtered_genelists"
iregulon_results_dir <- "filtered_genelists/graph_output"
iregulon_graphs_dir <- "filtered_genelists/iregulon_graphs"

target_table <- read.csv("reformatted_edges.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
U_nodes <- read.csv("intersected_upregulated_genes.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(U_nodes) <- "gene"
full_regulators <- read.table("regulators.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
colnames(full_regulators) <- "Regulator"
expression_data <- read.table("Expression_data.txt", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

#cliffs delta: P(target > rest) - P(target < rest), +1..-1, 0 = no separation
cliffs_delta <- function(a, b) {
    n <- length(a); m <- length(b)
    (2 * (sum(rank(c(a, b))[seq_len(n)]) - n * (n + 1) / 2) - n * m) / (n * m)
}

#symmetric diverging colour scale for cliff's delta: grey at 0, colour stops at
#the quartiles of the absolute values so +x and -x get the same intensity
make_cliff_scale <- function(v) {
    v <- v[!is.na(v)]
    L <- max(abs(v))
    absv <- abs(v[v != 0])
    qs <- if (length(absv) > 0) as.numeric(quantile(absv, c(0.25, 0.5, 0.75))) else numeric(0)
    stops <- unique(sort(c(-L, -rev(qs), 0, qs, L)))
    positions <- (stops + L) / (2 * L)
    cols <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(stops))
    scale_color_gradientn(
        colors = cols, values = positions, limits = c(-L, L),
        oob = scales::squish, na.value = "grey85",
        name = "Cliff's delta\n(n>=3, shrunk)",
        breaks = signif(pretty(c(-L, L)), 2)
    )
}

#build units (TFs or clusters) with their member TFs and target genes
if (use_clusters) {
    if (!exists("cl") || !exists("dominant_tf_per_cluster")) {
        stop("TF clusters not found in current scope. Run script 10 first.")
    }
    unit_tf_list <- setNames(
        lapply(dominant_tf_per_cluster$cluster, function(cs) intersect(names(cl[cl == cs]), full_regulators$Regulator)),
        dominant_tf_per_cluster$label
    )
} else {
    unit_tf_list <- setNames(as.list(full_regulators$Regulator), full_regulators$Regulator)
}
unit_names <- names(unit_tf_list)
unit_targets <- lapply(unit_tf_list, function(tfs) unique(target_table$`Target Gene`[target_table$`Regulator Gene` %in% tfs]))
unit_label <- if (use_clusters) "Cluster" else "TF"

all_coverage_tables <- list()
p1_list <- list(); p2_list <- list(); p3_list <- list()
for(csv in list.files(iregulon_results_dir, pattern = "\\.csv$", full.names = TRUE)) {
    name <- sub("\\_graph_output.csv$", "", basename(csv))
    print(name)
    regulators_in_G_ireg <- read.table(file.path(paste0(iregulon_graphs_dir, "/", name), "regulators.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    colnames(regulators_in_G_ireg) <- c("Regulator")
    not <- regulators_in_G_ireg %>%
        dplyr::filter(!Regulator %in% full_regulators$Regulator)

    print(paste0(name, " amount: ", nrow(regulators_in_G_ireg), " amount not in full list: ", nrow(not)))
    
    G_nodes = read.table(paste0(base_dir, "/", name, "_genes.txt"), header = FALSE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    colnames(G_nodes) <- c("name")
    U_genes <- unique(U_nodes$gene)
    G_genes <- unique(G_nodes$name)

    unotg_nodes <- setdiff(U_genes, G_genes)

    total_U <- length(U_genes)
    total_G <- length(G_genes)
    total_UnotG <- length(unotg_nodes)
    
    UnotG_intersection <- intersect(U_nodes$gene, G_nodes$name)
    #make master table with all calculations i need
    genelist_target_table <- data.frame(
        `Regulator Gene` = unit_names,
        Targets = sapply(unit_targets, paste, collapse = ", "),
        stringsAsFactors = FALSE,
        check.names = FALSE
    ) %>% 
        dplyr::mutate(
            num_targets = lengths(strsplit(Targets, ",\\s*")),
            num_targets_in_G = sapply(strsplit(Targets, ",\\s*"), function(x) sum(x %in% G_nodes$name)),
            targets_in_G = sapply(strsplit(Targets, ",\\s*"), function(x) paste(intersect(x, G_nodes$name), collapse = ", ")),
            targets_in_UnotG = sapply(strsplit(Targets, ",\\s*"), function(x) paste(intersect(x, unotg_nodes), collapse = ", ")),
            num_targets_in_UnotG = sapply(strsplit(Targets, ",\\s*"), function(x) sum(x %in% unotg_nodes)),
            U_coverage = num_targets / total_U,
            G_coverage = num_targets_in_G / total_G,
            specificity = G_coverage / U_coverage,
            log2spec = {
              val = log2(specificity)
              ifelse(is.infinite(val) | is.na(val) | val==0, NA, val)
            }
        )
    
    
    #FET + Wilcoxon on Exp data
    genelist_target_table$fisher_p_value <- NA_real_
    genelist_target_table$fisher_odds_ratio <- NA_real_ 
    genelist_target_table$wilcox_p_value <- NA_real_
    genelist_target_table$upregulated_in_G <- NA 
    genelist_target_table$degree_of_upregulation <- NA_real_
    genelist_target_table$cliff_delta <- NA_real_
    genelist_target_table$n_targets_expr <- NA_real_
  
    targets_expression_data <- expression_data %>%
        dplyr::select(symbol, log2FoldChange_MN, log2FoldChange_RGC, pvalue_MN, pvalue_RGC, EAE_Combined_Score)

    
    for(TF in unit_names){
        #FET
        #make 2x2 contingency table
        UnotG_target_list <- genelist_target_table %>%
            dplyr::filter(`Regulator Gene` == TF) %>%
            dplyr::pull(targets_in_UnotG) %>%
            strsplit(",\\s*") %>%
            unlist() %>%
            unique()
      
        G_target_list <- genelist_target_table %>%
            dplyr::filter(`Regulator Gene` == TF) %>%
            dplyr::pull(targets_in_G) %>%
            strsplit(",\\s*") %>%
            unlist()%>%
            unique()
      
        # Remove empty strings from strsplit("")
        G_target_list <- G_target_list[G_target_list != "" & !is.na(G_target_list)]
        if(length(G_target_list) == 0) {
          genelist_target_table <- genelist_target_table %>%
                dplyr::mutate(
                    fisher_p_value = ifelse(`Regulator Gene` == TF, NA_real_, fisher_p_value),
                    fisher_odds_ratio = ifelse(`Regulator Gene` == TF, NA_real_, fisher_odds_ratio),
                    wilcox_p_value = ifelse(`Regulator Gene` == TF, NA_real_, wilcox_p_value),
                    upregulated_in_G = ifelse(`Regulator Gene` == TF, NA, upregulated_in_G),
                    degree_of_upregulation = ifelse(`Regulator Gene` == TF, NA_real_, degree_of_upregulation)
                )
          next
        }
      
      
        num_targets_G <- length(G_target_list)
        num_targets_UnotG <- length(UnotG_target_list)
      
        
        num_non_targets_G <- total_G - num_targets_G
        num_non_targets_UnotG <- total_UnotG - num_targets_UnotG
      
        
        contingency_table <- matrix(c(num_targets_G, num_non_targets_G, num_targets_UnotG, num_non_targets_UnotG), nrow = 2, byrow = TRUE)
        fisher_result <- fisher.test(contingency_table)
        p_value <- fisher_result$p.value
        odds_ratio <- fisher_result$estimate
        genelist_target_table <- genelist_target_table %>%
            dplyr::mutate(
                fisher_p_value = ifelse(`Regulator Gene` == TF, p_value, fisher_p_value),
                fisher_odds_ratio = ifelse(`Regulator Gene` == TF, odds_ratio, fisher_odds_ratio)
            )
      
      
        #Wilcoxon test on expression data
        TF_list_expression <- targets_expression_data %>%
            dplyr::filter(symbol %in% G_target_list) %>%
            dplyr::pull(EAE_Combined_Score) %>%
            unlist()%>%
        na.omit()

        GnotTF_list_expression <- targets_expression_data %>%
            dplyr::filter(!symbol %in% G_target_list & symbol %in% G_genes) %>%
            dplyr::pull(EAE_Combined_Score) %>%
            unlist()%>%
            na.omit()

        n_expr_val <- length(TF_list_expression)
        expr_ok <- n_expr_val >= 3 && length(GnotTF_list_expression) >= 3

        w_p_value <- if (expr_ok) {
            wilcox.test(TF_list_expression, GnotTF_list_expression, alternative = "two.sided")$p.value
        } else {
            NA_real_
        }
        print(paste0("Wilcoxon p-value for ", TF, ": ", w_p_value))
        cliff_val <- if (expr_ok) cliffs_delta(TF_list_expression, GnotTF_list_expression) else NA_real_
        upreg_flag <- if (expr_ok && !is.na(w_p_value) && w_p_value < 0.05) {
            mean(TF_list_expression) > mean(GnotTF_list_expression)
        } else {
            NA
        }
        degree_val <- if (expr_ok) mean(TF_list_expression) - mean(GnotTF_list_expression) else NA_real_
        genelist_target_table <- genelist_target_table %>%
            dplyr::mutate(
                wilcox_p_value = ifelse(`Regulator Gene` == TF, w_p_value, wilcox_p_value),
                cliff_delta = ifelse(`Regulator Gene` == TF, cliff_val, cliff_delta),
                n_targets_expr = ifelse(`Regulator Gene` == TF, n_expr_val, n_targets_expr),
                upregulated_in_G = ifelse(`Regulator Gene` == TF, upreg_flag, upregulated_in_G),
                degree_of_upregulation = ifelse(`Regulator Gene` == TF, degree_val, degree_of_upregulation)
            )
      
    }
    assign(paste0(name, "_coverage_table"), genelist_target_table)
  
    
    
    #plots
    #region
  
    plot_df <- genelist_target_table %>%
    dplyr::mutate(
        TF = `Regulator Gene`,
        specificity = as.numeric(specificity),
        log2spec = as.numeric(log2spec),
        degree_of_upregulation = as.numeric(degree_of_upregulation),
        fisher_p_value = as.numeric(fisher_p_value),
        fisher_odds_ratio = as.numeric(fisher_odds_ratio),
        U_coverage_pct = U_coverage * 100,
        G_coverage_pct = G_coverage * 100,
        wilcox_p_value = as.numeric(wilcox_p_value),
        neglog10_wilcox = -log10(pmax(wilcox_p_value, .Machine$double.xmin)),
        neglog10_fisher = -log10(pmax(fisher_p_value, .Machine$double.xmin)),
        cliff_delta = as.numeric(cliff_delta),
        n_targets_expr = as.numeric(n_targets_expr),
        cliff_clean = cliff_delta * n_targets_expr / (n_targets_expr + 5),
        higher_coverage = factor(
            {
              ifelse(fisher_p_value < 0.05, G_coverage_pct > U_coverage_pct, NA)
            },
            levels = c(FALSE, TRUE, NA),
            labels = c("Not higher in G", "Higher in G", "NA"),
            exclude = NULL
        ),
        upregulated_in_G = factor(
            upregulated_in_G,
            levels = c(FALSE, TRUE, NA),
            labels = c("Not higher in G", "Higher in G", na = "NA"),
            exclude = NULL
        )
    )
    p1 <- ggplot(

    #p1 - scatter plot of U_coverage vs G_coverage, color = cleaned effect size
    plot_df,
    aes(
        x = G_coverage_pct,
        y = U_coverage_pct,
        color = cliff_clean
    )
    ) + geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        linewidth = 0.7
    ) + geom_point(
        size = 3.5,
        alpha = 0.85
    ) +
    geom_text_repel(
        aes(label = TF),
        size = 3,
        color = "grey30",
        max.overlaps = Inf,
        min.segment.length = 0,
        seed = 1
    ) +
    labs(
        title = paste0(name, ": ", unit_label, " Target Coverage"),
        x = paste0("% of Functional Group G Covered by ", unit_label),
        y = paste0("% of U Covered by ", unit_label)
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        axis.text = element_text(color = "black")
    )

    assign(paste0(name, "_p1"), p1)
    p1_list[[name]] <- p1
    #p2 - bar plot of ranking TF's by specificity
    p2 <- plot_df %>%
    dplyr::filter(is.finite(log2spec)) %>%
    mutate(
        TF = reorder(TF, log2spec)
    ) %>%
    ggplot(
        aes(
            x = TF,
            y = log2spec
        )
    ) +
    geom_col() +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.6
    ) +
    coord_flip() +
    labs(
        title = paste0(name, ": ", unit_label, " Functional Specificity"),
        x = unit_label,
        y = expression(log[2]~"Specificity")
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank()
    )

    assign(paste0(name, "_p2"), p2)
    p2_list[[name]] <- p2
    #p3 - x by y scatter plot where x axis is cliff's delta and y axis is specificity, static purple

    p3 <- ggplot(
    plot_df,
    aes(
        x = cliff_clean,
        y = log2spec
    )
) +
    geom_hline(
        yintercept = 0,
        linetype = "dashed",
        linewidth = 0.6
    ) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        linewidth = 0.6
    ) +
    geom_point(
        size = 4,
        alpha = 0.9,
        color = "purple"
    ) +
    geom_text_repel(
        aes(label = TF),
        size = 3,
        color = "purple",
        max.overlaps = Inf,
        min.segment.length = 0,
        seed = 1
    ) +
    labs(
        title = paste0(name, ": ", unit_label, " Specificity vs Effect"),
        x = "Cliff's delta (n>=3, shrunk)",
        y = expression(log[2]~"Specificity")
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank()
    )

    assign(paste0(name, "_p3"), p3)
    p3_list[[name]] <- p3
    #endregion
    all_coverage_tables[[name]] <- genelist_target_table
  
}

all_coverage <- bind_rows(
    all_coverage_tables,
    .id = "functional_group"
)

all_coverage <- all_coverage %>%
    mutate(
        TF = `Regulator Gene`,
        degree_of_upregulation = as.numeric(degree_of_upregulation),
        specificity = as.numeric(specificity),
        n_targets_expr = as.numeric(n_targets_expr),
        cliff_clean = cliff_delta * n_targets_expr / (n_targets_expr + 5)
    )

#consistent axis ranges across groups so panels are directly comparable;
#pad the ranges so repelled labels have room and points don't sit on the edge
pad_lims <- function(r, p = 0.15) {
    d <- diff(r)
    if (is.na(d) || d == 0) r + c(-p, p) else r + c(-d, d) * p
}
cov_xlim   <- pad_lims(range(all_coverage$G_coverage * 100))
cov_ylim   <- pad_lims(range(all_coverage$U_coverage * 100))
spec_ylim  <- pad_lims(range(all_coverage$log2spec, na.rm = TRUE))
cliff_xlim <- pad_lims(range(all_coverage$cliff_clean, na.rm = TRUE))

cliff_scale <- make_cliff_scale(all_coverage$cliff_clean)

for (nm in names(p1_list)) {
    p1_list[[nm]] <- p1_list[[nm]] + coord_cartesian(xlim = cov_xlim, ylim = cov_ylim) + cliff_scale
    p2_list[[nm]] <- p2_list[[nm]] + coord_flip(ylim = spec_ylim)
    p3_list[[nm]] <- p3_list[[nm]] + coord_cartesian(xlim = cliff_xlim, ylim = spec_ylim)
}

#p4 - bubble plot, size = specificity, colour = effect size, all groups
#size aesthetic is specificity^1.5; squish outliers above the 95th percentile
#of specificity to the cap so they don't compress the whole bulk into small sizes
size_range  <- if (use_clusters) c(2, 12) else c(0.2, 6.5)
size_breaks <- if (use_clusters) c(0.8, 1.1, 1.4, 1.5) else c(0.9, 1.1, 1.3, 1.4)
size_cap    <- as.numeric(quantile(all_coverage$specificity, 0.95, na.rm = TRUE))

p4 <- ggplot(
    all_coverage,
    aes(
        x = functional_group,
        y = TF,
        size = scales::oob_squish(specificity, c(0, size_cap))^1.5,
        color = cliff_clean
    )
) +
    geom_point(alpha = 0.85) +
    scale_size_continuous(
        range = size_range,
        limits = c(0.5, size_cap^1.5),
        breaks = size_breaks^1.5,
        labels = size_breaks,
        name = "Specificity"
    ) +
    cliff_scale +
    labs(
        title = paste0(unit_label, " Functional Specificity Across Disease Programs"),
        x = "Functional Group",
        y = unit_label
    ) +
    theme_bw(base_size = 11) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        ),
        panel.grid.minor = element_blank()
    )

#output one PDF: bubble plot (p4) first on a tall portrait page, then p1/p2/p3
#as 2x2 grids on landscape pages
grid_file <- tempfile(fileext = ".pdf")
p4_file   <- tempfile(fileext = ".pdf")

pdf(grid_file, width = 11, height = 8.5, onefile = TRUE)
gridExtra::grid.arrange(grobs = p1_list, ncol = 2, nrow = 2)
gridExtra::grid.arrange(grobs = p2_list, ncol = 2, nrow = 2)
gridExtra::grid.arrange(grobs = p3_list, ncol = 2, nrow = 2)
dev.off()

pdf(p4_file, width = 8.5, height = 12, onefile = TRUE)
print(p4)
dev.off()

outfile <- paste0("coverage_comparasion_plots_", unit_label, ".pdf")
qpdf::pdf_combine(c(p4_file, grid_file), outfile)
cat("Wrote", outfile, "\n")
