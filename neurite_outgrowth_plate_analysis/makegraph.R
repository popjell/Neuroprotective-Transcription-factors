library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ggsignif)

# ==========================================
# CONFIGURATION
# ==========================================
BASELINE_LABEL <- "NeuroBasal"  # Customize this label for baseline bars

# ==========================================
# 1. METADATA PARSING
# ==========================================

parse_metadata <- function(meta_path) {
  lines <- readLines(meta_path)
  lines <- trimws(lines)

  result <- list(
    baseline = character(0),
    drug_map  = list(),
    media_map = list()
  )

  section <- NULL

  for (line in lines) {
    if (line == "" || grepl("^#", line)) next

    if (grepl("^Baseline:", line)) {
      val <- trimws(sub("^Baseline:\\s*", "", line))
      if (val != "N/A" && val != "") {
        result$baseline <- trimws(unlist(strsplit(val, ",")))
      }
      next
    }

    if (grepl("^Drugged:", line)) { section <- "drugged"; next }
    if (grepl("^Media:", line))   { section <- "media";   next }

    if (section == "drugged" && grepl(":", line)) {
      parts <- strsplit(line, ":")[[1]]
      drug_name <- trimws(parts[1])
      cols <- as.numeric(trimws(unlist(strsplit(parts[2], ","))))
      for (col in cols) {
        result$drug_map[[as.character(col)]] <- drug_name
      }
    }

    if (section == "media" && grepl(":", line)) {
      parts <- strsplit(line, ":")[[1]]
      media_name <- trimws(parts[1])
      rows <- toupper(trimws(unlist(strsplit(parts[2], ","))))
      for (row in rows) {
        result$media_map[[tolower(row)]] <- media_name
      }
    }
  }

  return(result)
}

# ==========================================
# 2. PLATE DATA READING (summary xlsx)
# ==========================================
# Each "results" sheet has 3 measurement blocks stacked vertically:
#   Row i:   "Measurement" | "<name>"
#   Row i+1: header with plate column numbers (1-12)
#   Row i+2: row A (empty)
#   Rows i+3..i+8: rows B-G with values

read_plate_data <- function(file_path) {
  raw <- read_excel(file_path, sheet = "results", col_names = FALSE,
                    .name_repair = "minimal")
  mat <- as.matrix(raw)

  meas_rows <- which(mat[, 1] == "Measurement")
  blocks <- list()

  for (mr in meas_rows) {
    meas_name <- trimws(as.character(mat[mr, 2]))

    # Parse plate column numbers from header row (mr + 1)
    header <- mat[mr + 1, ]
    plate_cols <- rep(NA_integer_, ncol(mat))
    for (j in 2:ncol(mat)) {
      val <- trimws(as.character(header[j]))
      if (val != "" && !is.na(val)) {
        plate_cols[j] <- as.integer(round(as.numeric(val)))
      }
    }

    # Extract rows B-G (rows mr+3 through mr+8)
    row_letters <- c("b", "c", "d", "e", "f", "g")
    block_df <- expand.grid(Row = row_letters, Col = 1:12,
                            stringsAsFactors = FALSE)
    block_df$Value <- NA_real_

    for (i in 1:6) {
      data_row <- mat[mr + 2 + i, ]
      for (j in 2:ncol(mat)) {
        if (!is.na(plate_cols[j])) {
          pc <- plate_cols[j]
          val <- suppressWarnings(as.numeric(data_row[j]))
          idx <- which(block_df$Row == row_letters[i] & block_df$Col == pc)
          if (length(idx) > 0) block_df$Value[idx] <- val
        }
      }
    }

    block_df <- block_df[!is.na(block_df$Value), ]

    if (grepl("Total Cells", meas_name, ignore.case = TRUE)) {
      if (!"Total_Cells" %in% names(blocks)) blocks[["Total_Cells"]] <- block_df
    } else if (grepl("%.*Positive.*W2|Pct.*Positive", meas_name, ignore.case = TRUE)) {
      blocks[["Pct_Positive"]] <- block_df
    } else if (grepl("Positive W2", meas_name, ignore.case = TRUE)) {
      if (!"Positive_W2" %in% names(blocks)) blocks[["Positive_W2"]] <- block_df
    }
  }

  # If % Positive W2 is missing, calculate from counts
  if (!"Pct_Positive" %in% names(blocks) &&
      "Total_Cells" %in% names(blocks) &&
      "Positive_W2" %in% names(blocks)) {
    tc <- blocks[["Total_Cells"]] %>% rename(Total_Cells = Value)
    pw <- blocks[["Positive_W2"]] %>% rename(Positive_W2 = Value)
    calc <- left_join(tc, pw, by = c("Row", "Col")) %>%
      mutate(Value = (Positive_W2 / Total_Cells) * 100) %>%
      select(Row, Col, Value)
    blocks[["Pct_Positive"]] <- calc
  }

  return(blocks)
}

# ==========================================
# 3. BASELINE PARSING
# ==========================================

parse_baseline <- function(wells) {
  if (length(wells) == 0) {
    return(data.frame(Row = character(0), Col = integer(0),
                      stringsAsFactors = FALSE))
  }
  data.frame(
    Row = tolower(substr(wells, 1, 1)),
    Col = as.integer(substr(wells, 2, nchar(wells))),
    stringsAsFactors = FALSE
  )
}

# ==========================================
# 4. DRUG FACTOR ORDERING (DMSO first, then ascending conc.)
# ==========================================

order_drug_levels <- function(drug_names) {
  concs <- sapply(drug_names, function(d) {
    if (grepl("DMSO", d, ignore.case = TRUE)) return(0)
    m <- regmatches(d, regexpr("\\d+\\.?\\d*", d))
    if (length(m) > 0) as.numeric(m) else NA_real_
  })
  drug_names[order(concs)]
}

# ==========================================
# 5. PROCESS ALL PLATES
# ==========================================

meta_files <- list.files("Data", pattern = "metadata\\.txt$",
                         full.names = TRUE, recursive = TRUE)

if (length(meta_files) == 0) stop("No metadata files found in Data/")

all_well_data <- list()

for (mf in meta_files) {
  plate_id <- sub("_metadata\\.txt$", "", basename(mf))
  data_file <- file.path(dirname(mf), paste0(plate_id, "_data.xlsx"))

  if (!file.exists(data_file)) {
    warning("Data file not found: ", data_file, " -- skipping ", plate_id)
    next
  }

  cat("Processing:", plate_id, "\n")

  meta <- parse_metadata(mf)
  plate_blocks <- read_plate_data(data_file)

  if (!"Total_Cells" %in% names(plate_blocks) ||
      !"Pct_Positive" %in% names(plate_blocks)) {
    warning("Missing required measurement blocks in ", plate_id, " -- skipping")
    next
  }

  well_df <- plate_blocks[["Total_Cells"]] %>%
    rename(Total_Cells = Value) %>%
    left_join(plate_blocks[["Pct_Positive"]] %>% rename(Pct_Positive = Value),
              by = c("Row", "Col"))

  if ("Positive_W2" %in% names(plate_blocks)) {
    well_df <- well_df %>%
      left_join(plate_blocks[["Positive_W2"]] %>% rename(Positive_W2 = Value),
                by = c("Row", "Col"))
  }

  # Map conditions from metadata
  well_df <- well_df %>%
    mutate(
      Drug = sapply(as.character(Col), function(c) {
        if (c %in% names(meta$drug_map)) meta$drug_map[[c]] else NA_character_
      }),
      Media = sapply(Row, function(r) {
        if (r %in% names(meta$media_map)) meta$media_map[[r]] else NA_character_
      }),
      Plate = plate_id
    )

  # Mark baseline wells BEFORE filtering
  baseline_df <- parse_baseline(meta$baseline)
  well_df$Is_Baseline <- FALSE
  if (nrow(baseline_df) > 0) {
    well_df <- well_df %>%
      left_join(baseline_df %>% mutate(Is_Baseline = TRUE),
                by = c("Row", "Col"), suffix = c("", "_bl")) %>%
      mutate(Is_Baseline = coalesce(Is_Baseline_bl, Is_Baseline)) %>%
      select(-Is_Baseline_bl)
  }

  # Drop wells without a complete Drug + Media assignment (keep baseline even without Drug)
  well_df <- well_df %>% filter(!is.na(Media) & (!is.na(Drug) | Is_Baseline))
  # Assign "Baseline" to wells without a drug
  well_df <- well_df %>% mutate(Drug = ifelse(Is_Baseline & is.na(Drug), BASELINE_LABEL, Drug))

  all_well_data[[plate_id]] <- well_df
}

combined <- bind_rows(all_well_data)

# Drug levels: DMSO first, then ascending concentrations, Baseline last
drug_levels  <- order_drug_levels(unique(combined$Drug[combined$Drug != BASELINE_LABEL]))
drug_levels  <- c(drug_levels, BASELINE_LABEL)
media_levels <- sort(unique(combined$Media))

combined <- combined %>%
  mutate(
    Drug  = factor(Drug,  levels = drug_levels),
    Media = factor(Media, levels = media_levels)
  )

cat("\nCombined:", nrow(combined), "wells from", length(all_well_data), "plate(s)\n")
cat("Drugs:",  paste(levels(combined$Drug),  collapse = ", "), "\n")
cat("Media:",  paste(levels(combined$Media), collapse = ", "), "\n\n")

# ==========================================
# 6. NORMALIZATION (per-plate DMSO Control)
# ==========================================

dmso_means <- combined %>%
  filter(Drug == "DMSO", Media == "Control", !Is_Baseline) %>%
  group_by(Plate) %>%
  summarise(
    dmso_ctrl_mean_pw2   = mean(Positive_W2, na.rm = TRUE),
    dmso_ctrl_mean_pct   = mean(Pct_Positive, na.rm = TRUE),
    .groups = "drop"
  )

combined <- combined %>%
  left_join(dmso_means, by = "Plate") %>%
  mutate(
    PosW2_norm      = Positive_W2 / dmso_ctrl_mean_pw2,
    Pct_Positive_norm = Pct_Positive / dmso_ctrl_mean_pct
  ) %>%
  select(-dmso_ctrl_mean_pw2, -dmso_ctrl_mean_pct)

cat("Per-plate DMSO Control means:\n")
print(dmso_means)
cat("\n")

# ==========================================
# 7. METRICS DEFINITION
# ==========================================

metrics <- list(
  list(col = "Total_Cells",  label = "Total DAPI Cell Counts",
       ylab = "Mean Total Cells / Well (+ SD)",     fill = c("Control" = "#fb9a99", "PBMC" = "#b2df8a")),
  list(col = "Pct_Positive", label = "% Positive W2 (Beta-III Tubulin)",
       ylab = "Mean % Positive W2 (+ SD)",           fill = c("Control" = "#e31a1c", "PBMC" = "#33a02c")),
  list(col = "Positive_W2",  label = "Positive W2 Counts (Beta-III Tubulin)",
       ylab = "Mean Positive W2 Counts / Well (+ SD)",      fill = c("Control" = "#1f78b4", "PBMC" = "#6a3d9a")),
  list(col = "Pct_Positive_norm", label = "% Positive W2 Counts Normalized (vs Plate DMSO Ctrl)",
       ylab = "Fold Change vs DMSO Control (+ SD)",  fill = c("Control" = "#ff7f00", "PBMC" = "#cab2d6")),
  list(col = "PosW2_norm",   label = "Positive W2 Counts Normalized (vs Plate DMSO Ctrl)",
       ylab = "Fold Change vs DMSO Control (+ SD)",  fill = c("Control" = "#e5a43a", "PBMC" = "#984ea3"))
)

# ==========================================
# 8. MANUAL OUTLIER EXCLUSION
# ==========================================
# Add/remove rows here to flag wells as outliers.
# Reason is printed in the console log for traceability.

outliers <- tribble(
  ~Plate,     ~Row, ~Col, ~Reason,
  "plate_c2", "d",  5,    "Control 0.5uM - likely failed staining (21.4% vs 66% group mean, Grubbs G=2.51)",
  "plate_c1", "d",  6,    "Control 1.0uM - possible edge effect (38.5% vs 64.8% group mean, Grubbs G=2.39)",
  "plate_c2", "e",  4,    "PBMC 0.1uM - possible cell loss (11.4% vs 41% group mean, Grubbs G=2.25)",
  "plate_c2", "e",  5,    "PBMC 0.5uM - possible cell loss (20.7% vs 57.1% group mean, Grubbs G=2.36)"
)

if (nrow(outliers) > 0) {
  cat("Excluding", nrow(outliers), "outlier(s):\n")
  for (i in 1:nrow(outliers)) {
    cat(sprintf("  - %s row %s col %s: %s\n",
                outliers$Plate[i], outliers$Row[i], outliers$Col[i], outliers$Reason[i]))
  }
  cat("\n")

  combined <- combined %>%
    anti_join(outliers, by = c("Plate", "Row", "Col"))
}

# Separate baseline for plotting
baseline_data <- combined %>% filter(Is_Baseline)
combined_no_bl <- combined %>% filter(!Is_Baseline)

cat("Baseline wells:", nrow(baseline_data), "\n\n")

# ==========================================
# 9. TWO-WAY ANOVA FOR ALL METRICS
# ==========================================

run_anova <- function(formula, data) {
  fit <- tryCatch(aov(formula, data = data), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  tab <- summary(fit)[[1]]
  rn <- trimws(rownames(tab))
  if ("Media:Drug" %in% rn) {
    idx <- which(rn == "Media:Drug")
    round(tab[idx, "Pr(>F)"], 4)
  } else {
    NA_real_
  }
}

cat("Two-way ANOVA (Media x Drug) interaction p-values:\n")
for (m in metrics) {
  f <- as.formula(paste(m$col, "~ Media * Drug"))
  p <- run_anova(f, combined_no_bl)
  m$p_val <- p
  cat(sprintf("  %-35s p = %s\n", m$label, ifelse(is.na(p), "N/S", p)))
}
cat("\n")

# ==========================================
# 10. BAR PLOT GENERATOR
# ==========================================

make_bar_plot <- function(metric, combined, baseline_data, plot_data, anova_label) {
  col     <- metric$col
  label   <- metric$label
  ylab    <- metric$ylab
  fill    <- metric$fill

  p <- ggplot(plot_data, aes(x = Drug, y = Mean, fill = Media)) +
    geom_bar(stat = "identity", position = position_dodge(0.8),
             width = 0.7, color = "black") +
    geom_errorbar(aes(ymin = SD_Lower, ymax = SD_Upper),
                  position = position_dodge(0.8), width = 0.25) +
    geom_jitter(data = combined, aes(x = Drug, y = .data[[col]], group = Media),
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                size = 1.2, alpha = 0.6, shape = 21, color = "black", inherit.aes = FALSE)

  if (nrow(baseline_data) > 0) {
    bl_mean <- mean(baseline_data[[col]], na.rm = TRUE)
    bl_sd   <- sd(baseline_data[[col]], na.rm = TRUE)
    bl_df <- data.frame(Drug = BASELINE_LABEL, Mean = bl_mean,
                        SD_Lower = pmax(0, bl_mean - bl_sd),
                        SD_Upper = bl_mean + bl_sd)
    p <- p +
      geom_bar(data = bl_df, aes(x = Drug, y = Mean),
               stat = "identity", width = 0.5, fill = "gray60",
               color = "black", inherit.aes = FALSE) +
      geom_errorbar(data = bl_df, aes(x = Drug, ymin = SD_Lower, ymax = SD_Upper),
                    width = 0.2, inherit.aes = FALSE) +
      geom_jitter(data = baseline_data, aes(x = Drug, y = .data[[col]]),
                  position = position_jitter(width = 0.15), size = 1.2, alpha = 0.6, shape = 21,
                  color = "black", inherit.aes = FALSE)
  }

  p <- p +
    stat_compare_means(
      data = combined,
      aes(x = Drug, y = .data[[col]], group = Media),
      method = "t.test", label = "p.signif", vjust = -0.5, hide.ns = FALSE
    ) +
    labs(title = label, subtitle = anova_label,
         x = "Drug Treatment", y = ylab) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_fill_manual(values = fill) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major.x = element_blank(), legend.position = "top")
  p
}

# ==========================================
# 11. SURVIVAL % CALCULATOR
# ==========================================

calc_survival <- function(combined, col_name) {
  combined %>%
    group_by(Drug, Media) %>%
    summarise(Mean = mean(.data[[col_name]], na.rm = TRUE),
              SD   = sd(.data[[col_name]], na.rm = TRUE),
              N    = n(), .groups = "drop") %>%
    pivot_wider(names_from = Media, values_from = c(Mean, SD, N)) %>%
    mutate(
      Survival     = (Mean_PBMC / Mean_Control) * 100,
      SE_ratio     = Survival * sqrt(
        (SD_Control / Mean_Control)^2 / N_Control +
        (SD_PBMC    / Mean_PBMC)^2    / N_PBMC
      ),
      Survival_Lower = pmax(0, Survival - SE_ratio),
      Survival_Upper = Survival + SE_ratio
    )
}

calc_well_survival <- function(combined, col_name) {
  control_mean <- combined %>%
    filter(Media == "Control") %>%
    summarise(m = mean(.data[[col_name]], na.rm = TRUE)) %>%
    pull(m)
  combined %>%
    filter(Media == "PBMC") %>%
    mutate(Survival = (.data[[col_name]] / control_mean) * 100) %>%
    select(Drug, Survival)
}

calc_per_well_survival <- function(combined, col_name) {
  combined %>%
    filter(Media %in% c("Control", "PBMC")) %>%
    group_by(Drug) %>%
    mutate(ctrl_mean = mean(.data[[col_name]][Media == "Control"], na.rm = TRUE)) %>%
    ungroup() %>%
    filter(Media == "PBMC") %>%
    mutate(Survival = (.data[[col_name]] / ctrl_mean) * 100) %>%
    select(Drug, Survival)
}

get_survival_pairwise <- function(pw_surv) {
  fit <- aov(Survival ~ Drug, data = pw_surv)
  thsd <- TukeyHSD(fit, "Drug")$Drug
  comparisons <- list()
  labels <- c()
  for (cname in rownames(thsd)) {
    p <- thsd[cname, "p adj"]
    sig <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "n.s."
    if (sig != "n.s.") {
      parts <- strsplit(cname, "-")[[1]]
      comparisons <- c(comparisons, list(parts))
      labels <- c(labels, sig)
    }
  }
  list(comparisons = comparisons, labels = labels, anova_p = summary(fit)[[1]]["Drug", "Pr(>F)"])
}

make_survival_plot <- function(surv_data, metric, pairwise = NULL) {
  p <- ggplot(surv_data, aes(x = Drug, y = Survival)) +
    geom_bar(stat = "identity", width = 0.6, fill = "#4a90d9", color = "black") +
    geom_errorbar(aes(ymin = Survival_Lower, ymax = Survival_Upper),
                  width = 0.25) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "red", linewidth = 0.8)

  if (!is.null(pairwise) && !is.null(pairwise$anova_p)) {
    anova_lbl <- paste0("One-way ANOVA: p = ", round(pairwise$anova_p, 4))
    p <- p + labs(title = paste0(metric$label, " - Survival"),
                  subtitle = paste0("PBMC / Control (%)  |  ", anova_lbl))
  } else {
    p <- p + labs(title = paste0(metric$label, " - Survival"),
                  subtitle = "PBMC / Control (%)")
  }

  p <- p +
    xlab("Drug Treatment") + ylab("% Survival (PBMC / Control)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    theme_minimal(base_size = 14) +
    theme(panel.grid.major.x = element_blank())

  if (!is.null(pairwise) && length(pairwise$comparisons) > 0) {
    p <- p + geom_signif(
      comparisons = pairwise$comparisons,
      annotations = pairwise$labels,
      map_signif_level = FALSE,
      tip_length = 0.01,
      step_increase = 0.08
    )
  }
  p
}

# ==========================================
# 12. BUILD ALL PLOTS
# ==========================================

bar_plots  <- list()
surv_plots <- list()

for (m in metrics) {
  col <- m$col

  # Summary stats for bar plot (excluding baseline)
  summ <- combined_no_bl %>%
    group_by(Media, Drug) %>%
    summarise(
      Mean     = mean(.data[[col]], na.rm = TRUE),
      SD       = sd(.data[[col]], na.rm = TRUE),
      .groups  = "drop"
    ) %>%
    mutate(
      SD_Lower = pmax(0, Mean - SD),
      SD_Upper = Mean + SD
    )

  anova_label <- ifelse(is.na(m$p_val),
                        "Two-way ANOVA interaction: N/S",
                        paste("Two-way ANOVA interaction p =", m$p_val))

  bar_plots[[col]]  <- make_bar_plot(m, combined_no_bl, baseline_data, summ, anova_label)
  pw_surv <- calc_per_well_survival(combined_no_bl, col)
  pw_tests <- get_survival_pairwise(pw_surv)
  surv_plots[[col]] <- make_survival_plot(calc_survival(combined_no_bl, col), m, pairwise = pw_tests)
}

# ==========================================
# 13. OUTPUT
# ==========================================

# Print to plot viewer FIRST (before PDF device steals output)
for (m in metrics) {
  print(bar_plots[[m$col]])
  if (m$col != "Total_Cells") print(surv_plots[[m$col]])
}

# Combined PNG
combined_fig <- ggarrange(
  bar_plots[["Total_Cells"]],
  ggarrange(
    bar_plots[["Pct_Positive"]],       surv_plots[["Pct_Positive"]],
    bar_plots[["Positive_W2"]],        surv_plots[["Positive_W2"]],
    bar_plots[["Pct_Positive_norm"]],  surv_plots[["Pct_Positive_norm"]],
    bar_plots[["PosW2_norm"]],         surv_plots[["PosW2_norm"]],
    nrow = 4, ncol = 2
  ),
  nrow = 2, heights = c(1, 4)
)

ggsave("combined_anova_plot.png", plot = combined_fig,
       width = 16, height = 20, dpi = 300)
cat("Saved: combined_anova_plot.png\n")

# Individual plots PDF
pdf("all_plots.pdf", width = 10, height = 7)
for (m in metrics) {
  print(bar_plots[[m$col]])
  if (m$col != "Total_Cells") print(surv_plots[[m$col]])
}
dev.off()
cat("Saved: all_plots.pdf\n")