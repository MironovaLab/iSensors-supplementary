#!/usr/bin/env Rscript

# Refit the rice ARF model with the expanded local
# OSA ConjugationDeconjugation panel.
# Run with working directory = iSensors-supplementary/

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
})

analysis_dir <- "05-rice-conjugation-model"
input_dir <- file.path(analysis_dir, "in")
output_dir <- file.path(analysis_dir, "out")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

object_file <- "00-iSensors-objects/data/Wang_rice_iSensors.rds"
mapping_file <- paste0(
  "00-iSensors-objects/data/",
  "OSA-aux-trans-ConjugationDeconjugation_all_high_confidence.tsv"
)

message("Loading Wang rice atlas with existing iSensor scores")
rice <- readRDS(object_file)

required_metadata <- "tissue_cluster_names"
if (!required_metadata %in% colnames(rice[[]])) {
  stop("Missing metadata column: ", required_metadata)
}
required_sensors <- c(
  "OSA-aux-trans-ARF",
  "OSA-aux-trans-synthesis",
  "OSA-aux-trans-PAT",
  "OSA-aux-trans-IAA"
)
if (!all(required_sensors %in% rownames(rice[["iSensors_mean"]]))) {
  stop(
    "Missing existing iSensors: ",
    paste(setdiff(required_sensors, rownames(rice[["iSensors_mean"]])), collapse = ", ")
  )
}

mapping <- read.delim(mapping_file, stringsAsFactors = FALSE, check.names = FALSE)
required_mapping_columns <- c(
  "arabidopsis_gene", "rice_gene", "orthology_type", "orthology_confidence"
)
if (!all(required_mapping_columns %in% colnames(mapping))) {
  stop(
    "Mapping lacks columns: ",
    paste(setdiff(required_mapping_columns, colnames(mapping)), collapse = ", ")
  )
}

mapping <- mapping %>%
  filter(
    orthology_confidence == 1,
    orthology_type %in% c(
      "ortholog_one2one", "ortholog_one2many", "ortholog_many2many"
    ),
    nzchar(rice_gene)
  ) %>%
  distinct(arabidopsis_gene, rice_gene, .keep_all = TRUE) %>%
  mutate(present_in_Wang_RNA = rice_gene %in% rownames(rice[["RNA"]]))

conjugation_genes <- unique(mapping$rice_gene[mapping$present_in_Wang_RNA])
if (!length(conjugation_genes)) {
  stop("No expanded conjugation genes are present in Wang RNA")
}

message(
  "Expanded conjugation panel: ", length(conjugation_genes),
  " rice genes representing ",
  n_distinct(mapping$arabidopsis_gene[mapping$present_in_Wang_RNA]),
  " Arabidopsis genes"
)

write.table(
  mapping,
  file.path(input_dir, "expanded_conjugation_mapping_used.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# CalcSensors(mean) is exactly the column mean of the normalized RNA data for
# panel genes. Calculate only the replacement panel to avoid rescoring the
# entire 115k-cell atlas.
conjugation_matrix <- LayerData(
  rice[["RNA"]],
  layer = "data",
  features = conjugation_genes
)
new_conjugation <- Matrix::colMeans(conjugation_matrix)

sensor_matrix <- LayerData(
  rice[["iSensors_mean"]],
  layer = "data",
  features = required_sensors
)

cell_scores <- data.frame(
  group = as.character(rice[[required_metadata]][, 1]),
  ARF = as.numeric(sensor_matrix["OSA-aux-trans-ARF", ]),
  Synthesis = as.numeric(sensor_matrix["OSA-aux-trans-synthesis", ]),
  PAT = as.numeric(sensor_matrix["OSA-aux-trans-PAT", ]),
  Old_IAA = as.numeric(sensor_matrix["OSA-aux-trans-IAA", ]),
  Conjugation = as.numeric(new_conjugation),
  stringsAsFactors = FALSE
)

rm(rice, sensor_matrix, conjugation_matrix, new_conjugation)
gc()

# Average cell-level scores within the same 56 Organ--CellType groups used by
# the original Figure 7 analysis. Seurat::AverageExpression(layer = "data")
# exponentiates log-normalized values before averaging. Apply expm1 here to
# reproduce that convention exactly for both old and replacement panels.
df_wide <- cell_scores %>%
  group_by(group) %>%
  summarise(
    across(
      c(ARF, Synthesis, PAT, Old_IAA, Conjugation),
      \(x) mean(expm1(x), na.rm = TRUE)
    ),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  separate(
    group,
    into = c("Organ", "Tissue_Cleaned"),
    sep = "--",
    extra = "merge",
    fill = "right"
  ) %>%
  mutate(
    Organ = recode(
      Organ,
      Bud = "Bud",
      Flag = "Flag leaf",
      leaf = "Leaf",
      Root = "Root",
      SAM = "Shoot apical meristem",
      Seed = "Seed",
      SP = "Spikelet",
      ST = "Stem"
    )
  ) %>%
  drop_na(ARF, Synthesis, Conjugation, PAT)

m_new <- lm(ARF ~ 0 + Synthesis + Conjugation + PAT, data = df_wide)
m_old <- lm(ARF ~ 0 + Synthesis + Old_IAA + PAT, data = df_wide)

df_wide <- df_wide %>%
  mutate(
    ARF_hat = predict(m_new),
    residual = resid(m_new)
  )

write.csv(
  df_wide,
  file.path(input_dir, "rice_model_data_expanded_conjugation.csv"),
  row.names = FALSE
)

model_metrics <- function(model, model_name) {
  s <- summary(model)
  data.frame(
    model = model_name,
    n_groups = nobs(model),
    adjusted_R2 = s$adj.r.squared,
    RMSE = sqrt(mean(resid(model)^2)),
    residual_SE = s$sigma,
    AIC = AIC(model),
    stringsAsFactors = FALSE
  )
}

model_comparison <- bind_rows(
  model_metrics(m_old, "Old OSA IAA panel"),
  model_metrics(m_new, "Expanded ConjugationDeconjugation panel")
)
write.csv(
  model_comparison,
  file.path(input_dir, "old_vs_expanded_conjugation_model_metrics.csv"),
  row.names = FALSE
)

coefficient_table <- function(model, model_name) {
  tab <- summary(model)$coefficients
  data.frame(
    model = model_name,
    term = rownames(tab),
    estimate = tab[, "Estimate"],
    SE = tab[, "Std. Error"],
    t = tab[, "t value"],
    p_value = tab[, "Pr(>|t|)"],
    row.names = NULL
  )
}

coef_comparison <- bind_rows(
  coefficient_table(m_old, "Old OSA IAA panel"),
  coefficient_table(m_new, "Expanded ConjugationDeconjugation panel")
)
write.csv(
  coef_comparison,
  file.path(input_dir, "old_vs_expanded_conjugation_coefficients.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("EXPANDED CONJUGATION MODEL\n")
    print(summary(m_new))
    cat("\nOLD IAA MODEL RECOMPUTED ON IDENTICAL GROUPS\n")
    print(summary(m_old))
    cat("\nMODEL METRICS\n")
    print(model_comparison)
    cat("\nPANEL COVERAGE\n")
    cat("Rice genes used:", length(conjugation_genes), "\n")
    cat(
      "Arabidopsis source genes represented:",
      n_distinct(mapping$arabidopsis_gene[mapping$present_in_Wang_RNA]), "\n"
    )
  },
  file = file.path(input_dir, "rice_model_expanded_conjugation_summary.txt")
)

rmse <- sqrt(mean(df_wide$residual^2))
s_new <- summary(m_new)
fstat <- s_new$fstatistic
model_p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
model_label <- paste0(
  "adj. R² = ", format(round(s_new$adj.r.squared, 2), nsmall = 2),
  "\nRMSE = ", signif(rmse, 2),
  "\n", ifelse(model_p < 2.2e-16, "p < 2.2e-16", paste0("p = ", signif(model_p, 2)))
)

arf_range <- range(c(df_wide$ARF, df_wide$ARF_hat), na.rm = TRUE)
band <- tibble(x = seq(arf_range[1], arf_range[2], length.out = 200)) %>%
  mutate(lo = x - rmse, hi = x + rmse)

p_scatter <- ggplot(df_wide, aes(ARF, ARF_hat)) +
  geom_ribbon(
    data = band,
    aes(x = x, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.45
  ) +
  geom_point(size = 2, alpha = 0.8, colour = "#4d4d4d") +
  geom_abline(
    intercept = 0, slope = 1, linetype = "dashed",
    linewidth = 0.8
  ) +
  coord_equal(xlim = arf_range, ylim = arf_range) +
  annotate(
    "text", x = arf_range[1], y = arf_range[2],
    label = model_label, hjust = 0, vjust = 1,
    size = 3.2, colour = "grey30"
  ) +
  labs(
    x = "Observed ARF iSensor",
    y = "Predicted auxin level"
  ) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank())

coef_new <- coefficient_table(
  m_new, "Expanded ConjugationDeconjugation panel"
) %>%
  mutate(
    term = recode(term, Synthesis = "Synthesis", Conjugation = "Conjugation", PAT = "PAT"),
    term = factor(term, levels = rev(c("Synthesis", "Conjugation", "PAT"))),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    label_y = if_else(
      estimate >= 0,
      estimate + SE + 0.05 * max(abs(estimate + SE)),
      estimate - SE - 0.05 * max(abs(estimate - SE))
    )
  )

p_coefficients <- ggplot(coef_new, aes(term, estimate)) +
  geom_col(fill = "grey65", width = 0.6) +
  geom_errorbar(
    aes(ymin = estimate - SE, ymax = estimate + SE),
    width = 0.2, linewidth = 0.5
  ) +
  geom_text(aes(y = label_y, label = significance), size = 4) +
  labs(x = NULL, y = "Scaling coefficient") +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank())

df_wide <- df_wide %>%
  mutate(Tissue_Broad = case_when(
    str_detect(Tissue_Cleaned, regex("epidermis|guard.cell|atrichoblast|trichoblast", ignore_case = TRUE)) ~ "Epidermis",
    str_detect(Tissue_Cleaned, regex("vascular|procambium|xylem|phloem|stele", ignore_case = TRUE)) ~ "Vascular",
    str_detect(Tissue_Cleaned, regex("parenchyma|collenchymat", ignore_case = TRUE)) ~ "Parenchyma",
    str_detect(Tissue_Cleaned, regex("mesophyll", ignore_case = TRUE)) ~ "Mesophyll",
    str_detect(Tissue_Cleaned, regex("cortex", ignore_case = TRUE)) ~ "Cortex",
    str_detect(Tissue_Cleaned, regex("endodermis|exodermis", ignore_case = TRUE)) ~ "Endodermal",
    str_detect(Tissue_Cleaned, regex("columella|root.cap|lrc", ignore_case = TRUE)) ~ "Root cap / LRC",
    str_detect(Tissue_Cleaned, regex("meristem|initials|proliferating|primordium|apical|axis", ignore_case = TRUE)) ~ "Meristematic",
    str_detect(Tissue_Cleaned, regex("fiber|sclerenchyma", ignore_case = TRUE)) ~ "Fibers",
    str_detect(Tissue_Cleaned, regex("tapetum|pollen|ovule|ovary|stigmata|filament|lodicule|scutellum|plumule|mother.cell", ignore_case = TRUE)) ~ "Reproductive",
    str_detect(Tissue_Cleaned, regex("unknown", ignore_case = TRUE)) ~ "Unknown",
    TRUE ~ Tissue_Cleaned
  ))

rse <- s_new$sigma
make_forest <- function(summary_df, category, title) {
  plot_df <- summary_df %>%
    arrange(mean_residual) %>%
    mutate(category_plot = factor(.data[[category]], levels = .data[[category]]))

  ggplot(plot_df, aes(category_plot, mean_residual)) +
    geom_rect(
      ymin = -rse, ymax = rse, xmin = -Inf, xmax = Inf,
      fill = "grey90", alpha = 0.6, inherit.aes = FALSE
    ) +
    geom_rect(
      ymin = -0.3 * rse, ymax = 0.3 * rse, xmin = -Inf, xmax = Inf,
      fill = "grey75", alpha = 0.8, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_errorbar(
      aes(ymin = mean_residual - SE, ymax = mean_residual + SE),
      width = 0.25, linewidth = 0.5
    ) +
    geom_point(size = 2.5) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Mean residual (observed - predicted)",
      title = title
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      axis.text.y = element_text(size = 8),
      panel.grid = element_blank()
    )
}

organ_summary <- df_wide %>%
  group_by(Organ) %>%
  summarise(
    mean_residual = mean(residual),
    SE = sd(residual) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

tissue_summary <- df_wide %>%
  group_by(Tissue_Broad) %>%
  summarise(
    mean_residual = mean(residual),
    SE = sd(residual) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 2)

write.csv(
  organ_summary,
  file.path(input_dir, "rice_expanded_conjugation_organ_residuals.csv"),
  row.names = FALSE
)
write.csv(
  tissue_summary,
  file.path(input_dir, "rice_expanded_conjugation_tissue_residuals.csv"),
  row.names = FALSE
)

p_organ <- make_forest(organ_summary, "Organ", "Organ level")
p_tissue <- make_forest(tissue_summary, "Tissue_Broad", "Tissue level")
p_forest <- p_organ + p_tissue + plot_layout(widths = c(1, 1.25))

save_plot <- function(filename, plot, width, height) {
  ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    plot, width = width, height = height, dpi = 300, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    plot, width = width, height = height, dpi = 300, bg = "white"
  )
  ggsave(
    file.path(output_dir, paste0(filename, ".svg")),
    plot, width = width, height = height, dpi = 300, bg = "white"
  )
}

save_plot("Rice_ARF_scatter_expanded_conjugation", p_scatter, 4, 4)
save_plot("Rice_ARF_coefficients_expanded_conjugation", p_coefficients, 3.5, 3)
save_plot("Rice_ARF_residual_forests_expanded_conjugation", p_forest, 9, 5)

combined <- (p_scatter | p_coefficients) / p_forest +
  plot_layout(heights = c(1, 1.25)) +
  plot_annotation(tag_levels = "A")
save_plot("Rice_ARF_model_expanded_conjugation_combined", combined, 9, 9)

message("Completed expanded conjugation rice model.")
print(model_comparison)
