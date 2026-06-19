#!/usr/bin/env Rscript

# Linear ARF model for confidently reannotated Wang root cells using transferred
# Atlas_new subtypes and the expanded OSA ConjugationDeconjugation panel.
# Run with working directory = iSensors-supplementary/

suppressPackageStartupMessages({
  library(SeuratObject)
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

object_file <- paste0(
  "00-iSensors-objects/data/",
  "wang_rice_root_confident_Auxin_OSA_iSensors.rds"
)
min_cells_per_subtype <- 10L

message("Loading confidently reannotated Wang root object")
rice <- readRDS(object_file)

label_column <- "predicted.Atlas_new"
required_sensors <- c(
  ARF = "OSA-aux-trans-ARF",
  Synthesis = "OSA-aux-trans-synthesis",
  Conjugation = "OSA-aux-trans-ConjugationDeconjugation",
  PAT = "OSA-aux-trans-PAT"
)

if (!label_column %in% colnames(rice[[]])) {
  stop("Missing transferred annotation: ", label_column)
}
if (!all(required_sensors %in% rownames(rice[["iSensors_mean"]]))) {
  stop(
    "Missing iSensors: ",
    paste(
      setdiff(required_sensors, rownames(rice[["iSensors_mean"]])),
      collapse = ", "
    )
  )
}

sensor_matrix <- LayerData(
  rice[["iSensors_mean"]],
  layer = "data",
  features = unname(required_sensors)
)
rownames(sensor_matrix) <- names(required_sensors)

cell_scores <- as.data.frame(t(as.matrix(sensor_matrix)), check.names = FALSE) %>%
  mutate(
    subtype = as.character(rice[[label_column]][, 1]),
    cell = colnames(rice)
  )

rm(rice, sensor_matrix)
gc()

subtype_counts <- cell_scores %>%
  count(subtype, name = "n_cells", sort = TRUE) %>%
  mutate(in_primary_model = n_cells >= min_cells_per_subtype)

write.csv(
  subtype_counts,
  file.path(input_dir, "reannotated_root_subtype_counts.csv"),
  row.names = FALSE
)

# Match the Seurat::AverageExpression(layer = "data") convention used in
# Figure 7: exponentiate log-normalized scores before averaging.
df_wide <- cell_scores %>%
  group_by(subtype) %>%
  summarise(
    across(
      all_of(names(required_sensors)),
      \(x) mean(expm1(x), na.rm = TRUE)
    ),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  filter(n_cells >= min_cells_per_subtype)

if (nrow(df_wide) <= length(required_sensors) + 1L) {
  stop("Too few robust subtypes for the requested model")
}

m_root <- lm(ARF ~ 0 + Synthesis + Conjugation + PAT, data = df_wide)

df_wide <- df_wide %>%
  mutate(
    ARF_hat = predict(m_root),
    residual = resid(m_root),
    tissue = case_when(
      str_detect(subtype, regex("atrichoblast|trichoblast", ignore_case = TRUE)) ~ "Epidermis",
      str_detect(subtype, regex("cortex", ignore_case = TRUE)) ~ "Cortex",
      str_detect(subtype, regex("endoderm", ignore_case = TRUE)) ~ "Endodermis",
      str_detect(subtype, regex("xpp", ignore_case = TRUE)) ~ "Pericycle (XPP)",
      str_detect(subtype, regex("protoxylem|metaxylem", ignore_case = TRUE)) ~ "Xylem",
      str_detect(subtype, regex("lrc", ignore_case = TRUE)) ~ "Root cap",
      subtype == "LRP" ~ "Lateral root primordium",
      subtype %in% c("SI", "CEI", "LRCEI") ~ "Stem cell niche",
      str_detect(subtype, regex("procamb", ignore_case = TRUE)) ~ "Procambium",
      str_detect(subtype, regex("ppp", ignore_case = TRUE)) ~ "PPP",
      str_detect(subtype, regex("pcc", ignore_case = TRUE)) ~ "PCC",
      str_detect(subtype, regex("pse", ignore_case = TRUE)) ~ "PSE",
      TRUE ~ "Other"
    )
  )

write.csv(
  df_wide,
  file.path(input_dir, "reannotated_root_model_data.csv"),
  row.names = FALSE
)

s <- summary(m_root)
rmse <- sqrt(mean(df_wide$residual^2))
fstat <- s$fstatistic
model_p <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)

coef_table <- data.frame(
  term = rownames(s$coefficients),
  estimate = s$coefficients[, "Estimate"],
  SE = s$coefficients[, "Std. Error"],
  t = s$coefficients[, "t value"],
  p_value = s$coefficients[, "Pr(>|t|)"],
  row.names = NULL
)
write.csv(
  coef_table,
  file.path(input_dir, "reannotated_root_model_coefficients.csv"),
  row.names = FALSE
)

tissue_summary <- df_wide %>%
  group_by(tissue) %>%
  summarise(
    mean_residual = mean(residual),
    SE = if_else(n() > 1, sd(residual) / sqrt(n()), NA_real_),
    n_subtypes = n(),
    cells = sum(n_cells),
    .groups = "drop"
  )
write.csv(
  tissue_summary,
  file.path(input_dir, "reannotated_root_tissue_residuals.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("REANNOTATED WANG ROOT MODEL\n")
    print(summary(m_root))
    cat("\nMODEL DATA\n")
    cat("Cells:", sum(df_wide$n_cells), "\n")
    cat("Robust transferred subtypes:", nrow(df_wide), "\n")
    cat("Minimum cells per subtype:", min_cells_per_subtype, "\n")
    cat("RMSE:", rmse, "\n")
    cat("\nSUBTYPE RESIDUALS\n")
    print(df_wide %>% arrange(residual))
    cat("\nTISSUE RESIDUALS\n")
    print(tissue_summary %>% arrange(mean_residual))
  },
  file = file.path(input_dir, "reannotated_root_model_summary.txt")
)

model_label <- paste0(
  "adj. R² = ", format(round(s$adj.r.squared, 2), nsmall = 2),
  "\nRMSE = ", signif(rmse, 2),
  "\n", ifelse(
    model_p < 2.2e-16,
    "p < 2.2e-16",
    paste0("p = ", signif(model_p, 2))
  )
)

arf_range <- range(c(df_wide$ARF, df_wide$ARF_hat), na.rm = TRUE)
band <- tibble(x = seq(arf_range[1], arf_range[2], length.out = 200)) %>%
  mutate(lo = x - rmse, hi = x + rmse)

p_scatter <- ggplot(df_wide, aes(ARF, ARF_hat)) +
  geom_ribbon(
    data = band,
    aes(x = x, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    fill = "grey80", alpha = 0.45
  ) +
  geom_point(aes(size = n_cells), alpha = 0.78, colour = "#4d4d4d") +
  geom_abline(
    intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8
  ) +
  coord_equal(xlim = arf_range, ylim = arf_range) +
  annotate(
    "text", x = arf_range[1], y = arf_range[2],
    label = model_label, hjust = 0, vjust = 1,
    size = 3.2, colour = "grey30"
  ) +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  labs(
    x = "Observed ARF iSensor",
    y = "Predicted auxin level"
  ) +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank())

coef_plot_data <- coef_table %>%
  mutate(
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

p_coefficients <- ggplot(coef_plot_data, aes(term, estimate)) +
  geom_col(fill = "grey65", width = 0.6) +
  geom_errorbar(
    aes(ymin = estimate - SE, ymax = estimate + SE),
    width = 0.2, linewidth = 0.5
  ) +
  geom_text(aes(y = label_y, label = significance), size = 4) +
  labs(x = NULL, y = "Scaling coefficient") +
  theme_classic(base_size = 11) +
  theme(panel.grid = element_blank())

rse <- s$sigma
forest_background <- list(
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = -rse, ymax = rse,
    fill = "grey90", alpha = 0.6
  ),
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = -0.3 * rse, ymax = 0.3 * rse,
    fill = "grey75", alpha = 0.8
  )
)

subtype_plot_data <- df_wide %>%
  arrange(residual) %>%
  mutate(subtype = factor(subtype, levels = subtype))

p_subtype <- ggplot(subtype_plot_data, aes(subtype, residual)) +
  forest_background +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(aes(size = n_cells), colour = "black") +
  scale_size_continuous(range = c(2, 5), guide = "none") +
  coord_flip() +
  labs(
    x = NULL,
    y = "Residual (observed - predicted)",
    title = "Transferred subtype level"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank()
  )

tissue_plot_data <- tissue_summary %>%
  filter(n_subtypes >= 2) %>%
  arrange(mean_residual) %>%
  mutate(tissue = factor(tissue, levels = tissue))

p_tissue <- ggplot(tissue_plot_data, aes(tissue, mean_residual)) +
  forest_background +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbar(
    aes(ymin = mean_residual - SE, ymax = mean_residual + SE),
    width = 0.25
  ) +
  geom_point(size = 2.5) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Mean residual (observed - predicted)",
    title = "Broad tissue level"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

p_forest <- p_subtype + p_tissue + plot_layout(widths = c(1.6, 1))

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

save_plot("Reannotated_root_ARF_scatter", p_scatter, 4, 4)
save_plot("Reannotated_root_ARF_coefficients", p_coefficients, 3.5, 3)
save_plot("Reannotated_root_ARF_residual_forests", p_forest, 9, 6)

combined <- (p_scatter | p_coefficients) / p_forest +
  plot_layout(heights = c(1, 1.45)) +
  plot_annotation(tag_levels = "A")
save_plot("Reannotated_root_ARF_model_combined", combined, 9, 10)

message("Completed reannotated root model.")
message(
  "Model uses ", nrow(df_wide), " subtypes and ",
  sum(df_wide$n_cells), " cells; adjusted R2 = ",
  round(s$adj.r.squared, 3), ", RMSE = ", round(rmse, 4)
)
