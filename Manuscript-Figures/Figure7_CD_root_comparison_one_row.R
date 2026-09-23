#!/usr/bin/env Rscript

# Regenerate the Figure 7 C/D root-comparison data with the original
# Figure7_root_species_tissue_average.R pipeline, then export the two panels
# as a standalone one-row figure.
#
# Run from iSensors-supplementary/:
# Rscript Manuscript-Figures/Figure7_CD_root_comparison_one_row.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

source_script <- "Manuscript-Figures/Figure7_root_species_tissue_average.R"
output_dir <- "Manuscript-Figures/out"
input_dir <- "Manuscript-Figures/in"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

message("Regenerating Figure 7 C/D data via ", source_script)
script <- readLines(source_script, warn = FALSE)
script <- gsub(
  "Figure7D_root_Arabidopsis_vs_rice_tissue_mean_iSensors.csv",
  "Figure7D_root_Arabidopsis_vs_rice_tissue_mean_iSensors_regenerated.csv",
  script,
  fixed = TRUE
)
script <- gsub(
  "Figure7D_root_Arabidopsis_vs_rice_rank_statistics.csv",
  "Figure7D_root_Arabidopsis_vs_rice_rank_statistics_regenerated.csv",
  script,
  fixed = TRUE
)
script <- gsub(
  "Figure7D_root_epidermis_subtype_mean_iSensors.csv",
  "Figure7D_root_epidermis_subtype_mean_iSensors_regenerated.csv",
  script,
  fixed = TRUE
)
script <- gsub(
  "Figure7D_root_Arabidopsis_vs_rice_tissue_mean_4sensors",
  "Figure7D_root_Arabidopsis_vs_rice_tissue_mean_4sensors_regenerated",
  script,
  fixed = TRUE
)
script <- gsub(
  "Figure7D_root_Arabidopsis_vs_rice_ARF_mirrored",
  "Figure7D_root_Arabidopsis_vs_rice_ARF_mirrored_regenerated",
  script,
  fixed = TRUE
)
script <- gsub(
  "Figure7D_root_species_rank_difference_heatmap",
  "Figure7D_root_species_rank_difference_heatmap_regenerated",
  script,
  fixed = TRUE
)
script <- gsub(
  "Figure7D_root_epidermis_subtype_means",
  "Figure7D_root_epidermis_subtype_means_regenerated",
  script,
  fixed = TRUE
)
eval(parse(text = script), envir = new.env(parent = globalenv()))

tissue_csv <- file.path(
  input_dir,
  "Figure7D_root_Arabidopsis_vs_rice_tissue_mean_iSensors_regenerated.csv"
)
rank_csv <- file.path(
  input_dir,
  "Figure7D_root_Arabidopsis_vs_rice_rank_statistics_regenerated.csv"
)

plot_data <- read.csv(tissue_csv, check.names = FALSE)
rank_stats <- read.csv(rank_csv, check.names = FALSE)

sensor_display <- c(
  ARF = "ARF",
  Synthesis = "Synthesis",
  ConjugationDeconjugation = "Conjugation",
  PAT = "PAT"
)
sensor_order <- names(sensor_display)
tissue_order <- c(
  "Epidermis",
  "Cortex",
  "Endodermis",
  "Pericycle (XPP)",
  "Xylem",
  "Root cap",
  "Stem cell niche"
)

# Harmonized with Figure8.R tissue palette where labels overlap, with closely
# related tones for Figure7-specific root tissues.
tissue_colors <- c(
  "Epidermis" = "#E67E22",
  "Cortex" = "#27AE60",
  "Endodermis" = "#2980B9",
  "Pericycle (XPP)" = "#8E44AD",
  "Xylem" = "#6A3D9A",
  "Root cap" = "#C0392B",
  "Stem cell niche" = "#16A085"
)

tissue_axis_labels <- c(
  "Epidermis" = "Epidermis",
  "Cortex" = "Cortex",
  "Endodermis" = "Endodermis",
  "Pericycle (XPP)" = "Pericycle",
  "Xylem" = "Xylem",
  "Root cap" = "Root cap",
  "Stem cell niche" = "Stem cell niche"
)

arf_mirror <- plot_data %>%
  filter(sensor == "ARF") %>%
  filter(tissue %in% tissue_order) %>%
  mutate(
    species = factor(
      species,
      levels = c("Arabidopsis\n(Shahan)", "Rice\n(Wang)")
    ),
    tissue = factor(tissue, levels = tissue_order),
    direction = if_else(species == "Arabidopsis\n(Shahan)", 1, -1),
    mirrored_mean = direction * relative_mean,
    mirrored_se_low = if_else(
      direction > 0,
      pmax(relative_mean - relative_se, 0),
      -(relative_mean + relative_se)
    ),
    mirrored_se_high = if_else(
      direction > 0,
      relative_mean + relative_se,
      -pmax(relative_mean - relative_se, 0)
    )
  )

p_arf_mirror <- ggplot(
  arf_mirror,
  aes(x = tissue, y = mirrored_mean, fill = tissue)
) +
  geom_hline(yintercept = 0, linewidth = 0.45, color = "grey20") +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.25) +
  geom_errorbar(
    aes(ymin = mirrored_se_low, ymax = mirrored_se_high),
    width = 0.17,
    linewidth = 0.35
  ) +
  annotate(
    "text",
    x = 0.55,
    y = 1.04,
    label = "Primary root (Arabidopsis)",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4.8
  ) +
  annotate(
    "text",
    x = 0.55,
    y = -1.04,
    label = "Crown root (rice)",
    hjust = 0,
    vjust = 0,
    fontface = "bold",
    size = 4.8
  ) +
  scale_fill_manual(values = tissue_colors, guide = "none", drop = FALSE) +
  scale_x_discrete(labels = tissue_axis_labels, drop = FALSE) +
  scale_y_continuous(
    limits = c(-1.08, 1.08),
    breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Relative mean ARF\niSensor signal",
    title = NULL,
    subtitle = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 13),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(18, 10, 35, 30)
  )

p_rank <- rank_stats %>%
  filter(tissue %in% tissue_order) %>%
  mutate(
    sensor = factor(sensor, levels = sensor_order, labels = sensor_display),
    tissue = factor(tissue, levels = tissue_order)
  ) %>%
  ggplot(aes(x = tissue, y = sensor, fill = rice_minus_arabidopsis)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = significance), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,
    name = "Rice - Arabidopsis\nmean percentile"
  ) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = tissue_axis_labels, drop = FALSE) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    plot.margin = margin(7, 8, 35, 5)
  )

combined <- plot_grid(
  p_arf_mirror,
  p_rank,
  nrow = 1,
  rel_widths = c(1.1, 1)
)

output_stem <- file.path(output_dir, "Figure7_CD_root_comparison_one_row")
ggsave(
  paste0(output_stem, ".pdf"),
  combined,
  width = 14.4,
  height = 4.1,
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  paste0(output_stem, ".png"),
  combined,
  width = 14.4,
  height = 4.1,
  dpi = 350,
  bg = "white"
)
ggsave(
  paste0(output_stem, ".svg"),
  combined,
  width = 14.4,
  height = 4.1,
  bg = "white"
)

message("Saved one-row Figure 7 C/D panel: ", output_stem, ".pdf/.png/.svg")
