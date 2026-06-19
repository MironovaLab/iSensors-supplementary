#!/usr/bin/env Rscript

# Figure 7: mean root iSensor signal per tissue in Shahan Arabidopsis and
# label-transferred Wang rice. PCC, PSE, and PPP cells are excluded.
#
# Run from iSensors-supplementary/:
# Rscript Manuscript-Figures/Figure7_root_species_tissue_average.R

suppressPackageStartupMessages({
  library(Seurat)
  library(iSensors)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

set.seed(100)

shahan_rds <- paste0(
  "00-iSensors-objects/data/",
  "shahan_balanced_full_root_full_genes.rds"
)
rice_rds <- paste0(
  "00-iSensors-objects/data/",
  "wang_rice_root_confident_Auxin_OSA_iSensors.rds"
)
output_dir <- "Manuscript-Figures/out"
input_dir <- "Manuscript-Figures/in"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

output_stem <- file.path(
  output_dir,
  "Figure7D_root_Arabidopsis_vs_rice_tissue_mean_4sensors"
)
output_arf_mirror_stem <- file.path(
  output_dir,
  "Figure7D_root_Arabidopsis_vs_rice_ARF_mirrored"
)
output_csv <- file.path(
  input_dir,
  "Figure7D_root_Arabidopsis_vs_rice_tissue_mean_iSensors.csv"
)
output_stats_csv <- file.path(
  input_dir,
  "Figure7D_root_Arabidopsis_vs_rice_rank_statistics.csv"
)
output_epidermis_csv <- file.path(
  input_dir,
  "Figure7D_root_epidermis_subtype_mean_iSensors.csv"
)

min_cells_per_tissue <- 10L

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
  "Lateral root primordium",
  "Stem cell niche"
)

tissue_colors <- c(
  "Epidermis" = "#009E73",
  "Cortex" = "#D55E00",
  "Endodermis" = "#E69F00",
  "Pericycle (XPP)" = "#56B4E9",
  "Xylem" = "#0072B2",
  "Procambium" = "#6A3D9A",
  "Root cap" = "#8C510A",
  "Lateral root primordium" = "#F46D43",
  "Stem cell niche" = "#D73027"
)

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

classify_tissue <- function(labels) {
  lower <- tolower(as.character(labels))

  excluded <- grepl("(^|_)pcc($|_)|(^|_)pse($|_)|(^|_)ppp($|_)", lower)
  tissue <- case_when(
    excluded ~ NA_character_,
    grepl("^(qc|csc|si|cei|lrcei)$|initial", lower) ~ "Stem cell niche",
    grepl("^lrp$", lower) ~ "Lateral root primordium",
    grepl("lrc|columella|^cscd$|^c1$|root.?cap|stressed rcc", lower) ~
      "Root cap",
    grepl("atrichoblast|trichoblast|epiderm", lower) ~ "Epidermis",
    grepl("cortex", lower) ~ "Cortex",
    grepl("endoderm", lower) ~ "Endodermis",
    grepl("xpp|pericycle", lower) ~ "Pericycle (XPP)",
    grepl("protoxylem|metaxylem|xylem", lower) ~ "Xylem",
    grepl("procamb", lower) ~ NA_character_,
    TRUE ~ NA_character_
  )
  tissue
}

polish_sensor <- function(x) {
  x <- sub("^ATH-aux-trans-", "", x)
  x <- sub("^OSA-aux-trans-", "", x)
  x <- sub("^OS-aux-trans-", "", x)
  recode(
    x,
    synthesis = "Synthesis",
    Synthesis = "Synthesis",
    PolarAuxinTransport = "PAT",
    PAT = "PAT",
    ARF = "ARF",
    ConjugationDeconjugation = "ConjugationDeconjugation"
  )
}

calculate_shahan_sensors <- function() {
  log_step("Loading balanced Shahan root object")
  shahan <- readRDS(shahan_rds)
  panels <- LoadSensors(
    setName = "Auxin",
    species = "ATH",
    hormone = "aux",
    customPanels = FALSE,
    random = FALSE
  )
  wanted <- c(
    "ATH-aux-trans-ARF",
    "ATH-aux-trans-Synthesis",
    "ATH-aux-trans-ConjugationDeconjugation",
    "ATH-aux-trans-PolarAuxinTransport"
  )
  panels[["panels"]] <- panels[["panels"]][
    intersect(wanted, names(panels[["panels"]]))
  ]

  genes <- unique(unlist(panels[["panels"]], use.names = FALSE))
  genes <- intersect(genes, rownames(shahan))
  counts <- LayerData(shahan[["RNA"]], layer = "counts")[genes, , drop = FALSE]

  small <- CreateSeuratObject(
    counts,
    assay = "RNA",
    meta.data = shahan[[]],
    project = "Shahan_Figure7D_tissue_mean"
  )
  small <- NormalizeData(small, verbose = FALSE)
  CalcSensors(
    small,
    seurLayer = "data",
    panelSet = panels,
    signals = "mean"
  )
}

extract_cell_values <- function(obj, labels, species) {
  sensor_mat <- as.matrix(
    LayerData(obj[["iSensors_mean"]], layer = "counts")
  )
  rownames(sensor_mat) <- polish_sensor(rownames(sensor_mat))
  sensor_mat <- sensor_mat[
    intersect(sensor_order, rownames(sensor_mat)), , drop = FALSE
  ]

  data <- as.data.frame(t(sensor_mat), check.names = FALSE) %>%
    tibble::rownames_to_column("cell") %>%
    mutate(
      label = as.character(labels),
      tissue = classify_tissue(labels)
    ) %>%
    pivot_longer(
      cols = all_of(rownames(sensor_mat)),
      names_to = "sensor",
      values_to = "value"
    ) %>%
    filter(!is.na(tissue)) %>%
    mutate(species = species)
}

summarize_tissues <- function(data) {
  data %>%
    group_by(tissue, sensor) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      sd_value = sd(value, na.rm = TRUE),
      cells = n(),
      se_value = sd_value / sqrt(cells),
      .groups = "drop"
    ) %>%
    filter(cells >= min_cells_per_tissue)
}

shahan <- calculate_shahan_sensors()
shahan_cells <- extract_cell_values(
  shahan,
  shahan$Atlas_new,
  "Arabidopsis\n(Shahan)"
)
shahan_summary <- summarize_tissues(shahan_cells) %>%
  mutate(species = "Arabidopsis\n(Shahan)")
rm(shahan)
invisible(gc())

log_step("Loading Wang rice root iSensors object")
rice <- readRDS(rice_rds)
rice_cells <- extract_cell_values(
  rice,
  rice$predicted.Atlas_new,
  "Rice\n(Wang)"
)
rice_summary <- summarize_tissues(rice_cells) %>%
  mutate(species = "Rice\n(Wang)")

all_cells <- bind_rows(shahan_cells, rice_cells)

plot_data <- bind_rows(shahan_summary, rice_summary) %>%
  group_by(species, sensor) %>%
  mutate(
    panel_max = max(mean_value, na.rm = TRUE),
    relative_mean = if_else(panel_max > 0, mean_value / panel_max, 0),
    relative_se = if_else(panel_max > 0, se_value / panel_max, 0)
  ) %>%
  ungroup() %>%
  mutate(
    species = factor(
      species,
      levels = c("Arabidopsis\n(Shahan)", "Rice\n(Wang)")
    ),
    sensor = factor(sensor, levels = sensor_order, labels = sensor_display),
    tissue = factor(tissue, levels = tissue_order)
  ) %>%
  arrange(species, sensor, tissue)

write.csv(plot_data, output_csv, row.names = FALSE)

# Exploratory cross-species comparison ---------------------------------------
# ATH and OSA panels have different gene composition and absolute ranges.
# Therefore each cell is converted to a percentile within its own
# species-sensor distribution before comparing tissue positions.
ranked_cells <- all_cells %>%
  group_by(species, sensor) %>%
  mutate(
    percentile = (rank(value, ties.method = "average") - 0.5) / n()
  ) %>%
  ungroup()

rank_stats <- ranked_cells %>%
  group_by(sensor, tissue) %>%
  group_modify(~ {
    ath <- .x$percentile[.x$species == "Arabidopsis\n(Shahan)"]
    rice_values <- .x$percentile[.x$species == "Rice\n(Wang)"]
    if (length(ath) < min_cells_per_tissue ||
        length(rice_values) < min_cells_per_tissue) {
      return(tibble())
    }
    test <- suppressWarnings(wilcox.test(
      rice_values,
      ath,
      exact = FALSE
    ))
    difference <- mean(rice_values) - mean(ath)
    difference_se <- sqrt(
      var(rice_values) / length(rice_values) +
        var(ath) / length(ath)
    )
    tibble(
      arabidopsis_cells = length(ath),
      rice_cells = length(rice_values),
      arabidopsis_mean_percentile = mean(ath),
      rice_mean_percentile = mean(rice_values),
      rice_minus_arabidopsis = difference,
      ci_low = difference - 1.96 * difference_se,
      ci_high = difference + 1.96 * difference_se,
      cliffs_delta = 2 * as.numeric(test$statistic) /
        (length(rice_values) * length(ath)) - 1,
      p_value = test$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adjusted < 0.001 ~ "***",
      p_adjusted < 0.01 ~ "**",
      p_adjusted < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

write.csv(rank_stats, output_stats_csv, row.names = FALSE)

p_difference <- rank_stats %>%
  mutate(
    sensor = factor(sensor, levels = sensor_order, labels = sensor_display),
    tissue = factor(tissue, levels = tissue_order)
  ) %>%
  ggplot(aes(x = tissue, y = sensor, fill = rice_minus_arabidopsis)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = significance), size = 4, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,
    name = "Rice \u2212 Arabidopsis\nmean percentile"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL,
    subtitle = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 50, hjust = 1, size = 8),
    axis.text.y = element_text(face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_text(size = 8),
    panel.grid = element_blank()
  )

ggsave(
  file.path(output_dir, "Figure7D_root_species_rank_difference_heatmap.pdf"),
  p_difference,
  width = 9,
  height = 4.5,
  device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Figure7D_root_species_rank_difference_heatmap.png"),
  p_difference,
  width = 9,
  height = 4.5,
  dpi = 350,
  bg = "white"
)

# Epidermis diagnostic -------------------------------------------------------
epidermis_data <- bind_rows(shahan_cells, rice_cells) %>%
  filter(tissue == "Epidermis") %>%
  mutate(
    epidermis_subtype = case_when(
      str_detect(label, regex("atrichoblast", ignore_case = TRUE)) ~
        "Atrichoblast",
      str_detect(label, regex("trichoblast", ignore_case = TRUE)) ~
        "Trichoblast",
      TRUE ~ "Other epidermis"
    )
  ) %>%
  group_by(species, sensor, epidermis_subtype) %>%
  summarise(
    mean_value = mean(value),
    cells = n(),
    se_value = sd(value) / sqrt(cells),
    .groups = "drop"
  ) %>%
  group_by(species, sensor) %>%
  mutate(relative_mean = mean_value / max(mean_value)) %>%
  ungroup()

write.csv(epidermis_data, output_epidermis_csv, row.names = FALSE)

p_epidermis <- epidermis_data %>%
  mutate(
    species = factor(
      species,
      levels = c("Arabidopsis\n(Shahan)", "Rice\n(Wang)")
    ),
    sensor = factor(sensor, levels = sensor_order, labels = sensor_display)
  ) %>%
  ggplot(aes(epidermis_subtype, relative_mean, fill = epidermis_subtype)) +
  geom_col(width = 0.72, color = "grey25", linewidth = 0.2) +
  facet_grid(sensor ~ species) +
  scale_fill_manual(
    values = c(
      Atrichoblast = "#009E73",
      Trichoblast = "#56B4E9",
      `Other epidermis` = "#999999"
    ),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    x = NULL,
    y = "Relative epidermal subtype mean",
    title = "Epidermal contribution to auxin-related iSensor signals"
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

ggsave(
  file.path(output_dir, "Figure7D_root_epidermis_subtype_means.pdf"),
  p_epidermis,
  width = 8,
  height = 9,
  device = cairo_pdf
)
ggsave(
  file.path(output_dir, "Figure7D_root_epidermis_subtype_means.png"),
  p_epidermis,
  width = 8,
  height = 9,
  dpi = 350,
  bg = "white"
)

p <- ggplot(
  plot_data,
  aes(x = tissue, y = relative_mean, fill = tissue)
) +
  geom_col(width = 0.76, color = "grey25", linewidth = 0.2) +
  geom_errorbar(
    aes(
      ymin = pmax(relative_mean - relative_se, 0),
      ymax = relative_mean + relative_se
    ),
    width = 0.18,
    linewidth = 0.35
  ) +
  facet_grid(
    rows = vars(sensor),
    cols = vars(species),
    scales = "fixed"
  ) +
  scale_fill_manual(values = tissue_colors, guide = "none", drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1.08),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  labs(
    x = NULL,
    y = "Relative mean iSensor signal",
    title = "Mean auxin-related iSensor activity across root tissues",
    subtitle = paste0(
      "PCC, PSE and PPP excluded; each species\u2013sensor panel scaled ",
      "to its highest tissue; mean \u00b1 SE"
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 7),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(0.7, "cm"),
    panel.spacing.y = grid::unit(0.35, "cm")
  )

ggsave(paste0(output_stem, ".pdf"), p, width = 12, height = 11, device = cairo_pdf)
ggsave(paste0(output_stem, ".svg"), p, width = 12, height = 11)
ggsave(
  paste0(output_stem, ".png"),
  p,
  width = 12,
  height = 11,
  dpi = 350,
  bg = "white"
)

# Mirrored ARF-only comparison -----------------------------------------------
arf_mirror <- plot_data %>%
  filter(sensor == "ARF") %>%
  mutate(
    direction = if_else(
      species == "Arabidopsis\n(Shahan)",
      1,
      -1
    ),
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
    label = "Arabidopsis (Shahan)",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4
  ) +
  annotate(
    "text",
    x = 0.55,
    y = -1.04,
    label = "Rice (Wang)",
    hjust = 0,
    vjust = 0,
    fontface = "bold",
    size = 4
  ) +
  scale_fill_manual(values = tissue_colors, guide = "none", drop = FALSE) +
  scale_y_continuous(
    limits = c(-1.08, 1.08),
    breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    x = NULL,
    y = "Relative mean ARF iSensor signal",
    title = "Mirrored ARF iSensor activity across root tissues",
    subtitle = "Arabidopsis shown above zero; rice shown below zero; mean \u00b1 SE"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 9.5, color = "grey30"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    axis.title.y = element_text(size = 10),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.3),
    panel.grid = element_blank(),
    plot.margin = margin(8, 10, 8, 10)
  )

ggsave(
  paste0(output_arf_mirror_stem, ".pdf"),
  p_arf_mirror,
  width = 9,
  height = 6.5,
  device = cairo_pdf
)
ggsave(
  paste0(output_arf_mirror_stem, ".svg"),
  p_arf_mirror,
  width = 9,
  height = 6.5
)
ggsave(
  paste0(output_arf_mirror_stem, ".png"),
  p_arf_mirror,
  width = 9,
  height = 6.5,
  dpi = 350,
  bg = "white"
)

log_step("Saved tissue-mean Figure 7 comparison: ", output_stem, ".pdf")
log_step("Saved mirrored ARF comparison: ", output_arf_mirror_stem, ".pdf")
log_step("Saved unscaled and relative means: ", output_csv)
log_step("Saved exploratory rank statistics: ", output_stats_csv)
