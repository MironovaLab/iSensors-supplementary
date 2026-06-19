#!/usr/bin/env Rscript

# Figure 7 root comparison:
# Four stacked iSensor barplots comparing Shahan Arabidopsis root (left) with
# Wang rice root (right). Root tissues are on the x-axis; developmental zones
# and anatomical subtypes form the stacked segments.
#
# Run from iSensors-supplementary/:
# Rscript Manuscript-Figures/Figure7_root_species_stacked_barplot.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(iSensors)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
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

output_pdf <- file.path(
  output_dir,
  "Figure7D_root_Arabidopsis_vs_rice_stacked_barplot_4sensors.pdf"
)
output_svg <- file.path(
  output_dir,
  "Figure7D_root_Arabidopsis_vs_rice_stacked_barplot_4sensors.svg"
)
output_png <- file.path(
  output_dir,
  "Figure7D_root_Arabidopsis_vs_rice_stacked_barplot_4sensors.png"
)
output_csv <- file.path(
  input_dir,
  "Figure7D_root_Arabidopsis_vs_rice_iSensor_means.csv"
)

min_cells_per_group <- 10L

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

sensor_display <- c(
  ARF = "ARF",
  Synthesis = "Synthesis",
  ConjugationDeconjugation = "Conjugation",
  PAT = "PAT"
)
sensor_order <- names(sensor_display)

tissue_order <- c(
  "Atrichoblast",
  "Trichoblast",
  "Cortex",
  "Endodermis",
  "PPP",
  "XPP",
  "PSE",
  "PCC",
  "Protoxylem",
  "Metaxylem",
  "Procambium",
  "Root cap",
  "Lateral root primordium",
  "Stem cell niche"
)

section_order <- c(
  "m1", "m2", "t", "e1", "e2", "d",
  "LRC", "Young LRC", "Dying LRC", "Columella", "CSCD", "C1",
  "LRP", "Stem cell niche"
)

parse_root_annotation <- function(labels) {
  labels <- as.character(labels)
  lower <- tolower(labels)

  is_scn <- grepl(
    "^(qc|csc|si|cei|lrcei)$|initial",
    lower
  )
  is_lrp <- grepl("^lrp$", lower)
  is_root_cap <- grepl(
    "lrc|columella|^cscd$|^c1$|root.?cap|stressed rcc",
    lower
  )

  tissue <- case_when(
    is_scn ~ "Stem cell niche",
    is_lrp ~ "Lateral root primordium",
    is_root_cap ~ "Root cap",
    grepl("atrichoblast", lower) ~ "Atrichoblast",
    grepl("trichoblast", lower) ~ "Trichoblast",
    grepl("cortex", lower) ~ "Cortex",
    grepl("endoderm", lower) ~ "Endodermis",
    grepl("ppp", lower) ~ "PPP",
    grepl("xpp", lower) ~ "XPP",
    grepl("pse", lower) ~ "PSE",
    grepl("pcc", lower) ~ "PCC",
    grepl("protoxylem", lower) ~ "Protoxylem",
    grepl("metaxylem", lower) ~ "Metaxylem",
    grepl("procamb", lower) ~ "Procambium",
    TRUE ~ NA_character_
  )

  zone <- str_match(lower, "_(m1|m2|t|e1|e2|d)$")[, 2L]
  subtype <- case_when(
    is_scn ~ "Stem cell niche",
    is_lrp ~ "LRP",
    grepl("dying lrc", lower) ~ "Dying LRC",
    grepl("young lrc", lower) ~ "Young LRC",
    grepl("^lrc$", lower) ~ "LRC",
    grepl("columella", lower) ~ "Columella",
    grepl("^cscd$", lower) ~ "CSCD",
    grepl("^c1$", lower) ~ "C1",
    !is.na(zone) ~ zone,
    TRUE ~ labels
  )

  data.frame(
    label = labels,
    tissue = tissue,
    subtype = subtype,
    stringsAsFactors = FALSE
  )
}

polish_sensor <- function(x) {
  x <- sub("^ATH-aux-trans-", "", x)
  x <- sub("^OSA-aux-trans-", "", x)
  x <- sub("^OS-aux-trans-", "", x)
  x <- recode(
    x,
    synthesis = "Synthesis",
    Synthesis = "Synthesis",
    `PolarAuxinTransport` = "PAT",
    PAT = "PAT",
    ARF = "ARF",
    ConjugationDeconjugation = "ConjugationDeconjugation"
  )
  x
}

calculate_shahan_sensors <- function() {
  log_step("Loading balanced Shahan full-root object")
  shahan <- readRDS(shahan_rds)
  ath_panels <- LoadSensors(
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
  wanted <- intersect(wanted, names(ath_panels[["panels"]]))
  ath_panels[["panels"]] <- ath_panels[["panels"]][wanted]

  panel_genes <- unique(unlist(ath_panels[["panels"]], use.names = FALSE))
  genes_found <- intersect(panel_genes, rownames(shahan))
  log_step("Shahan panel genes found: ", length(genes_found), "/", length(panel_genes))

  counts <- LayerData(shahan[["RNA"]], layer = "counts")[
    genes_found, , drop = FALSE
  ]
  meta <- shahan[[]]
  shahan_small <- CreateSeuratObject(
    counts = counts,
    assay = "RNA",
    meta.data = meta,
    project = "Shahan_root_Figure7D"
  )
  shahan_small <- NormalizeData(shahan_small, verbose = FALSE)
  shahan_small <- CalcSensors(
    shahan_small,
    seurLayer = "data",
    panelSet = ath_panels,
    signals = "mean"
  )
  shahan_small
}

extract_sensor_data <- function(obj, labels, species) {
  parsed <- parse_root_annotation(labels)
  sensor_mat <- as.matrix(
    LayerData(obj[["iSensors_mean"]], layer = "counts")
  )
  rownames(sensor_mat) <- polish_sensor(rownames(sensor_mat))
  sensor_mat <- sensor_mat[
    intersect(sensor_order, rownames(sensor_mat)), , drop = FALSE
  ]
  sensor_features <- rownames(sensor_mat)

  values <- as.data.frame(t(sensor_mat), check.names = FALSE) %>%
    tibble::rownames_to_column("cell") %>%
    bind_cols(parsed[, c("label", "tissue", "subtype")]) %>%
    pivot_longer(
      cols = all_of(sensor_features),
      names_to = "sensor",
      values_to = "value"
    ) %>%
    filter(!is.na(tissue), !is.na(subtype))

  values %>%
    group_by(tissue, subtype, sensor) %>%
    summarise(
      value = mean(value, na.rm = TRUE),
      cells = n(),
      .groups = "drop"
    ) %>%
    filter(cells >= min_cells_per_group) %>%
    mutate(species = species)
}

shahan_isensors <- calculate_shahan_sensors()
shahan_summary <- extract_sensor_data(
  shahan_isensors,
  shahan_isensors$Atlas_new,
  "Arabidopsis\n(Shahan)"
)
rm(shahan_isensors)
invisible(gc())

log_step("Loading confidently transferred Wang rice iSensors object")
rice <- readRDS(rice_rds)
rice_summary <- extract_sensor_data(
  rice,
  rice$predicted.Atlas_new,
  "Rice\n(Wang)"
)

plot_data <- bind_rows(shahan_summary, rice_summary) %>%
  mutate(
    species = factor(
      species,
      levels = c("Arabidopsis\n(Shahan)", "Rice\n(Wang)")
    ),
    sensor = factor(sensor, levels = sensor_order, labels = sensor_display),
    tissue = factor(tissue, levels = tissue_order),
    subtype = factor(
      subtype,
      levels = c(
        section_order,
        setdiff(sort(unique(subtype)), section_order)
      )
    )
  ) %>%
  group_by(species, sensor) %>%
  mutate(
    tissue_total = ave(value, tissue, FUN = sum),
    panel_max = max(tissue_total, na.rm = TRUE),
    relative_value = if_else(panel_max > 0, value / panel_max, 0)
  ) %>%
  ungroup() %>%
  arrange(species, sensor, tissue, subtype)

write.csv(plot_data, output_csv, row.names = FALSE)

subtype_levels <- levels(droplevels(plot_data$subtype))
zone_base <- c(
  m1 = "#542788",
  m2 = "#8073AC",
  t = "#2C7BB6",
  e1 = "#80CDC1",
  e2 = "#ABDDA4",
  d = "#D9EF8B"
)
subtype_colors <- setNames(
  grDevices::hcl.colors(length(subtype_levels), palette = "Dynamic"),
  subtype_levels
)
subtype_colors[names(zone_base)] <- zone_base
subtype_colors["Stem cell niche"] <- "#D73027"
subtype_colors["LRP"] <- "#F46D43"
subtype_colors["LRC"] <- "#8C510A"
subtype_colors["Young LRC"] <- "#BF812D"
subtype_colors["Dying LRC"] <- "#DFC27D"
subtype_colors["Columella"] <- "#01665E"
subtype_colors["CSCD"] <- "#35978F"
subtype_colors["C1"] <- "#80CDC1"

p <- ggplot(
  plot_data,
  aes(x = tissue, y = relative_value, fill = subtype)
) +
  geom_col(position = "stack", width = 0.82) +
  facet_grid(
    rows = vars(sensor),
    cols = vars(species),
    scales = "fixed"
  ) +
  scale_fill_manual(
    values = subtype_colors,
    name = "Zone / subtype",
    drop = TRUE
  ) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = "Relative stacked iSensor signal",
    title = "Root auxin-related iSensors across tissues and developmental subtypes",
    subtitle = paste0(
      "Each species\u2013sensor panel scaled to its highest tissue; ",
      "groups with \u2265 ", min_cells_per_group, " cells"
    )
  ) +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey30"),
    strip.background = element_blank(),
    strip.text.x = element_text(face = "bold", size = 10),
    strip.text.y = element_text(face = "bold", size = 9),
    axis.text.x = element_text(
      angle = 55,
      hjust = 1,
      vjust = 1,
      size = 7
    ),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.3),
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(0.7, "cm"),
    panel.spacing.y = grid::unit(0.35, "cm"),
    legend.position = "bottom",
    legend.text = element_text(size = 6),
    legend.title = element_text(face = "bold", size = 8),
    legend.key.size = grid::unit(0.28, "cm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave(output_pdf, p, width = 14, height = 12, device = cairo_pdf)
ggsave(output_svg, p, width = 14, height = 12)
ggsave(output_png, p, width = 14, height = 12, dpi = 350, bg = "white")

log_step("Saved Figure 7 root comparison: ", output_pdf)
log_step("Saved plotted summary data: ", output_csv)
