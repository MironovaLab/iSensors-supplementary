# FigureS6_Guo_atlas_heatmap.R
# Supplementary Figure 6: iSensor expression across the Guo 2025 multi-organ atlas
#
# Panel A: heatmap  — organs (y) × cell types (x), independent colour scale per sensor
# Panel B: stacked barplot — cell types (x), organs stacked (fill), no-root, 4 sensors
#
# Input:  00-iSensors-objects/data/Guo/guo_avg_exp_iSensors_mean.rds
# Output: Manuscript-Figures/out/FigureS6_Guo_atlas_heatmap.pdf/.svg
#
# Run from: iSensors-supplementary/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(scales)
})

# --------------------------------------------------------------------------- #
#  Shared helpers                                                              #
# --------------------------------------------------------------------------- #
polish_name <- function(x) {
  x <- gsub("^ATH-aux-cistrans-|^ATH-aux-cis-|^ATH-aux-trans-|^ATH-aux-reg-", "", x)
  x <- gsub("^AT-aux-cistrans-|^AT-aux-cis-|^AT-aux-trans-|^AT-aux-reg-", "", x)
  gsub("PolarAuxinTransport", "PAT", x)
}

organ_rename <- c(
  cauline        = "Cauline 42 DAG",
  D0_silique     = "Silique 0 DPA",
  D2_silique     = "Silique 2 DPA",
  D3_silique     = "Silique 3 DPA",
  D4_silique     = "Silique 4 DPA",
  D5_silique     = "Silique 5 DPA",
  D6_root        = "Root 6 DAG",
  D11_root       = "Root 11 DAG",
  Early_flower   = "Early flower",
  Middle_flower  = "Middle flower",
  Late_flower    = "Late flower",
  S1_leaf        = "Rosette 14 DAG",
  S2_leaf        = "Rosette 21 DAG",
  S3_leaf        = "Rosette 28 DAG",
  S4_leaf        = "Rosette 35 DAG",
  S5_leaf        = "Rosette 42 DAG",
  S6_leaf        = "Rosette 49 DAG",
  stem           = "Stem 42 DAG"
)

organ_order_all <- c(
  "Root 6 DAG", "Root 11 DAG",
  "Rosette 14 DAG", "Rosette 21 DAG", "Rosette 28 DAG",
  "Rosette 35 DAG", "Rosette 42 DAG", "Rosette 49 DAG",
  "Early flower", "Middle flower", "Late flower",
  "Silique 0 DPA", "Silique 2 DPA", "Silique 3 DPA",
  "Silique 4 DPA", "Silique 5 DPA",
  "Cauline 42 DAG", "Stem 42 DAG"
)

organ_order_noroot <- c(
  "Rosette 14 DAG", "Rosette 21 DAG", "Rosette 28 DAG",
  "Rosette 35 DAG", "Rosette 42 DAG", "Rosette 49 DAG",
  "Early flower", "Middle flower", "Late flower",
  "Silique 0 DPA", "Silique 2 DPA", "Silique 3 DPA",
  "Silique 4 DPA", "Silique 5 DPA",
  "Cauline 42 DAG", "Stem 42 DAG"
)

standardize_tissue <- function(x) {
  x <- str_trim(str_to_lower(x))
  x <- gsub("−|‑|‐|­", "-", x)   # normalise unicode dashes
  x <- str_remove(x, "-\\d+$")
  # Collapse near-duplicates before group rules
  x <- if_else(str_detect(x, "^guard.?cells?$"),                                 "guard cell", x)
  x <- if_else(str_detect(x, "companion.?cell"),                                  "companion cell", x)
  x <- if_else(str_detect(x, "bundle.?sheath"),                                   "bundle sheath", x)
  x <- if_else(str_detect(x, "lateral.?root.*(initial|iniational)"),              "lateral root initial cell", x)
  # Biological groupings
  x <- if_else(str_detect(x, "^(phloem|phloem[- ]parenchyma.*|protophloem)$"),   "phloem", x)
  x <- if_else(str_detect(x, "^(xylem|protoxylem.*|xylem/procambium)$"),          "xylem", x)
  x <- if_else(str_detect(x, "atrichoblast|trichoblast"),                          "trichoblast/atrichoblast", x)
  x <- if_else(str_detect(x, "epidermis") & !str_detect(x, "atrichoblast|trichoblast"), "epidermis", x)
  x <- if_else(str_detect(x, "seed coat"),                                         "seed coat", x)
  x <- if_else(str_detect(x, "^(pollen|sperm|tapetum|microsporocyte)$"),           "male gametophyte", x)
  x <- if_else(str_detect(x, "embryo|endosperm"),                                  "endosperm/embryo", x)
  x <- if_else(x %in% c("vascular", "vascular cells", "interfascicular fibers"),   "vascular (unspecified)", x)
  x <- if_else(x %in% c("g2m stage", "s stage"),                                   "dividing cell", x)
  x
}

sensors_order <- c("ARF", "Synthesis", "PAT", "ConjugationDeconjugation")

sensor_label <- c(
  ARF                     = "ARF",
  Synthesis               = "Synthesis",
  PAT                     = "PAT",
  ConjugationDeconjugation = "Conjugation /\nDeconjugation"
)

# --------------------------------------------------------------------------- #
#  Load and assemble data                                                      #
# --------------------------------------------------------------------------- #
guo_avg_exp <- readRDS("00-iSensors-objects/data/Guo/guo_avg_exp_iSensors_mean.rds")

all_data <- lapply(names(guo_avg_exp), function(org_key) {
  mat     <- guo_avg_exp[[org_key]]
  org_lab <- organ_rename[org_key]
  if (is.na(org_lab)) org_lab <- org_key
  rownames(mat) <- polish_name(rownames(mat))
  as.data.frame(as.matrix(mat)) |>
    tibble::rownames_to_column("sensor") |>
    pivot_longer(-sensor, names_to = "tissue", values_to = "value") |>
    mutate(organ = org_lab, tissue_simple = standardize_tissue(tissue))
}) |>
  bind_rows()

# Aggregate: max value per organ × tissue × sensor
all_agg <- all_data |>
  group_by(organ, tissue_simple, sensor) |>
  summarise(value = max(value, na.rm = TRUE), .groups = "drop")

# --------------------------------------------------------------------------- #
#  Panel A: heatmap                                                            #
# --------------------------------------------------------------------------- #
all_agg_A <- all_agg |>
  mutate(organ = factor(organ, levels = organ_order_all))

# Tissue order by prevalence across ALL organs
tissue_order_A <- all_agg_A |>
  distinct(organ, tissue_simple) |>
  count(tissue_simple) |>
  arrange(desc(n)) |>
  pull(tissue_simple)

all_agg_A <- mutate(all_agg_A,
  tissue_simple = factor(tissue_simple, levels = tissue_order_A)
)

sensors_A <- intersect(sensors_order, unique(all_agg_A$sensor))
message("Panel A sensors: ", paste(sensors_A, collapse = ", "))

make_heatmap_panel <- function(df, sensor_name, show_x = FALSE) {
  lim <- quantile(df$value, 0.99, na.rm = TRUE)
  ggplot(df, aes(x = tissue_simple, y = organ, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.25) +
    scale_fill_distiller(
      palette   = "YlOrRd",
      direction = 1,
      limits    = c(0, lim),
      oob       = squish,
      na.value  = "grey90",
      name      = "Mean\nexpression"
    ) +
    labs(title = sensor_name, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid      = element_blank(),
      axis.text.x     = if (show_x)
                          element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)
                        else
                          element_blank(),
      axis.ticks.x    = if (show_x) element_line(linewidth = 0.3) else element_blank(),
      axis.text.y     = element_text(size = 8),
      plot.title      = element_text(face = "bold", size = 10),
      legend.position = "right",
      legend.key.height = unit(0.55, "cm"),
      legend.title    = element_text(size = 7),
      legend.text     = element_text(size = 7)
    )
}

plots_A <- lapply(seq_along(sensors_A), function(i) {
  s <- sensors_A[i]
  make_heatmap_panel(
    filter(all_agg_A, sensor == s),
    sensor_name = s,
    show_x      = (i == length(sensors_A))
  )
})

p_A <- wrap_plots(plots_A, ncol = 1)

# --------------------------------------------------------------------------- #
#  Panel B: stacked barplot, tissues on x, organs as fill (no-root)           #
# --------------------------------------------------------------------------- #
all_agg_B <- all_agg |>
  filter(!organ %in% c("Root 6 DAG", "Root 11 DAG")) |>
  mutate(organ = factor(organ, levels = organ_order_noroot))

# Tissue order by prevalence across no-root organs
tissue_order_B <- all_agg_B |>
  distinct(organ, tissue_simple) |>
  count(tissue_simple) |>
  arrange(desc(n)) |>
  pull(tissue_simple)

all_agg_B <- mutate(all_agg_B,
  tissue_simple = factor(tissue_simple, levels = tissue_order_B)
)

sensors_B <- intersect(sensors_order, unique(all_agg_B$sensor))
message("Panel B sensors: ", paste(sensors_B, collapse = ", "))

organ_colors <- c(
  "Rosette 14 DAG" = "#c7e9b4",
  "Rosette 21 DAG" = "#7fcdbb",
  "Rosette 28 DAG" = "#41b6c4",
  "Rosette 35 DAG" = "#1d91c0",
  "Rosette 42 DAG" = "#225ea8",
  "Rosette 49 DAG" = "#0c2c84",
  "Early flower"   = "#fcc5c0",
  "Middle flower"  = "#f768a1",
  "Late flower"    = "#ae017e",
  "Silique 0 DPA"  = "#ffeda0",
  "Silique 2 DPA"  = "#feb24c",
  "Silique 3 DPA"  = "#f03b20",
  "Silique 4 DPA"  = "#bd0026",
  "Silique 5 DPA"  = "#7f0027",
  "Cauline 42 DAG" = "#b8a4c9",
  "Stem 42 DAG"    = "#4d4d4d"
)

make_bar_panel <- function(sensor_name, show_x = FALSE) {
  df <- all_agg_B |>
    filter(sensor == sensor_name) |>
    mutate(tissue_simple = factor(tissue_simple, levels = tissue_order_B))

  ggplot(df, aes(x = tissue_simple, y = value, fill = organ)) +
    geom_col(position = "stack", width = 0.85) +
    scale_fill_manual(
      values = organ_colors,
      name   = "Organ",
      drop   = FALSE
    ) +
    labs(x = NULL, y = sensor_label[sensor_name]) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x  = if (show_x)
                       element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)
                     else
                       element_blank(),
      axis.ticks.x = if (show_x) element_line(linewidth = 0.3) else element_blank(),
      axis.line    = element_line(linewidth = 0.4),
      axis.ticks.y = element_line(linewidth = 0.3),
      axis.title.y = element_text(size = 9, face = "bold"),
      panel.grid   = element_blank()
    )
}

plots_B <- lapply(seq_along(sensors_B), function(i) {
  make_bar_panel(sensors_B[i], show_x = (i == length(sensors_B)))
})

p_B <- wrap_plots(plots_B, ncol = 1, guides = "collect") &
  theme(legend.position = "right")

# --------------------------------------------------------------------------- #
#  Assemble and save                                                           #
# --------------------------------------------------------------------------- #
w_A    <- 8    # inches — narrower heatmap
w_B    <- 9    # inches — barplot with organ legend
height <- 3.5 * length(sensors_A)

combined <- (p_A | p_B) +
  plot_layout(widths = c(w_A, w_B)) +
  plot_annotation(
    tag_levels = list(c("A", rep("", length(sensors_A) - 1),
                        "B", rep("", length(sensors_B) - 1))),
    theme = theme(plot.tag = element_text(face = "bold", size = 13))
  )

out_dir <- "Manuscript-Figures/out"
ggsave(file.path(out_dir, "FigureS6_Guo_atlas_heatmap.pdf"),
       plot = combined, width = w_A + w_B, height = height,
       dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(out_dir, "FigureS6_Guo_atlas_heatmap.svg"),
       plot = combined, width = w_A + w_B, height = height,
       dpi = 300, bg = "white", limitsize = FALSE)

message(sprintf("Saved FigureS6_Guo_atlas_heatmap.pdf / .svg  (%g x %g in)",
                w_A + w_B, height))
message("Done.")
