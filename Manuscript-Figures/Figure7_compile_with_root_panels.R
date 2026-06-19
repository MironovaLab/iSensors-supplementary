#!/usr/bin/env Rscript

# Compile the existing Figure7_combined.pdf as the first row, with mirrored
# root ARF (C) and title-free rank-difference heatmap (D) below.
#
# Run from iSensors-supplementary/:
# Rscript Manuscript-Figures/Figure7_compile_with_root_panels.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(png)
})

output_dir <- "Manuscript-Figures/out"
input_dir <- "Manuscript-Figures/in"

top_pdf <- file.path(output_dir, "Figure7_combined.pdf")
arf_png <- file.path(
  output_dir,
  "Figure7D_root_Arabidopsis_vs_rice_ARF_mirrored.png"
)
stats_csv <- file.path(
  input_dir,
  "Figure7D_root_Arabidopsis_vs_rice_rank_statistics.csv"
)

output_stem <- file.path(output_dir, "Figure7_combined_with_root_comparison")
top_png_prefix <- file.path(tempdir(), "Figure7_top")
top_png <- paste0(top_png_prefix, ".png")

pdftoppm <- Sys.which("pdftoppm")
if (!nzchar(pdftoppm)) {
  stop("pdftoppm is required to rasterize Figure7_combined.pdf")
}

status <- system2(
  pdftoppm,
  args = c(
    "-png", "-r", "300", "-singlefile",
    shQuote(normalizePath(top_pdf, winslash = "/")),
    shQuote(top_png_prefix)
  )
)
if (status != 0 || !file.exists(top_png)) {
  stop("Failed to rasterize ", top_pdf)
}

top_raster <- readPNG(top_png)
arf_raster <- readPNG(arf_png)

rank_stats <- read.csv(stats_csv, check.names = FALSE)

sensor_display <- c(
  ARF = "ARF",
  Synthesis = "Synthesis",
  ConjugationDeconjugation = "Conjugation",
  PAT = "PAT"
)
sensor_order <- names(sensor_display)
tissue_order <- c(
  "Stem cell niche",
  "Epidermis",
  "Cortex",
  "Endodermis",
  "Pericycle (XPP)",
  "Xylem",
  "Root cap",
  "Lateral root primordium"
)

p_rank <- rank_stats |>
  transform(
    sensor = factor(sensor, levels = sensor_order, labels = sensor_display),
    tissue = factor(tissue, levels = tissue_order)
  ) |>
  ggplot(aes(x = tissue, y = sensor, fill = rice_minus_arabidopsis)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_text(aes(label = significance), size = 3.5, fontface = "bold") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    oob = scales::squish,
    name = "Rice \u2212 Arabidopsis\nmean percentile"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
    axis.text.y = element_text(face = "bold", size = 8.5),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7),
    plot.margin = margin(7, 8, 5, 5)
  )

raster_panel <- function(image) {
  ggplot() +
    annotation_custom(
      grid::rasterGrob(image, width = grid::unit(1, "npc"),
                       height = grid::unit(1, "npc"), interpolate = TRUE),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}

p_top <- raster_panel(top_raster)
p_arf <- raster_panel(arf_raster)

bottom <- plot_grid(
  p_arf,
  p_rank,
  nrow = 1,
  rel_widths = c(1.1, 1),
  labels = c("C", "D"),
  label_fontface = "bold",
  label_size = 15,
  label_x = c(0.01, 0.01),
  label_y = c(0.99, 0.99),
  hjust = 0,
  vjust = 1
)

combined <- plot_grid(
  p_top,
  bottom,
  ncol = 1,
  rel_heights = c(1, 0.9)
)

ggsave(
  paste0(output_stem, ".pdf"),
  combined,
  width = 12,
  height = 10.5,
  device = cairo_pdf,
  bg = "white"
)
ggsave(
  paste0(output_stem, ".png"),
  combined,
  width = 12,
  height = 10.5,
  dpi = 350,
  bg = "white"
)

message("Saved combined Figure 7: ", output_stem, ".pdf")
