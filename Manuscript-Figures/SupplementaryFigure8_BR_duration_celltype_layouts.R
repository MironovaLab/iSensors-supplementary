# SupplementaryFigure8_BR_duration_celltype_layouts.R
# Supplementary Figure S8: all-11-BR-sensor validation detail behind Figure 8B,
# plus realistic root layouts for the main BR sensor set.
#
#   A  BR duration response (BL 0.5-8h time course, single 100nM dose),
#      per-sample z-scored heatmap, all 11 BR panels.
#      = bottom panel of BRiSensors/out/BR_concentration_duration_celltype_figure_allsensors.pdf
#      (script 19), with the redundant/outlier BL_2h replicate "sc_5 (2h)"
#      column removed (BL_2h is still represented by sc_46 and sc_49).
#   B  BR cell-type response (BL 8h vs Control), per cell_type, heatmap,
#      all 11 BR panels.
#      = top panel of the same source figure.
#   C  Realistic longitudinal root layouts for 6 BR sensors (5 trans + reg-4hr-up),
#      recolored to the same green palette used for BR in Figure 8
#      (white -> #C7E9C0 -> #74C476 -> #238B45 -> #00441B), independent scale
#      per sensor.
#      = BRiSensors/out/BR_all_sensors_combined_row.pdf (script 11), recolored.
#
# Run with working directory = iSensors-supplementary/
# Input:  ../BRiSensors/out/GSE212230_iSensors_obj.rds
#         00-iSensors-objects/data/shahan_balanced_full_root_full_genes.rds
# Output: Manuscript-Figures/out/SupplementaryFigure8_BR_duration_celltype_layouts.pdf / .png / .svg

suppressPackageStartupMessages({
  library(Seurat); library(iSensors)
  library(dplyr); library(tidyr); library(tibble); library(stringr); library(forcats)
  library(ggplot2); library(patchwork); library(broom); library(scales)
  library(ggRootCellAtlas); library(ggPlantmap)
})

output_dir <- "Manuscript-Figures/out"
dir.create(output_dir, showWarnings = FALSE)
BR_DIR <- "../BRiSensors"

BR_GREEN <- c("white", "#C7E9C0", "#74C476", "#238B45", "#00441B")

# Ordered best -> worst (mean |r| BL-vs-BRZ, GSE212230) - same as script 19
SENSORS <- c(
  "ATH-br-reg-4hr-down", "ATH-br-reg-2hr-down", "ATH-br-reg-4hr-up",
  "ATH-br-trans-NegativeSignaling", "ATH-br-reg-2hr-up", "ATH-br-trans-TF-induced",
  "ATH-br-trans-Homeostasis", "ATH-br-trans-TF-repressed", "ATH-br-trans-TF",
  "ATH-br-trans-PositiveSignaling", "ATH-br-trans-Biosynthesis"
)
SENSOR_LABELS <- c(
  "ATH-br-reg-4hr-down"             = "reg-late-down",
  "ATH-br-reg-2hr-down"             = "reg-early-down",
  "ATH-br-reg-4hr-up"               = "reg-late-up",
  "ATH-br-trans-NegativeSignaling"  = "trans-NegSignaling",
  "ATH-br-reg-2hr-up"               = "reg-early-up",
  "ATH-br-trans-TF-induced"         = "trans-TF-induced",
  "ATH-br-trans-Homeostasis"        = "trans-Homeostasis",
  "ATH-br-trans-TF-repressed"       = "trans-TF-repressed",
  "ATH-br-trans-TF"                 = "trans-TF",
  "ATH-br-trans-PositiveSignaling"  = "trans-PosSignaling",
  "ATH-br-trans-Biosynthesis"       = "trans-Biosynthesis"
)

cat("Loading GSE212230 iSensors object...\n")
obj <- readRDS(file.path(BR_DIR, "out/GSE212230_iSensors_obj.rds"))
DefaultAssay(obj) <- "iSensors_mean"

meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, sample, cell_type, condition)

sens <- FetchData(obj, vars = SENSORS) %>% tibble::rownames_to_column("cell")

df <- meta %>% inner_join(sens, by = "cell")
rm(obj); invisible(gc())

# ══════════════════════════════════════════════════════════════════════════
# Panel B - Cell-type response (BL 8h vs Control), all 11 sensors
# ══════════════════════════════════════════════════════════════════════════
ct_df <- df %>% filter(condition %in% c("Control", "BL_8h"))

ct_avg <- ct_df %>%
  group_by(cell_type, condition) %>%
  summarise(across(all_of(SENSORS), ~ mean(.x, na.rm = TRUE)), n = dplyr::n(), .groups = "drop") %>%
  filter(n >= 10)

ct_wide <- ct_avg %>%
  select(cell_type, condition, all_of(SENSORS)) %>%
  pivot_longer(all_of(SENSORS), names_to = "sensor", values_to = "value") %>%
  pivot_wider(names_from = condition, values_from = value) %>%
  filter(!is.na(Control), !is.na(BL_8h))

rng <- range(c(ct_wide$Control, ct_wide$BL_8h), na.rm = TRUE)
use_ratio <- rng[1] >= 0
pc <- 1e-6
ct_wide <- ct_wide %>%
  mutate(effect = if (use_ratio) log2((BL_8h + pc) / (Control + pc)) else (BL_8h - Control),
         sensor_label = factor(SENSOR_LABELS[sensor], levels = SENSOR_LABELS))

lims_ct <- quantile(ct_wide$effect[is.finite(ct_wide$effect)], probs = c(0.02, 0.98), na.rm = TRUE)

p_celltype <- ggplot(ct_wide, aes(x = cell_type, y = sensor_label, fill = effect)) +
  geom_tile() +
  scale_fill_gradient2(low = "#008837", mid = "white", high = "#7b3294", midpoint = 0,
                        limits = lims_ct, oob = scales::squish,
                        name = if (use_ratio) expression(log[2]*"(BL 8h / Control)") else "BL 8h - Control") +
  scale_y_discrete(limits = rev(levels(ct_wide$sensor_label))) +
  labs(x = NULL, y = NULL, title = "BR cell-type response (BL 8h vs Control)") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_blank(), axis.line = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right", legend.key.width = unit(0.35, "cm"), legend.key.height = unit(0.9, "cm"),
        plot.title = element_text(size = 10, face = "bold"))

# ══════════════════════════════════════════════════════════════════════════
# Panel A - BR duration (BL 0.5-8h time course), all 11 sensors,
# "sc_5 (2h)" replicate removed (BL_2h still covered by sc_46, sc_49)
# ══════════════════════════════════════════════════════════════════════════
dur_df <- df %>%
  filter(condition %in% c("BL_0.5h", "BL_1h", "BL_2h", "BL_4h", "BL_8h")) %>%
  filter(sample != "sc_5") %>%
  mutate(hours = case_when(
    condition == "BL_0.5h" ~ 0.5, condition == "BL_1h" ~ 1, condition == "BL_2h" ~ 2,
    condition == "BL_4h" ~ 4, condition == "BL_8h" ~ 8))

dur_pb <- dur_df %>%
  group_by(sample, condition, hours) %>%
  summarise(across(all_of(SENSORS), ~ mean(.x, na.rm = TRUE)), n_cells = dplyr::n(), .groups = "drop")

dur_hm <- dur_pb %>%
  pivot_longer(all_of(SENSORS), names_to = "sensor", values_to = "value") %>%
  group_by(sensor) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup() %>%
  mutate(sensor_label = factor(SENSOR_LABELS[sensor], levels = SENSOR_LABELS))

samp_order <- dur_pb %>% distinct(sample, hours) %>% arrange(hours, sample) %>% pull(sample)
dur_hm <- dur_hm %>% mutate(sample = factor(sample, levels = samp_order))

hour_lab <- dur_pb %>% distinct(sample, hours) %>% arrange(hours) %>% mutate(lab = paste0(sample, " (", hours, "h)"))

p_duration <- ggplot(dur_hm, aes(x = sample, y = sensor_label, fill = value_z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0, name = "z-score") +
  scale_x_discrete(labels = setNames(hour_lab$lab, hour_lab$sample)) +
  scale_y_discrete(limits = rev(levels(dur_hm$sensor_label))) +
  labs(x = NULL, y = NULL, title = "BR duration response (BL 0.5-8h time course)") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 11),
        axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
        legend.position = "right", legend.key.width = unit(0.35, "cm"), legend.key.height = unit(0.9, "cm"),
        plot.title = element_text(size = 10, face = "bold"))

rm(df); invisible(gc())

# ══════════════════════════════════════════════════════════════════════════
# Panel C - Realistic longitudinal root layouts, green palette (BR),
# 5 trans panels + reg-4hr-up (= BR_all_sensors_combined_row.pdf selection)
# ══════════════════════════════════════════════════════════════════════════
cat("Loading shahan_balanced_full_root_full_genes.rds...\n")
sc_obj <- readRDS("00-iSensors-objects/data/shahan_balanced_full_root_full_genes.rds")
DefaultAssay(sc_obj) <- "RNA"
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)

# LoadSensors(customPanels=TRUE) looks for an "iSensors/" folder in the
# current working directory - the BR custom panels live in BRiSensors/iSensors/,
# so switch there temporarily just for this call.
.orig_wd <- getwd()
setwd(BR_DIR)
BRpanels <- LoadSensors(setName = "BR", defaultPanels = FALSE, customPanels = TRUE)
setwd(.orig_wd)
stopifnot(length(BRpanels$panels) > 0)
cat("Loaded", length(BRpanels$panels), "BR custom panels\n")
cat("Running CalcSensors for all BR panels on", ncol(sc_obj), "cells...\n")
iobj <- CalcSensors(sc_obj, seurLayer = "data", panelSet = BRpanels, signals = "mean")
rm(sc_obj); invisible(gc())

DefaultAssay(iobj) <- "iSensors_mean"
sensor_data <- LayerData(iobj, assay = "iSensors_mean", layer = "data")
clusters <- unique(Idents(iobj))

avg_mat <- data.frame()
for (clust in clusters) {
  cells_in_cluster <- WhichCells(iobj, idents = clust)
  if (length(cells_in_cluster) > 0) {
    cluster_means <- rowMeans(sensor_data[, cells_in_cluster, drop = FALSE])
    avg_mat <- rbind(avg_mat, data.frame(t(cluster_means)))
  }
}
avg_mat <- as.matrix(avg_mat)
colnames(avg_mat) <- rownames(sensor_data)
rownames(avg_mat) <- clusters
avg_mat <- t(avg_mat)
rm(iobj); invisible(gc())

br_trans <- rownames(avg_mat)[grep("^ATH-br-trans-", rownames(avg_mat))]
br_trans <- br_trans[!grepl("random|majortrend|TF-induced|TF-repressed", br_trans)]
br_sensors <- br_trans
reg_late_up <- "ATH-br-reg-4hr-up"
cat("Panel C sensors:", paste(br_sensors, collapse = ", "), "\n")
cat("Panel D sensor:", reg_late_up, "\n")

data("ggPm.At.longroot.longitudinal", package = "ggRootCellAtlas", envir = .GlobalEnv)

make_row_panel <- function(sensor_name) {
  short_name <- gsub("ATH-br-trans-|ATH-br-reg-", "", sensor_name)
  short_name <- gsub("4hr", "late", short_name)
  short_name <- gsub("2hr", "early", short_name)
  score_df <- data.frame(Atlas = colnames(avg_mat), Gene = as.numeric(avg_mat[sensor_name, ]), stringsAsFactors = FALSE)
  merged <- ggPlantmap.merge(ggPm.At.longroot.longitudinal, score_df, "Atlas")
  lims <- quantile(score_df$Gene, c(0.02, 0.98), na.rm = TRUE)
  ggPlantmap.heatmap(merged, Gene, linewidth = 0.15) +
    scale_fill_gradientn(colours = BR_GREEN, limits = lims, oob = squish, na.value = "grey92", name = "BR\nactivity") +
    labs(title = short_name) +
    theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"),
          plot.margin = margin(t = 3, r = 3, b = 3, l = 3, unit = "pt"),
          legend.position = "right", legend.key.height = unit(0.7, "cm"),
          legend.text = element_text(size = 6), legend.title = element_text(size = 7))
}

sensor_plots <- lapply(br_sensors, make_row_panel)
p_layouts <- wrap_plots(sensor_plots, ncol = length(sensor_plots))

# ══════════════════════════════════════════════════════════════════════════
# Panel D - reg-late-up: longitudinal view + all 6 transverse (radial)
# cross-sections (m1, m2, t, e1, e2, d), one shared color scale (same max)
# across all 7 views, styled like BRiSensors/07-shahan-BR-root-layout.R
# ══════════════════════════════════════════════════════════════════════════
data("ggPm.At.root.crosssection.m1", package = "ggRootCellAtlas", envir = .GlobalEnv)
data("ggPm.At.root.crosssection.m2", package = "ggRootCellAtlas", envir = .GlobalEnv)
data("ggPm.At.root.crosssection.t",  package = "ggRootCellAtlas", envir = .GlobalEnv)
data("ggPm.At.root.crosssection.e1", package = "ggRootCellAtlas", envir = .GlobalEnv)
data("ggPm.At.root.crosssection.e2", package = "ggRootCellAtlas", envir = .GlobalEnv)
data("ggPm.At.root.crosssection.d",  package = "ggRootCellAtlas", envir = .GlobalEnv)

reg_late_up_row <- avg_mat[reg_late_up, ]
score_df_D <- data.frame(Atlas = names(reg_late_up_row), Gene = as.numeric(reg_late_up_row), stringsAsFactors = FALSE)

layouts_D <- list(
  longitudinal = ggPm.At.longroot.longitudinal,
  m1 = ggPm.At.root.crosssection.m1, m2 = ggPm.At.root.crosssection.m2,
  t  = ggPm.At.root.crosssection.t,  e1 = ggPm.At.root.crosssection.e1,
  e2 = ggPm.At.root.crosssection.e2, d  = ggPm.At.root.crosssection.d
)

# One shared scale (same max) across the longitudinal view and all 6 cross-sections
matched_vals_D <- unlist(lapply(layouts_D, function(layout) {
  ggPlantmap.merge(layout, score_df_D, "Atlas")$Gene
}))
lims_D <- quantile(matched_vals_D, c(0.02, 0.98), na.rm = TRUE)

make_D_panel <- function(layout_data, subtitle_label) {
  merged <- ggPlantmap.merge(layout_data, score_df_D, "Atlas")
  ggPlantmap.heatmap(merged, Gene, linewidth = 0.15) +
    labs(subtitle = subtitle_label) +
    theme(plot.subtitle = element_text(size = 7, hjust = 0.5))
}

p_long_D <- make_D_panel(layouts_D$longitudinal, "longitudinal")
p_m1_D <- make_D_panel(layouts_D$m1, "m1")
p_m2_D <- make_D_panel(layouts_D$m2, "m2")
p_t_D  <- make_D_panel(layouts_D$t,  "transition")
p_e1_D <- make_D_panel(layouts_D$e1, "elongation 1")
p_e2_D <- make_D_panel(layouts_D$e2, "elongation 2")
p_d_D  <- make_D_panel(layouts_D$d,  "differentiated")

design_D <- "
1#
17
17
16
16
15
15
14
14
13
13
12
12
1#
1#
"
p_panelD <- (p_long_D + p_m1_D + p_m2_D + p_t_D + p_e1_D + p_e2_D + p_d_D) +
  plot_layout(design = design_D, guides = "collect") +
  plot_annotation(title = "reg-late-up",
                   theme = theme(plot.title = element_text(size = 9, hjust = 0.5, face = "bold"))) &
  scale_fill_gradientn(colours = BR_GREEN, limits = lims_D, oob = squish, na.value = "grey92", name = "BR\nactivity") &
  theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════
# Combine: (A | B) / (C | D), 2 rows, C and D same size, and save
# ══════════════════════════════════════════════════════════════════════════
p_row1 <- p_duration + p_celltype
p_row2 <- (wrap_elements(full = p_layouts) + wrap_elements(full = p_panelD)) +
  plot_layout(widths = c(1, 1))
p_final <- (p_row1 / p_row2) +
  plot_layout(heights = c(1, 1.3)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(output_dir, "SupplementaryFigure8_BR_duration_celltype_layouts.pdf"), plot = p_final, width = 16, height = 18, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "SupplementaryFigure8_BR_duration_celltype_layouts.png"), plot = p_final, width = 16, height = 18, dpi = 200, bg = "white")
ggsave(file.path(output_dir, "SupplementaryFigure8_BR_duration_celltype_layouts.svg"), plot = p_final, width = 16, height = 18, dpi = 300, bg = "white")
cat("Saved SupplementaryFigure8_BR_duration_celltype_layouts.pdf/.png/.svg\n")
