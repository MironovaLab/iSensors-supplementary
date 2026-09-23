# 04-GSE212230-sc-analysis.R
# Runs BR + Auxin iSensors on GSE212230 (Nolan lab, 79,982 cells)
# BL time course: Control | BRZ | BL 0.5h | 1h | 2h | 4h | 8h
#
# Dataset: Seurat v3.1.5 object — use SeuratObject-only loading to avoid
#          S4 dispatch crash ("no slot of name 'images'") when full Seurat
#          is loaded against an old object.
#
# Input:  iSensors-supplementary/00-iSensors-objects/data/GSE212230_BR_time_course_inner.rds.gz
# Output: out/GSE212230_line_plots.pdf         — BR + Aux sensors across time
#         out/GSE212230_zone_heatmap.pdf        — sensors × spatial zone × condition
#         out/GSE212230_celltype_heatmap.pdf    — sensors × cell type × condition
#         out/GSE212230_iSensors_obj.rds        — full scored Seurat v5 object

suppressPackageStartupMessages(library(SeuratObject))

setwd("C:/Users/Victoria Mironova/Code/DigitalSensor-Toolbox/BRiSensors")
dir.create("out", showWarnings = FALSE)

# ── 1. Extract from Seurat v3.1.5 ────────────────────────────────────────────
cat("Loading GSE212230 (SeuratObject only — avoids v3.1.5 S4 dispatch crash)...\n")
f_in <- "00-iSensors-objects/data/GSE212230_BR_time_course_inner.rds.gz"
old_obj    <- readRDS(f_in)
rna_counts <- old_obj@assays$RNA@counts
meta_data  <- old_obj@meta.data
cat("  Cells:", nrow(meta_data), " | Genes:", nrow(rna_counts), "\n")
cat("  Meta columns:", paste(colnames(meta_data), collapse = ", "), "\n")
cat("  Samples:\n"); print(table(meta_data$sample))
cat("  Treatment × time_trt:\n")
print(table(meta_data$treatment, meta_data$time_trt, useNA = "ifany"))
rm(old_obj); invisible(gc())
cat("  Old object removed. RAM freed.\n")

# ── 2. Build clean condition label from treatment + time_trt ─────────────────
# time_trt values: "" (Control), "BRZ", "0.5_hour_BL", "1_hour_BL",
#                  "2_hour_BL", "4_hour_BL", "8_hour_BL"
meta_data$condition <- dplyr::case_when(
  meta_data$treatment == "Control"                  ~ "Control",
  meta_data$treatment == "BRZ"                      ~ "BRZ",
  meta_data$time_trt  == "0.5_hour_BL"              ~ "BL_0.5h",
  meta_data$time_trt  == "1_hour_BL"                ~ "BL_1h",
  meta_data$time_trt  == "2_hour_BL"                ~ "BL_2h",
  meta_data$time_trt  == "4_hour_BL"                ~ "BL_4h",
  meta_data$time_trt  == "8_hour_BL"                ~ "BL_8h",
  TRUE                                               ~ "Other"
)

condition_order <- c("Control", "BRZ", "BL_0.5h", "BL_1h", "BL_2h", "BL_4h", "BL_8h")
condition_order <- intersect(condition_order, unique(meta_data$condition))

cat("\n  Condition counts:\n")
print(table(meta_data$condition, meta_data$sample, useNA = "ifany"))

# ── 3. Build Seurat v5 object ────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat); library(iSensors); library(ggplot2)
  library(dplyr); library(tidyr); library(tibble)
  library(forcats); library(patchwork); library(scales)
})

cat("\nCreating Seurat v5 object...\n")
sc_obj <- CreateSeuratObject(counts = rna_counts, meta.data = meta_data)
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
rm(rna_counts); invisible(gc())
cat("  Object built and normalized. RAM freed.\n")

# Ordered factors
zone_levels <- c(
  "Proliferation Domain", "Transition Domain",
  "Proximal Columella",   "Distal Columella",
  "Proximal Lateral Root Cap", "Distal Lateral Root Cap",
  "Elongation", "Maturation"
)
zone_labels <- c(
  "Proliferation Domain"      = "Prolif.",
  "Transition Domain"         = "Transition",
  "Proximal Columella"        = "Prox. Col.",
  "Distal Columella"          = "Dist. Col.",
  "Proximal Lateral Root Cap" = "Prox. LRC",
  "Distal Lateral Root Cap"   = "Dist. LRC",
  "Elongation"                = "Elongation",
  "Maturation"                = "Maturation"
)
celltype_levels <- c(
  "Quiescent Center", "Columella", "Lateral Root Cap",
  "Atrichoblast", "Trichoblast", "Cortex", "Endodermis",
  "Phloem", "Procambium", "Xylem", "Pericycle"
)

sc_obj$condition <- factor(sc_obj$condition, levels = condition_order)
sc_obj$treatment <- factor(sc_obj$treatment,
                           levels = intersect(c("Control", "BRZ", "BL"), unique(meta_data$treatment)))
if ("time_zone" %in% colnames(meta_data)) {
  sc_obj$time_zone <- factor(sc_obj$time_zone,
                              levels = intersect(zone_levels, unique(meta_data$time_zone)))
}
if ("cell_type" %in% colnames(meta_data)) {
  sc_obj$cell_type <- factor(sc_obj$cell_type,
                              levels = intersect(celltype_levels, unique(meta_data$cell_type)))
}

# ── 4. Subsample cells to keep memory manageable ─────────────────────────────
# CalcSensors coerces sparse→dense (genes × cells), which requires ~18 GB for
# 80K cells. Subsample to ≤3000 cells per condition (stratified), keeping all
# conditions equally represented. Total ≈ 21K cells → dense mat ≈ 4.7 GB.
set.seed(42)
max_per_cond <- 3000
sub_cells <- sc_obj@meta.data %>%
  rownames_to_column("cell") %>%
  group_by(condition) %>%
  slice_sample(n = max_per_cond) %>%
  pull(cell)
cat("Subsampling to", length(sub_cells), "cells (",
    max_per_cond, "max per condition)...\n")
sc_sub <- subset(sc_obj, cells = sub_cells)
rm(sc_obj); invisible(gc())
cat("  Subsampled: condition counts:\n")
print(table(sc_sub$condition))

# ── 5. Load panels: BR (custom) + Auxin trans-panels only ─────────────────────
# Using only trans-panels (not cis/reg) to keep the analysis focused and RAM low
cat("Loading BR custom panels...\n")
BRpanels <- LoadSensors(setName = "BR", defaultPanels = FALSE, customPanels = TRUE)
cat("  BR panels:", paste(names(BRpanels$panels), collapse = ", "), "\n")

cat("Loading Auxin default panels (ATH, trans-panels only)...\n")
auxpanels <- LoadSensors(setName = "Auxin", species = "ATH", hormone = "aux",
                         customPanels = FALSE)
# Keep only trans-panels (ARF, IAA, Synthesis, Transport, PAT, Receptors, etc.)
aux_trans_names <- grep("^ATH-aux-trans", names(auxpanels$panels), value = TRUE)
auxpanels$panels <- auxpanels$panels[aux_trans_names]
cat("  Aux trans panels:", paste(names(auxpanels$panels), collapse = ", "), "\n")

# Combine panel sets: append auxin panels into BR panel object
combined_panels <- BRpanels
combined_panels$panels <- c(BRpanels$panels, auxpanels$panels)
cat("  Total panels:", length(combined_panels$panels), "\n")

# ── 6. Run CalcSensors ───────────────────────────────────────────────────────
cat("Running CalcSensors on", ncol(sc_sub), "cells with",
    length(combined_panels$panels), "panels...\n")
iSensors_obj <- CalcSensors(sc_sub, seurLayer = "data",
                             panelSet = combined_panels, signals = "mean")
cat("  Done. Assays:", paste(Assays(iSensors_obj), collapse = ", "), "\n")

saveRDS(iSensors_obj, "00-BR-objects/out/GSE212230_iSensors_obj.rds")
cat("  Saved: out/GSE212230_iSensors_obj.rds\n")

rm(sc_sub); invisible(gc())

# ── 6. Extract scores ─────────────────────────────────────────────────────────
DefaultAssay(iSensors_obj) <- "iSensors_mean"
all_rows <- rownames(iSensors_obj[["iSensors_mean"]])

br_panel_order <- c(
  "ATH-br-trans-Biosynthesis",
  "ATH-br-trans-Homeostasis",
  "ATH-br-trans-PositiveSignaling",
  "ATH-br-trans-NegativeSignaling",
  "ATH-br-trans-TF",
  "ATH-br-trans-TF-induced",
  "ATH-br-trans-TF-repressed",
  "ATH-br-reg-2hr-up",
  "ATH-br-reg-2hr-down",
  "ATH-br-reg-4hr-up",
  "ATH-br-reg-4hr-down"
)
aux_panel_order <- grep("^ATH-aux-trans", all_rows, value = TRUE)
panel_order     <- intersect(c(br_panel_order, aux_panel_order), all_rows)

br_labels <- c(
  "ATH-br-trans-Biosynthesis"      = "BR Biosynthesis (inv.)",
  "ATH-br-trans-Homeostasis"       = "BR Homeostasis",
  "ATH-br-trans-PositiveSignaling" = "BR Pos. Signaling",
  "ATH-br-trans-NegativeSignaling" = "BR Neg. Signaling",
  "ATH-br-trans-TF"                = "BR TF (all)",
  "ATH-br-trans-TF-induced"        = "BR TF induced",
  "ATH-br-trans-TF-repressed"      = "BR TF repressed (inv.)",
  "ATH-br-reg-2hr-up"              = "BR 2hr DEG up",
  "ATH-br-reg-2hr-down"            = "BR 2hr DEG down",
  "ATH-br-reg-4hr-up"              = "BR 4hr DEG up",
  "ATH-br-reg-4hr-down"            = "BR 4hr DEG down"
)
# Auto-labels for auxin panels (use panel name minus the ATH-aux-trans- prefix)
aux_labels <- setNames(
  sub("ATH-aux-trans-", "Aux ", aux_panel_order),
  aux_panel_order
)
sensor_labels <- c(br_labels, aux_labels)

cat("  Sensors to plot:", length(panel_order), "\n")
cat("  Panels found:", paste(panel_order, collapse = "\n    "), "\n")

scores_mat <- as.data.frame(t(as.matrix(
  GetAssayData(iSensors_obj, assay = "iSensors_mean", layer = "data")[panel_order, ]
)))
scores_mat$cell_barcode <- rownames(scores_mat)
scores_mat$condition    <- as.character(iSensors_obj$condition)
scores_mat$cell_type    <- as.character(iSensors_obj$cell_type)
scores_mat$time_zone    <- as.character(iSensors_obj$time_zone)
scores_mat$sample       <- as.character(iSensors_obj$sample)

scores_long <- scores_mat %>%
  pivot_longer(cols = all_of(panel_order),
               names_to = "sensor", values_to = "score") %>%
  mutate(
    sensor    = factor(sensor, levels = rev(panel_order)),
    condition = factor(condition, levels = condition_order),
    cell_type = factor(cell_type, levels = intersect(celltype_levels, unique(cell_type))),
    time_zone = factor(time_zone, levels = intersect(zone_levels, unique(time_zone)))
  )

rm(iSensors_obj); invisible(gc())

# Condition colour palette
cond_colors <- c(
  "Control"  = "#2980B9",
  "BRZ"      = "#8E44AD",
  "BL_0.5h"  = "#ABEBC6",
  "BL_1h"    = "#52BE80",
  "BL_2h"    = "#1E8449",
  "BL_4h"    = "#E59866",
  "BL_8h"    = "#CA6F1E"
)
cond_colors <- cond_colors[names(cond_colors) %in% condition_order]

# ── 7. Line plot: mean score per sensor × condition ───────────────────────────
cat("Plotting line plots...\n")

cond_means <- scores_long %>%
  group_by(sensor, condition) %>%
  summarise(mean_score = mean(score, na.rm = TRUE),
            se_score   = sd(score, na.rm = TRUE) / sqrt(n()),
            .groups = "drop") %>%
  mutate(
    sensor_label = sensor_labels[as.character(sensor)],
    sensor_label = ifelse(is.na(sensor_label), as.character(sensor), sensor_label),
    sensor_label = factor(sensor_label,
                          levels = rev(sensor_labels[levels(sensor)]))
  )

# Determine which group each sensor belongs to (for faceting)
cond_means <- cond_means %>%
  mutate(group = ifelse(grepl("^ATH-br", as.character(sensor)), "BR sensors", "Auxin sensors"))

p_line <- ggplot(cond_means,
                 aes(x = condition, y = mean_score, colour = condition, group = sensor_label)) +
  geom_line(colour = "grey60", linewidth = 0.5) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - se_score, ymax = mean_score + se_score),
                width = 0.25, linewidth = 0.5) +
  facet_grid(sensor_label ~ group, scales = "free", space = "free_x") +
  scale_colour_manual(values = cond_colors) +
  labs(x = NULL, y = "Mean iSensor score",
       title = "BR + Auxin iSensors — GSE212230 BL time course") +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1),
    strip.text.y    = element_text(size = 7, angle = 0, hjust = 0),
    strip.text.x    = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

pdf("out/GSE212230_line_plots.pdf", width = 10, height = 14)
print(p_line)
dev.off()
cat("  Saved: out/GSE212230_line_plots.pdf\n")

# ── 8. Heatmap: sensors × condition, faceted by cell type ─────────────────────
cat("Plotting cell type heatmap...\n")

ct_means_raw <- scores_long %>%
  filter(!is.na(cell_type)) %>%
  group_by(sensor, condition, cell_type) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

write.csv(ct_means_raw, "out/GSE212230_celltype_sensor_means.csv", row.names = FALSE)
cat("  Saved: out/GSE212230_celltype_sensor_means.csv\n")

ct_means <- ct_means_raw %>%
  group_by(sensor) %>%
  mutate(z_score = (mean_score - mean(mean_score, na.rm = TRUE)) /
                     (sd(mean_score, na.rm = TRUE) + 1e-9)) %>%
  ungroup() %>%
  mutate(
    sensor_label = sensor_labels[as.character(sensor)],
    sensor_label = ifelse(is.na(sensor_label), as.character(sensor), sensor_label),
    sensor_label = factor(sensor_label,
                          levels = sensor_labels[levels(sensor)])
  )

p_ct_heatmap <- ggplot(ct_means, aes(x = condition, y = sensor_label, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(. ~ cell_type, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Z-score") +
  labs(x = NULL, y = NULL,
       title = "iSensor activity by cell type — GSE212230 BL time course") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    strip.text      = element_text(size = 6, angle = 90, hjust = 0),
    panel.spacing   = unit(0.1, "cm"),
    legend.position = "right"
  )

pdf("out/GSE212230_celltype_heatmap.pdf", width = 16, height = 7)
print(p_ct_heatmap)
dev.off()
cat("  Saved: out/GSE212230_celltype_heatmap.pdf\n")

# ── 9. Heatmap: sensors × condition, faceted by spatial zone ──────────────────
cat("Plotting spatial zone heatmap...\n")

zone_means_raw <- scores_long %>%
  filter(!is.na(time_zone)) %>%
  group_by(sensor, condition, time_zone) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

write.csv(zone_means_raw, "out/GSE212230_zone_sensor_means.csv", row.names = FALSE)
cat("  Saved: out/GSE212230_zone_sensor_means.csv\n")

zone_means <- zone_means_raw %>%
  group_by(sensor) %>%
  mutate(z_score = (mean_score - mean(mean_score, na.rm = TRUE)) /
                     (sd(mean_score, na.rm = TRUE) + 1e-9)) %>%
  ungroup() %>%
  mutate(
    zone_label = zone_labels[as.character(time_zone)],
    zone_label = ifelse(is.na(zone_label), as.character(time_zone), zone_label),
    zone_label = factor(zone_label, levels = zone_labels[zone_levels]),
    sensor_label = sensor_labels[as.character(sensor)],
    sensor_label = ifelse(is.na(sensor_label), as.character(sensor), sensor_label),
    sensor_label = factor(sensor_label,
                          levels = sensor_labels[levels(sensor)])
  )

p_zone_heatmap <- ggplot(zone_means, aes(x = condition, y = sensor_label, fill = z_score)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  facet_grid(. ~ zone_label, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "Z-score") +
  labs(x = NULL, y = NULL,
       title = "iSensor activity by spatial zone — GSE212230 BL time course") +
  theme_bw(base_size = 8) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y     = element_text(size = 7),
    strip.text      = element_text(size = 7, angle = 90, hjust = 0),
    panel.spacing   = unit(0.1, "cm"),
    legend.position = "right"
  )

pdf("out/GSE212230_zone_heatmap.pdf", width = 14, height = 7)
print(p_zone_heatmap)
dev.off()
cat("  Saved: out/GSE212230_zone_heatmap.pdf\n")

# ── 10. Save per-condition mean scores to CSV ─────────────────────────────────
cond_means %>%
  mutate(sensor = as.character(sensor),
         condition = as.character(condition)) %>%
  select(sensor, sensor_label, group, condition, mean_score, se_score) %>%
  write.csv("out/GSE212230_condition_means.csv", row.names = FALSE)
cat("  Saved: out/GSE212230_condition_means.csv\n")

cat("\nAll done!\n")
