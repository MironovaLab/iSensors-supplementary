# 03-iSensors-AuxinRoot-metasensor-evaluation.R
# Adds 2 ARF/IAA meta-sensors (product, sum) to the Martin-Arevalillo
# GSE241573 root Auxin dataset and reproduces the replica heatmap + forest
# + cell-type fold-change figure from 02-iSensors-AuxinRoot-global-evaluation.R
# (iSensors_replica_heatmap_plus_forest.pdf), with the meta-sensors included
# in the same ranking/plot so their performance is directly comparable to
# ARF, IAA, and all other Auxin panels.
#
# Meta-panel mechanism (iSensors::LoadSensors(metaPanels=...) + CalcSensors):
# per-cell score = rule(ARF_score, IAA_score), applied after normal scoring.
#   ATH-aux-meta-ARF-IAA-product = ARF_score(cell) * IAA_score(cell)
#   ATH-aux-meta-ARF-IAA-sum     = ARF_score(cell) + IAA_score(cell)
#
# One fix vs the original reference script (verified against the raw
# object's actual metadata before writing this):
#   keep_celltypes used hyphens ("Young-Atrichoblast-1") but Cell.Ident
#   actually uses underscores ("Young_Atrichoblast_1") - only 3 of 6
#   requested cell types survived the original intersect(). Fixed to
#   match actual Cell.Ident values.
#
# Condition regex matches the original reference script exactly: only
# samples ending in exactly "_T"/"_C" are kept (DR15_T/_C, DR5_T/_C,
# ER8_T/_C, IR13_T/_C = 8 samples); the "_2"/"_3" replicate-suffixed
# samples (DR5_C_2/_3, DR5_T_2/_3, IR13_C_2, IR13_T_2) are intentionally
# dropped, per request.
#
# Input:  in/MartinArevalillo2025.rds (raw, unscored Seurat object)
# Output: out/iSensors_replica_heatmap_plus_forest_with_ARFIAA_meta.pdf
#         out/metasensor_evaluation_stats.csv

library(Seurat)
library(iSensors)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(stringr)
library(broom)
library(patchwork)

set.seed(100)
setwd("D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/02-single-cell-data-analysis")
dir.create("out", showWarnings = FALSE)

# ── 1. Load raw object, define meta-panels, score ────────────────────────────
# (Reuses the already-scored object from the prior run if present - the
# meta-panel definitions/scoring don't depend on the sample regex below, so
# no need to recompute the ~20+ min CalcSensors step just to change which
# samples are kept for the pseudobulk lm.)
meta_panels_def <- list(
  "ATH-aux-meta-ARF-IAA-product" = list(
    srcPanels = c("ATH-aux-trans-ARF", "ATH-aux-trans-IAA"),
    rule = prod
  ),
  "ATH-aux-meta-ARF-IAA-sum" = list(
    srcPanels = c("ATH-aux-trans-ARF", "ATH-aux-trans-IAA"),
    rule = sum
  )
)
meta_names <- names(meta_panels_def)

scored_cache <- "in/iSensors-Martin-Arevalillo-auxin-root2025_with_ARFIAA_meta.rds"
if (file.exists(scored_cache)) {
  cat("Loading cached scored object...\n")
  iSensors_obj <- readRDS(scored_cache)
} else {
  cat("Loading raw Martin-Arevalillo object...\n")
  seurat_obj <- readRDS("in/MartinArevalillo2025.rds")
  cat("  Cells:", ncol(seurat_obj), " Genes:", nrow(seurat_obj), "\n")

  cat("Loading all Auxin panels + ARF/IAA meta-panels...\n")
  customTransPanel <- LoadSensors(
    setName = 'Auxin', species = 'ATH', hormone = 'aux',
    customPanels = FALSE,
    randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), majortrend = TRUE),
    metaPanels = meta_panels_def
  )
  cat("  Panels loaded (excluding meta):", length(customTransPanel$panels), "\n")

  cat("Calculating sensor scores (mean only - mean_normed is not used)...\n")
  iSensors_obj <- CalcSensors(
    seurat_obj, seurLayer = 'data', panelSet = customTransPanel,
    signals = "mean"
  )
  rm(seurat_obj); invisible(gc())
  saveRDS(iSensors_obj, scored_cache)
}

# NOTE: signals="mean" only - mean_normed is not used for this analysis.
DefaultAssay(iSensors_obj) <- "iSensors_mean"
isensors_list <- unlist(iSensors_obj@assays$iSensors_mean@counts@Dimnames[1])
cat("  Total sensors (incl. meta + random/majortrend controls):", length(isensors_list), "\n")
cat("  Meta panels present:", paste(intersect(meta_names, isensors_list), collapse = ", "), "\n\n")

# ── 2. Sample metadata: sample_id/condition (original regex - exact _T/_C only,
#    "_2"/"_3" replicate suffixes intentionally dropped per request) ─────────
meta <- iSensors_obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  transmute(
    cell,
    sample_id = orig.ident,
    condition = case_when(
      str_ends(sample_id, "_T") ~ "Auxin",
      str_ends(sample_id, "_C") ~ "Control",
      TRUE ~ NA_character_
    ),
    condition = factor(condition, levels = c("Control", "Auxin"))
  ) %>%
  filter(!is.na(condition))

cat("Samples included (original regex - exact _T/_C only):\n")
print(table(meta$sample_id, meta$condition, useNA = "ifany"))
cat("\nTotal samples:", n_distinct(meta$sample_id), "(was 8/14 with the original regex)\n\n")

sens <- FetchData(iSensors_obj, vars = isensors_list) %>%
  tibble::rownames_to_column("cell")
df <- meta %>% inner_join(sens, by = "cell")

sensor_meta <- tibble(iSensor = isensors_list) %>%
  mutate(
    sensor_class = case_when(
      iSensor %in% meta_names               ~ "meta",
      str_detect(iSensor, "^ATH-aux-reg-")   ~ "reg",
      str_detect(iSensor, "^ATH-aux-cis-")   ~ "cis",
      str_detect(iSensor, "^ATH-aux-trans-") ~ "trans",
      TRUE                                    ~ "other"
    ),
    sensor_label = iSensor %>%
      str_remove("^ATH-aux-meta-") %>%
      str_remove("^ATH-aux-reg-") %>%
      str_remove("^ATH-aux-cis-") %>%
      str_remove("^ATH-aux-trans-")
  )
sensor_meta$sensor_class <- factor(sensor_meta$sensor_class,
                                   levels = c("cis", "trans", "reg", "meta", "other"))

# ── 3. Pseudobulk per sample (replicate-aware), weighted lm per iSensor ──────
pb_wide <- df %>%
  group_by(sample_id, condition) %>%
  summarise(across(all_of(isensors_list), ~ mean(.x, na.rm = TRUE)),
            n_cells = dplyr::n(), .groups = "drop") %>%
  filter(n_cells > 0)

pb_long <- pb_wide %>%
  pivot_longer(cols = all_of(isensors_list), names_to = "iSensor", values_to = "score")

hm_df <- pb_wide %>%
  select(sample_id, condition, all_of(isensors_list)) %>%
  pivot_longer(cols = all_of(isensors_list), names_to = "iSensor", values_to = "value") %>%
  group_by(iSensor) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup()

cat("Fitting weighted replicate-aware lm per iSensor (score ~ condition, weights=n_cells)...\n")
res_global <- pb_long %>%
  group_by(iSensor) %>%
  group_modify(~{
    fit <- lm(score ~ condition, data = .x, weights = n_cells)
    broom::tidy(fit) %>% filter(term == "conditionAuxin")
  }) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    effect = estimate,
    lo = effect - 1.96 * std.error,
    hi = effect + 1.96 * std.error,
    # Classify by significance alone (p_adj<0.05); the original reference
    # script's additional |effect|>0.05 floor was calibrated for the
    # mean_normed scale and mislabels real effects as "n.s." on the much
    # smaller raw "mean" scale (e.g. ARF: effect=0.047, p_adj=0.003).
    class = case_when(
      p_adj < 0.05 & effect > 0 ~ "Up",
      p_adj < 0.05 & effect < 0 ~ "Down",
      TRUE ~ "n.s."
    )
  ) %>%
  arrange(p_adj, desc(abs(effect)))

write.csv(res_global, "out/metasensor_evaluation_stats.csv", row.names = FALSE)
cat("Saved: out/metasensor_evaluation_stats.csv\n\n")

cat("=== ARF, IAA, and meta-panel results ===\n")
print(res_global %>% filter(iSensor %in% c("ATH-aux-trans-ARF", "ATH-aux-trans-IAA", meta_names)) %>%
        select(iSensor, effect, lo, hi, p_adj, class) %>% as.data.frame())
cat("\n")

# ── 4. Plot layout (matches original script exactly) ─────────────────────────
res_plot <- res_global %>% mutate(iSensor = fct_reorder(iSensor, effect))
is_order <- levels(res_plot$iSensor)

sensor_meta <- sensor_meta %>%
  filter(iSensor %in% is_order) %>%
  mutate(iSensor = factor(iSensor, levels = is_order),
         sensor_label = factor(sensor_label, levels = sensor_label[match(is_order, iSensor)]))

hm_df <- hm_df %>%
  left_join(sensor_meta %>% select(iSensor, sensor_label, sensor_class), by = "iSensor") %>%
  mutate(sensor_label = factor(sensor_label, levels = levels(sensor_meta$sensor_label)))

res_plot <- res_plot %>%
  left_join(sensor_meta %>% select(iSensor, sensor_label, sensor_class), by = "iSensor") %>%
  mutate(sensor_label = factor(sensor_label, levels = levels(sensor_meta$sensor_label)))

sample_order <- hm_df %>% distinct(sample_id, condition) %>%
  arrange(condition, sample_id) %>% pull(sample_id)
hm_df <- hm_df %>% mutate(sample_id = factor(sample_id, levels = sample_order))

class_cols <- c(cis = "#fdbf6f", trans = "#ff7f00", reg = "#b15928",
                meta = "#1f78b4", other = "gray", control = "gray")

p_class_strip <- ggplot(sensor_meta, aes(x = 1, y = sensor_label, fill = sensor_class)) +
  geom_tile() +
  scale_fill_manual(values = class_cols, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "on") + theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

p_forest <- ggplot(res_plot, aes(x = sensor_label, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.5, color = "grey50") +
  geom_point(aes(color = class, shape = sensor_class == "meta"), size = 2) +
  scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "n.s." = "grey75"), guide = "none") +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16), guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Effect (Auxin - Control)") +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin(4, 6, 4, 2))

p_heat <- ggplot(hm_df, aes(x = sample_id, y = sensor_label, fill = value_z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0, name = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6, hjust = 0), axis.line = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.3, "cm"),
        plot.margin = margin(2, 4, 4, 2))

# ── 5. Cell-type fold-change heatmap (FIXED cell-type names) ─────────────────
condition_col <- "orig.ident2"
celltype_col  <- "Cell.Ident"

iSensors_obj[[condition_col]][, 1] <- recode(iSensors_obj[[condition_col]][, 1], "Ctr" = "Control", "Aux" = "Auxin")
iSensors_obj[[condition_col]][, 1] <- factor(iSensors_obj[[condition_col]][, 1], levels = c("Control", "Auxin"))

keep_celltypes <- c("Columella", "Young_Atrichoblast_1", "Young_Trichoblast_1",
                     "Cortex", "Endodermis", "Stele_Pericycle")

assay_use <- "iSensors_mean"
DefaultAssay(iSensors_obj) <- assay_use

cat("Computing per-cell-type average expression...\n")
avg_list <- AverageExpression(iSensors_obj, assays = assay_use, features = is_order,
                              group.by = c(condition_col, celltype_col), slot = "data", verbose = FALSE)
avg_mat <- avg_list[[assay_use]]
colnames(avg_mat) <- gsub("\\.", "_", colnames(avg_mat))
ctr_cols <- grep("^Control_", colnames(avg_mat), value = TRUE)
aux_cols <- grep("^Auxin_",  colnames(avg_mat), value = TRUE)
ctr_mat <- avg_mat[, ctr_cols, drop = FALSE]; colnames(ctr_mat) <- sub("^Control_", "", colnames(ctr_mat))
aux_mat <- avg_mat[, aux_cols, drop = FALSE]; colnames(aux_mat) <- sub("^Auxin_",  "", colnames(aux_mat))
common_ct <- intersect(colnames(ctr_mat), colnames(aux_mat))
keep_ct <- intersect(keep_celltypes, common_ct)
cat("  Cell types matched:", paste(keep_ct, collapse = ", "), "(", length(keep_ct), "/", length(keep_celltypes), "requested)\n\n")
stopifnot(length(keep_ct) > 0)

ctr_mat <- ctr_mat[, keep_ct, drop = FALSE]
aux_mat <- aux_mat[, keep_ct, drop = FALSE]

rng <- range(avg_mat, na.rm = TRUE)
use_ratio <- rng[1] >= 0
pseudocount <- 1e-6
if (use_ratio) {
  log2_fc <- log2((aux_mat + pseudocount) / (ctr_mat + pseudocount))
} else {
  log2_fc <- aux_mat - ctr_mat
}
log2_fc <- as.matrix(log2_fc)[is_order, , drop = FALSE]

vals <- as.numeric(log2_fc); vals <- vals[is.finite(vals)]
lims <- quantile(vals, probs = c(0.01, 0.99), na.rm = TRUE)

hm_long <- as.data.frame(log2_fc) %>%
  tibble::rownames_to_column("iSensor") %>%
  pivot_longer(-iSensor, names_to = "celltype", values_to = "effect") %>%
  mutate(iSensor = factor(iSensor, levels = is_order), celltype = factor(celltype, levels = keep_ct))

p_fc_heatmap <- ggplot(hm_long, aes(x = celltype, y = iSensor, fill = effect)) +
  geom_tile() +
  scale_fill_gradient2(low = "#008837", mid = "white", high = "#7b3294", midpoint = 0,
                       limits = lims, oob = scales::squish, name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), panel.grid = element_blank(),
        axis.text.y = element_blank(), axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
        legend.position = "bottom", legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.35, "cm"),
        plot.margin = margin(4, 6, 4, 6))

# ── 6. Compose and save ───────────────────────────────────────────────────────
p_heat   <- p_heat   + theme(plot.margin = margin(2, 2, 2, 2))
p_forest <- p_forest + theme(plot.margin = margin(2, 2, 2, 2))

p_final <- (p_fc_heatmap | p_heat | p_class_strip | p_forest) +
  plot_layout(widths = c(1, 1.25, 0.05, 0.60)) &
  theme(plot.margin = margin(0, 6, 6, 6))

out_file <- "out/iSensors_replica_heatmap_plus_forest_with_ARFIAA_meta.pdf"
ggsave(out_file, plot = p_final, width = 15, height = 9, dpi = 300, useDingbats = FALSE)
cat("Saved:", out_file, "\n")

png_file <- sub("\\.pdf$", ".png", out_file)
ggsave(png_file, plot = p_final, width = 15, height = 9, dpi = 200)
cat("Saved:", png_file, "\n")
