library(Seurat)
library(iSensors)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(stringr)
library(broom)
library(ggrepel)

set.seed(100)
iSensors_obj <- readRDS(file = "C:/!Victoria/!Projects/DigitalSensors/Data/Martin-Arevalillo-2025/iSensors-Martin-Arevalillo-auxin-root2025.rds")
iSensors_obj <- readRDS(file = "02-single-cell-data-analysis/in/iSensors-Martin-Arevalillo-auxin-root2025.rds")

DefaultAssay(iSensors_obj) <-"iSensors_mean_normed"

# exploring the global trend based on replicas----

# Pull sample IDs from meta.data (guaranteed)
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

meta
table(meta$sample_id, meta$condition, useNA = "ifany")
table(meta$condition, useNA = "ifany")
isensors_list <- unlist(iSensors_obj@assays$iSensors_mean_normed@counts@Dimnames[1])


# Pull sensor values (works whether they are in meta.data or an assay, as long as FetchData can see them)
sens <- FetchData(iSensors_obj, vars = isensors_list) %>%
  tibble::rownames_to_column("cell")

df <- meta %>%
  inner_join(sens, by = "cell")


sensor_meta <- tibble(iSensor = isensors_list) %>%
  mutate(
    sensor_class = case_when(
      str_detect(iSensor, "^ATH-aux-reg-") ~ "reg",
      str_detect(iSensor, "^ATH-aux-cis-")      ~ "cis",
      str_detect(iSensor, "^ATH-aux-trans-")    ~ "trans",
      TRUE                                     ~ "other"
    ),
    sensor_label = iSensor %>%
      str_remove("^ATH-aux-reg-") %>%
      str_remove("^ATH-aux-cis-") %>%
      str_remove("^ATH-aux-trans-")
  )

sensor_meta$sensor_class <- factor(sensor_meta$sensor_class,
                                   levels = c("cis","trans","reg","other"))

# Ensure consistent condition naming/order
df <- df %>%
  mutate(condition = recode(condition, "Ctr" = "Control", "Aux" = "Auxin")) %>%
  mutate(condition = factor(condition, levels = c("Control", "Auxin")))

# Pseudobulk per sample (replicate-aware)
pb_wide <- df %>%
  group_by(sample_id, condition) %>%
  summarise(
    across(all_of(isensors_list), ~ mean(.x, na.rm = TRUE)),
    n_cells = dplyr::n(),
    .groups = "drop"
  ) %>%
  filter(n_cells > 0)

pb_long <- pb_wide %>%
  pivot_longer(
    cols = all_of(isensors_list),
    names_to = "iSensor",
    values_to = "score"
  )

hm_df <- pb_wide %>%
  select(sample_id, condition, all_of(isensors_list)) %>%
  pivot_longer(
    cols = all_of(isensors_list),
    names_to = "iSensor",
    values_to = "value"
  ) %>%
#  filter(iSensor %in% is_order) %>%
#  mutate(iSensor = factor(iSensor, levels = is_order)) %>%
  group_by(iSensor) %>%
  mutate(value_z = as.numeric(scale(value))) %>%   # z-score per sensor
  ungroup()

# Weighted replicate-aware linear model per iSensor
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
    class = case_when(
      p_adj < 0.05 & effect > 0.05 ~ "Up",
      p_adj < 0.05 & effect < -0.05 ~ "Down",
      TRUE                     ~ "n.s."
    )
  ) %>%
  arrange(p_adj, desc(abs(effect)))



# Choose ordering; effect ordering is intuitive for a forest plot
res_plot <- res_global %>%
  mutate(iSensor = fct_reorder(iSensor, effect))




# Enforce same iSensor order as forest plot
is_order <- levels(res_plot$iSensor)  # from your forest plot ordering

sensor_meta <- sensor_meta %>%
  filter(iSensor %in% is_order) %>%
  mutate(
    iSensor = factor(iSensor, levels = is_order),
    sensor_label = factor(sensor_label, levels = sensor_label[match(is_order, iSensor)])
  )


hm_df <- hm_df %>%
  left_join(sensor_meta %>% select(iSensor, sensor_label, sensor_class),
            by = "iSensor") %>%
  mutate(sensor_label = factor(sensor_label, levels = levels(sensor_meta$sensor_label)))

res_plot <- res_plot %>%
  left_join(sensor_meta %>% select(iSensor, sensor_label, sensor_class),
            by = "iSensor") %>%
  mutate(sensor_label = factor(sensor_label, levels = levels(sensor_meta$sensor_label)))


# Order samples: Control first, then Auxin (within each, alphabetical sample_id)
sample_order <- hm_df %>%
  distinct(sample_id, condition) %>%
  arrange(condition, sample_id) %>%
  pull(sample_id)
hm_df <- hm_df %>%
  mutate(sample_id = factor(sample_id, levels = sample_order))

sample_order
class_cols <- c(
  cis      = "#fdbf6f",
  trans    = "#ff7f00",
  reg = "#b15928",
  control    = "gray"
)

p_class_strip <- ggplot(sensor_meta, aes(x = 1, y = sensor_label, fill = sensor_class)) +
  geom_tile() +
  scale_fill_manual(values = class_cols, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "on") +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0)
  )

p_class_strip


p_forest <- ggplot(res_plot, aes(x = sensor_label, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.18, linewidth = 0.5, color = "grey50") +
  geom_point(aes(color = class), size = 2) +
  scale_color_manual(values = c("Up"="#d73027","Down"="#4575b4","n.s."="grey75"),
                     guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Effect (Auxin – Control)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(4, 6, 4, 2)
  )


p_forest


p_heat <- ggplot(hm_df, aes(x = sample_id, y = sensor_label, fill = value_z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027", midpoint = 0,
                       name = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 8, hjust = 0),
#    axis.text.y = element_blank(),
    axis.line    = element_blank(),   # <<< removes axis lines
    axis.ticks.y = element_blank(),
    legend.position = "bottom",
    legend.key.width = unit(1.2, "cm"),
    legend.key.height = unit(0.3, "cm"),
    plot.margin = margin(2, 4, 4, 2)
  )

p_heat

library(patchwork)

p_heat   <- p_heat   + theme(plot.margin = margin(2, 2, 2, 2))
p_forest <- p_forest + theme(plot.margin = margin(2, 2, 2, 2))


p_final <- (p_heat | p_class_strip | p_forest) +
  plot_layout(widths = c(1.25, 0.05, 0.60))&
  theme(plot.margin = margin(6, 6, 6, 6))

p_final

ggsave(
  "02-single-cell-data-analysis/out/iSensors_replica_heatmap_plus_forest_2.pdf",
  plot = p_final,
  width = 10,
  height = 8,
  dpi = 300,
  useDingbats = FALSE
)


#vizualizing fold change heatmap

condition_col <- "orig.ident2"
celltype_col  <- "Cell.Ident"   # change to "Cell.Ident" if that is what you use

levels(iSensors_obj@meta.data$Cell.Ident)

#cell-type specific heatmap

keep_celltypes <- c(
 "Columella", "Young-Atrichoblast-1", "Young-Trichoblast-1", "Cortex", "Endodermis", "Stele-Pericycle")

iSensors_obj[[condition_col]][,1] <- recode(
  iSensors_obj[[condition_col]][,1],
  "Ctr" = "Control",
  "Aux" = "Auxin"
)

iSensors_obj[[condition_col]][,1] <- factor(
  iSensors_obj[[condition_col]][,1],
  levels = c("Control", "Auxin")
)

assay_use <- "iSensors_mean_normed"  # change if your p_final uses another assay
DefaultAssay(iSensors_obj) <- assay_use

avg_list <- AverageExpression(
  iSensors_obj,
  assays   = assay_use,
  features = is_order,                       # exact same sensors/order as p_final
  group.by = c(condition_col, celltype_col),
  slot     = "data",
  verbose  = FALSE
)


avg_mat <- avg_list[[assay_use]]  # rows=iSensors, cols="Control_<celltype>" etc.
colnames(avg_mat) <- gsub("\\.", "_", colnames(avg_mat))  # normalize separator
colnames(avg_mat)
ctr_cols <- grep("^Control_", colnames(avg_mat), value = TRUE)
aux_cols <- grep("^Auxin_",  colnames(avg_mat), value = TRUE)

stopifnot(length(ctr_cols) > 0, length(aux_cols) > 0)

ctr_mat <- avg_mat[, ctr_cols, drop = FALSE]
aux_mat <- avg_mat[, aux_cols, drop = FALSE]

colnames(ctr_mat) <- sub("^Control_", "", colnames(ctr_mat))
colnames(aux_mat) <- sub("^Auxin_",  "", colnames(aux_mat))

common_ct <- intersect(colnames(ctr_mat), colnames(aux_mat))

# apply your requested subset of cell types, safely
keep_ct <- intersect(keep_celltypes, common_ct)
stopifnot(length(keep_ct) > 0)

ctr_mat <- ctr_mat[, keep_ct, drop = FALSE]
aux_mat <- aux_mat[, keep_ct, drop = FALSE]

#compute fold change

rng <- range(avg_mat, na.rm = TRUE)

use_ratio <- rng[1] >= 0  # crude but effective rule

pseudocount <- 1e-6

if (use_ratio) {
  log2_fc <- log2((aux_mat + pseudocount) / (ctr_mat + pseudocount))
  fc_label <- expression(log[2]*"(Auxin/Control)")
} else {
  # centered scores: use difference; still interpretable as "effect"
  log2_fc <- aux_mat - ctr_mat
  fc_label <- "Auxin – Control"
}

log2_fc <- as.matrix(log2_fc)

# enforce iSensor order exactly as p_final
log2_fc <- log2_fc[is_order, , drop = FALSE]

vals <- as.numeric(log2_fc)
vals <- vals[is.finite(vals)]
lims <- quantile(vals, probs = c(0.01, 0.99), na.rm = TRUE)

hm_long <- as.data.frame(log2_fc) %>%
  tibble::rownames_to_column("iSensor") %>%
  pivot_longer(-iSensor, names_to = "celltype", values_to = "effect") %>%
  mutate(
    iSensor  = factor(iSensor, levels = is_order),
    celltype = factor(celltype, levels = keep_ct)
  )

p_fc_heatmap <- ggplot(hm_long, aes(x = celltype, y = iSensor, fill = effect)) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#008837", mid = "white", high = "#7b3294",
    midpoint = 0,
    limits = lims,
    oob = scales::squish,
    name = NULL
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_blank(),   # <<< remove x-axis labels
#    axis.text.x = element_text(size = 8),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.key.width  = unit(1.5, "cm"),
    legend.key.height = unit(0.35, "cm"),
    plot.margin = margin(4, 6, 4, 6)
  )


p_fc_heatmap


p_final2 <- (p_fc_heatmap | p_heat | p_class_strip | p_forest) +
  plot_layout(widths = c(1, 1.25, 0.05, 0.60))&
  theme(plot.margin = margin(0, 6, 6, 6))


p_final2
ggsave(
  "02-single-cell-data-analysis/out/iSensors_replica_heatmap_plus_forest.pdf",
  plot = p_final2,
  width = 14,
  height = 8,
  dpi = 300,
  useDingbats = FALSE
)
