library(Seurat)
library(tidyverse)
library(limma)
library(dplyr)
library(tibble)
library(stringr)
library(statmod)
library(iSensors)
set.seed(100)

#install.packages("statmod")
in_path <- "01-bulk-data-obtain-and-analyze/Data/"
out_path <- "01-bulk-data-obtain-and-analyze/Out/"
pseudocells <- readRDS(paste0(in_path, "seurat_auxin_microarrays.rds"))

DimPlot(pseudocells)
meta_pseudo <- pseudocells@meta.data

head(meta_pseudo)

#Loading custom panels
customTransPanel <- LoadSensors(setName = 'Auxin', species = 'AT', hormone = 'aux',
                                customPanels = FALSE,
                                randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), 
                                                  majortrend = TRUE),
                                metaPanels = list('Response1' = list('srcPanels' = c('AT_aux_trans_A_ARF', 'AT_aux_trans_IAA'),
                                                                     'rule' = prod),
                                                  'Response2' = list('srcPanels' = c('AT_aux_trans_ARF', 'AT_aux_trans_IAA'),
                                                                     'rule' = prod),
                                                  'Response3' = list('srcPanels' = c('AT_aux_trans_A_ARF', 'AT_aux_trans_IAA'),
                                                                     'rule' = mean),
                                                  'Response4' = list('srcPanels' = c('AT_aux_trans_ARF', 'AT_aux_trans_IAA'),
                                                                     'rule' = mean))
)

iSensors_obj2 <- CalcSensors(
  pseudocells,
  seurLayer = 'data',
  panelSet = customTransPanel,
  signals = c("mean_normed")
)

DimPlot(
  iSensors_obj2,
  reduction = "umap",
  split.by = "orig.ident",
  combine = TRUE
)

DefaultAssay(iSensors_obj2) <-"iSensors_mean_normed"

# evaluate iSensors

# --- Extract metadata ---
meta <- iSensors_obj2@meta.data %>%
  rownames_to_column("cell") %>%
  transmute(
    cell,
    experiment_ID = as.factor(Author),
    condition_raw = as.character(orig.ident)
  )

# --- Map condition  ---
meta <- meta %>%
  mutate(
    condition = case_when(
      condition_raw %in% c("Ctr", "control", "C") ~ "Ctr",
      condition_raw %in% c("Aux", "Treatment", "treat") ~ "Aux",
      TRUE ~ NA_character_
    ),
    condition = factor(condition, levels = c("Ctr", "Aux"))
  ) %>%
  filter(!is.na(condition), !is.na(experiment_ID))

# --- Extract iSensors matrix (cells x sensors) ---
isensors_list <- unlist(iSensors_obj2@assays$iSensors_mean_normed@counts@Dimnames[1])

sens <- FetchData(iSensors_obj2, vars = isensors_list) %>%
  rownames_to_column("cell")

df <- meta %>% inner_join(sens, by = "cell")

# Sanity checks
stopifnot(nrow(df) == ncol(iSensors_obj2))  # if all pseudo-cells are included; remove if you subsetted
table(df$condition)
table(df$experiment_ID)

#prpaper matrix for limma
Y <- df %>%
  select(all_of(isensors_list)) %>%
  as.matrix()

# Y is samples x sensors; transpose to sensors x samples
E <- t(Y)

# Ensure column names match samples
colnames(E) <- df$cell

#Fit limma
design <- model.matrix(~ condition, data = df)
# This encodes an intercept + conditionAux effect (Aux vs Ctr)

# Estimate within-experiment correlation
corfit <- duplicateCorrelation(E, design = design, block = df$experiment_ID)

# Fit with blocking and correlation
fit <- lmFit(E, design = design, block = df$experiment_ID, correlation = corfit$consensus)
fit <- eBayes(fit)

#Extract results

res_limma <- topTable(
  fit,
  coef = "conditionAux",
  number = Inf,
  sort.by = "P"
) %>%
  rownames_to_column("iSensor") %>%
  as_tibble() %>%
  rename(
    effect = logFC,      # Aux – Ctr on your iSensor scale
    p_value = P.Value,
    p_adj = adj.P.Val,
    t = t
  ) %>%
  mutate(neglog10_fdr = -log10(p_adj + 1e-300))

res_limma

write.csv(res_limma, file = paste0(out_path, "04-bulk-limma-results.csv"))

# per experiment concordance test
pb_exp <- df %>%
  select(cell, experiment_ID, condition, all_of(isensors_list)) %>%
  tidyr::pivot_longer(cols = all_of(isensors_list), names_to = "iSensor", values_to = "score") %>%
  group_by(experiment_ID, condition, iSensor) %>%
  summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = mean_score) %>%
  mutate(delta = Aux - Ctr)

consistency <- pb_exp %>%
  group_by(iSensor) %>%
  summarise(
    n_experiments = sum(!is.na(delta)),
    frac_positive = mean(delta > 0, na.rm = TRUE),
    sd_delta = sd(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    .groups = "drop"
  )

# Join to limma results
res_final <- res_limma %>%
  left_join(consistency, by = "iSensor")

res_final

write.csv(res_final, file = paste0(out_path, "04-bulk-limma-results-final.csv"))


#Vizaulization on volcano plots (not impressive)
library(ggplot2)
library(ggrepel)

top <- res_final %>% slice_min(p_adj, n = 60)

ggplot(res_final, aes(effect, neglog10_fdr)) +
  geom_point(aes(color = p_adj < 0.05), alpha = 0.8, size = 1.6) +
  scale_color_manual(values = c(`TRUE` = "#D73027", `FALSE` = "grey70"),
                     labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."),
                     name = NULL) +
  geom_text_repel(data = top, aes(label = iSensor), size = 3, max.overlaps = 50) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic()


library(dplyr)

res_forest <- res_limma %>%
  mutate(
    se = abs(effect / t),
    lo = effect - 1.96 * se,
    hi = effect + 1.96 * se,
    sig = p_adj < 0.05
  )


ggplot(res_forest, aes(y = iSensor, x = effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(
    aes(xmin = lo, xmax = hi),
    height = 0.25,
    color = "grey50"
  ) +
  geom_point(
    aes(color = sig),
    size = 2.5
  ) +
  scale_color_manual(
    values = c(`TRUE` = "#D73027", `FALSE` = "grey70"),
    labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "n.s."),
    name = NULL
  ) +
  labs(
    x = "Effect size (Aux – Control)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9)
  )
ggsave(plot = get_last_plot(), filename = paste0(out_path, "iSensors_microarray_effect_size.jpg"), dpi = 300)

