library(Seurat)
library(iSensors)
library(tidyverse)
library(RColorBrewer)
library(scales)

iSensors_obj <- readRDS(file = "02-single-cell-data-analysis/in/iSensors-Martin-Arevalillo-auxin-root2025_mean.rds")


#For supplementary Figure 4A
p <- DimPlot(
  iSensors_obj,
  reduction = "umap",
  split.by = "orig.ident2",
  combine = TRUE,
  #  label = TRUE,
  label.size = 2.2,        # smaller cluster labels
  repel = TRUE
) +
  scale_color_manual(
    values = Seurat::DiscretePalette(
      length(unique(Idents(iSensors_obj))),
      palette = "polychrome"   # colorblind-friendly, high contrast
    ) 
  )+  guides(color = guide_legend(override.aes = list(size = 2)))

p
p+theme(
  legend.text  = element_text(size = 8),   # smaller cluster names
  axis.line        = element_blank(),
  axis.ticks       = element_blank(),
  axis.text        = element_blank(),
  axis.title       = element_blank(),
  legend.key.height = unit(0.4, "cm"),  # 🔑 reduces vertical spacing
  
)
ggsave("Manuscript-Figures/out/Supplementary_AuxinTreated_DimPlot_splitted.pdf", last_plot(), width = 8, height = 4.0, dpi = 300)


# ── Figure S3: log2(AUX/CTR) per tissue, non-significant values masked ────────
# Source: iSensors-supplementary/02-single-cell-data-analysis/02-SupplementaryFigure2.R
# Uses the same iSensors_obj loaded at the top of this file (iSensors_mean assay).

library(pheatmap)

assay_use    <- "iSensors_mean"
de_csv_path  <- "iSensors-supplementary/02-single-cell-data-analysis/in/iSensors_DE_aux_vs_ctr_by_tissue.csv"
isensors_list <- rownames(iSensors_obj[[assay_use]])

# ── Step 1: Differential expression per tissue (Wilcoxon, per cell type) ──────
# Re-use the pre-computed CSV if it exists; otherwise run FindMarkers.
if (file.exists(de_csv_path)) {
  message("Reading pre-computed DE table from: ", de_csv_path)
  de_all <- read.csv(de_csv_path, stringsAsFactors = FALSE)
} else {
  message("DE CSV not found — running FindMarkers per tissue (this may take a while) …")

  DefaultAssay(iSensors_obj) <- assay_use
  Idents(iSensors_obj)       <- iSensors_obj$Cell.Ident

  de_list <- list()
  for (tissue in levels(Idents(iSensors_obj))) {
    cells_t <- WhichCells(iSensors_obj, idents = tissue)
    if (length(cells_t) < 20) next

    obj_t <- subset(iSensors_obj, cells = cells_t)
    conds <- unique(as.character(obj_t$orig.ident2))
    if (!all(c("Control", "Auxin") %in% conds)) next

    de <- FindMarkers(
      object          = obj_t,
      ident.1         = "Auxin",
      ident.2         = "Control",
      group.by        = "orig.ident2",
      features        = isensors_list,
      test.use        = "wilcox",
      logfc.threshold = 0,
      min.pct         = 0
    )
    if (nrow(de) == 0) next

    de_list[[tissue]] <- de %>%
      rownames_to_column("iSensor") %>%
      mutate(tissue     = tissue,
             p_val_adj  = p.adjust(p_val, method = "BH"))
  }

  de_all <- bind_rows(de_list)
  write.csv(de_all, de_csv_path, row.names = FALSE)
  message("DE table saved to: ", de_csv_path)
}

# ── Step 2: Average iSensor signal per condition × tissue ────────────────────
DefaultAssay(iSensors_obj) <- assay_use

avg_mat <- AverageExpression(
  iSensors_obj,
  assays   = assay_use,
  features = isensors_list,
  group.by = c("orig.ident2", "Cell.Ident"),
  slot     = "data",
  verbose  = FALSE
)[[assay_use]]

ctr_cols <- grep("^Control_", colnames(avg_mat), value = TRUE)
aux_cols <- grep("^Auxin_",   colnames(avg_mat), value = TRUE)

ctr_mat  <- avg_mat[, ctr_cols, drop = FALSE]
aux_mat  <- avg_mat[, aux_cols, drop = FALSE]
colnames(ctr_mat) <- sub("^Control_", "", colnames(ctr_mat))
colnames(aux_mat) <- sub("^Auxin_",   "", colnames(aux_mat))

common_tissues <- intersect(colnames(ctr_mat), colnames(aux_mat))
ctr_mat <- ctr_mat[, common_tissues, drop = FALSE]
aux_mat <- aux_mat[, common_tissues, drop = FALSE]

# ── Step 3: log2(AUX/CTR) ratio ──────────────────────────────────────────────
pseudocount  <- 1e-6
log2_ratio_m <- log2((aux_mat + pseudocount) / (ctr_mat + pseudocount))
log2_ratio_m <- as.matrix(log2_ratio_m)

# Canonical name conversion so DE table and heatmap rows/cols align
canon <- function(x) gsub("-", "_", trimws(x))
colnames(log2_ratio_m) <- canon(colnames(log2_ratio_m))
de_all$tissue          <- canon(de_all$tissue)

stopifnot(all(colnames(log2_ratio_m) %in% unique(de_all$tissue)))
stopifnot(all(rownames(log2_ratio_m) %in% unique(de_all$iSensor)))

# ── Step 4: Significance mask ─────────────────────────────────────────────────
padj_df  <- de_all %>%
  select(iSensor, tissue, p_val_adj) %>%
  distinct() %>%
  pivot_wider(names_from = tissue, values_from = p_val_adj)

padj_mat              <- as.data.frame(padj_df)
rownames(padj_mat)    <- padj_mat$iSensor
padj_mat$iSensor      <- NULL
padj_mat              <- as.matrix(padj_mat)
padj_mat              <- padj_mat[rownames(log2_ratio_m),
                                   colnames(log2_ratio_m), drop = FALSE]

sig_mat <- padj_mat < 0.05
sig_mat[is.na(sig_mat)] <- FALSE

log2_ratio_masked          <- log2_ratio_m
log2_ratio_masked[!sig_mat] <- NA   # grey out non-significant cells

# ── Step 5: Draw and save the heatmap ─────────────────────────────────────────
vals_m        <- as.numeric(log2_ratio_masked)
vals_m        <- vals_m[is.finite(vals_m)]
lims_m        <- quantile(vals_m, probs = c(0.01, 0.99), na.rm = TRUE)
breaks_masked <- seq(lims_m[1], lims_m[2], length.out = 101)

pheatmap(
  log2_ratio_masked,
  breaks       = breaks_masked,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  main         = "log2(AUX/CTR) per tissue (non-significant masked)",
  border_color = NA,
  na_col       = "grey90",
  filename     = "Manuscript-Figures/out/FigureS3_log2_AUX_vs_CTR_masked.pdf",
  width        = 10,
  height       = 12
)
message("Saved Figure S3: Manuscript-Figures/out/FigureS3_log2_AUX_vs_CTR_masked.pdf")


# ── Supplementary Figure 6: iSensor heatmap across the Guo 2025 multi-organ atlas
# Sensors × cell types, one panel per selected sensor, organs on Y axis.
# Input: 00-iSensors-objects/data/Guo/guo_avg_exp_iSensors_mean.rds
#        (generated by 00-iSensors-objects/guo-atlas-iSensors-objects.R)

library(patchwork)

guo_avg_exp_s6 <- readRDS("00-iSensors-objects/data/Guo/guo_avg_exp_iSensors_mean.rds")

polish_name_s6 <- function(x) {
  x <- gsub("^ATH-aux-cistrans-|^ATH-aux-cis-|^ATH-aux-trans-|^ATH-aux-reg-", "", x)
  x <- gsub("^AT-aux-cistrans-|^AT-aux-cis-|^AT-aux-trans-|^AT-aux-reg-",     "", x)
  gsub("PolarAuxinTransport", "PAT", x)
}

organ_rename_s6 <- c(
  cauline = "Cauline 42 DAG", D0_silique = "Silique 0 DPA",
  D2_silique = "Silique 2 DPA", D3_silique = "Silique 3 DPA",
  D4_silique = "Silique 4 DPA", D5_silique = "Silique 5 DPA",
  D6_root = "Root 6 DAG", D11_root = "Root 11 DAG",
  Early_flower = "Early flower", Middle_flower = "Middle flower",
  Late_flower = "Late flower",
  S1_leaf = "Rosette 14 DAG", S2_leaf = "Rosette 21 DAG",
  S3_leaf = "Rosette 28 DAG", S4_leaf = "Rosette 35 DAG",
  S5_leaf = "Rosette 42 DAG", S6_leaf = "Rosette 49 DAG",
  stem = "Stem 42 DAG"
)

organ_order_s6 <- c(
  "Root 6 DAG", "Root 11 DAG",
  "Rosette 14 DAG", "Rosette 21 DAG", "Rosette 28 DAG",
  "Rosette 35 DAG", "Rosette 42 DAG", "Rosette 49 DAG",
  "Early flower", "Middle flower", "Late flower",
  "Silique 0 DPA", "Silique 2 DPA", "Silique 3 DPA",
  "Silique 4 DPA", "Silique 5 DPA",
  "Cauline 42 DAG", "Stem 42 DAG"
)

standardize_tissue_s6 <- function(x) {
  x <- stringr::str_trim(stringr::str_to_lower(x))
  x <- stringr::str_remove(x, "-\\d+$")
  x <- dplyr::if_else(stringr::str_detect(x, "^(phloem|phloem[- ]parenchyma.*|protophloem)$"), "phloem", x)
  x <- dplyr::if_else(stringr::str_detect(x, "^(xylem|protoxylem.*|xylem/procambium)$"),        "xylem", x)
  x <- dplyr::if_else(stringr::str_detect(x, "atrichoblast|trichoblast|epidermis"),              "epidermis", x)
  x <- dplyr::if_else(stringr::str_detect(x, "seed coat.*"),                                     "seed coat", x)
  x <- dplyr::if_else(stringr::str_detect(x, "^(pollen|sperm|tapetum|microsporocyte)$"),         "male gametophyte", x)
  x <- dplyr::if_else(stringr::str_detect(x, "embryo|endosperm"),                               "endosperm/embryo", x)
  x <- dplyr::if_else(x %in% c("vascular", "vascular cells", "interfascicular fibers"),         "vascular (unspecified)", x)
  x <- dplyr::if_else(x %in% c("g2m stage", "s stage"),                                         "dividing cell", x)
  x <- dplyr::if_else(x %in% c("companion cell"),                                               "phloem", x)
  x
}

all_data_s6 <- lapply(names(guo_avg_exp_s6), function(org_key) {
  mat     <- guo_avg_exp_s6[[org_key]]
  org_lab <- organ_rename_s6[org_key]
  if (is.na(org_lab)) org_lab <- org_key
  rownames(mat) <- polish_name_s6(rownames(mat))
  as.data.frame(mat) %>%
    tibble::rownames_to_column("sensor") %>%
    tidyr::pivot_longer(-sensor, names_to = "tissue", values_to = "value") %>%
    dplyr::mutate(organ = org_lab, tissue_simple = standardize_tissue_s6(tissue))
}) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(organ = factor(organ, levels = organ_order_s6))

all_data_s6_agg <- all_data_s6 %>%
  dplyr::group_by(organ, tissue_simple, sensor) %>%
  dplyr::summarise(value = max(value, na.rm = TRUE), .groups = "drop")

tissue_order_s6 <- all_data_s6_agg %>%
  dplyr::distinct(organ, tissue_simple) %>%
  dplyr::count(tissue_simple) %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::pull(tissue_simple)

all_data_s6_agg <- dplyr::mutate(
  all_data_s6_agg,
  tissue_simple = factor(tissue_simple, levels = tissue_order_s6)
)

sensors_s6 <- c("ARF", "Synthesis", "PAT", "ConjugationDeconjugation")
sensors_s6  <- intersect(sensors_s6, unique(all_data_s6_agg$sensor))

# Each panel gets its own independent scale (different sensors have different ranges)
make_heatmap_s6 <- function(df, sensor_name, show_x = FALSE) {
  lim <- quantile(df$value, 0.99, na.rm = TRUE)
  ggplot(df, aes(x = tissue_simple, y = organ, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_distiller(palette = "YlOrRd", direction = 1,
                         limits = c(0, lim), oob = scales::squish,
                         na.value = "grey90", name = sensor_name) +
    labs(title = sensor_name, x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid   = element_blank(),
      axis.text.x  = if (show_x) element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)
                     else element_blank(),
      axis.text.y  = element_text(size = 8),
      plot.title   = element_text(face = "bold", size = 10),
      legend.position = "right"
    )
}

plots_s6 <- lapply(seq_along(sensors_s6), function(i) {
  s <- sensors_s6[i]
  make_heatmap_s6(dplyr::filter(all_data_s6_agg, sensor == s), s,
                  show_x = (i == length(sensors_s6)))
})

combined_s6 <- wrap_plots(plots_s6, ncol = 1)   # each panel keeps its own legend

ggsave("Manuscript-Figures/out/FigureS6_Guo_atlas_heatmap.pdf",
       plot = combined_s6,
       width = 14, height = 3.5 * length(sensors_s6), dpi = 300, bg = "white")
ggsave("Manuscript-Figures/out/FigureS6_Guo_atlas_heatmap.svg",
       plot = combined_s6,
       width = 14, height = 3.5 * length(sensors_s6), dpi = 300, bg = "white")
message("Saved Supplementary Figure 6 heatmap")
