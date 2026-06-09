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
