library(Seurat)
library(iSensors)
library(tidyverse)
library(RColorBrewer)
library(scales)
library(pheatmap)


iSensors_obj <- readRDS(file = "02-single-cell-data-analysis/in/iSensors-Martin-Arevalillo-auxin-root2025.rds")

#izualizing fold change heatmap

condition_col <- "orig.ident2"
celltype_col  <- "Cell.Ident"   # change to "Cell.Ident" if that is what you use

levels(iSensors_obj@meta.data$Cell.Ident)


tissues <- levels(Idents(iSensors_obj))

de_list <- list()



for (tissue in tissues) {
  cells_t <- WhichCells(iSensors_obj, idents = tissue)
  
  if (length(cells_t) < 20) next
  
  obj_t <- subset(iSensors_obj, cells = cells_t)
  
  conds <- unique(as.character(obj_t$orig.ident2))
  if (!all(c("Control", "Auxin") %in% conds)) next
  
  de <- FindMarkers(
    object = obj_t,
    ident.1 = "Auxin",
    ident.2 = "Control",
    group.by = "orig.ident2",
    features = isensors_list,
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0
  )
  
  if (nrow(de) == 0) next
  
  de <- de %>%
    rownames_to_column("iSensor") %>%
    mutate(
      tissue = tissue,
      p_val_adj = p.adjust(p_val, method = "BH"),
      neglog10_padj = -log10(p_val_adj + 1e-300)
    )
  
  de_list[[tissue]] <- de
}

de_all <- bind_rows(de_list)

# Save table
write.csv(de_all, "02-single-cell-data-analysis/in/iSensors_DE_aux_vs_ctr_by_tissue.csv", row.names = FALSE)



# Vizualize auxin levels and responses on heatmaps ----
# ---------------------------
avg_list <- AverageExpression(
  iSensors_obj,
  assays   = "iSensors_mean_normed",
  features = isensors_list,
  group.by = c("orig.ident2", "Cell.Ident"),
  slot     = "data",   # typically use normalized/log-normalized values
  verbose  = FALSE
)


avg_mat <- avg_list[["iSensors_mean_normed"]]  # genes x (orig.ident2_tissue)

# Split into ctr and aux matrices with tissues as columns
# Column names created by Seurat look like "ctr_tissueA", "aux_tissueA", etc.
ctr_cols <- grep("^Control_", colnames(avg_mat), value = TRUE)
aux_cols <- grep("^Auxin_", colnames(avg_mat), value = TRUE)

ctr_mat <- avg_mat[, ctr_cols, drop = FALSE]
aux_mat <- avg_mat[, aux_cols, drop = FALSE]

# Clean column names to just tissue names
colnames(ctr_mat) <- sub("^Control_", "", colnames(ctr_mat))
colnames(aux_mat) <- sub("^Auxin_", "", colnames(aux_mat))

# Align tissues (keep common tissues, same order)
common_tissues <- intersect(colnames(ctr_mat), colnames(aux_mat))
ctr_mat <- ctr_mat[, common_tissues, drop = FALSE]
aux_mat <- aux_mat[, common_tissues, drop = FALSE]

# Optional: order tissues in a meaningful way (currently uses common_tissues order)
# common_tissues <- sort(common_tissues)
# ctr_mat <- ctr_mat[, common_tissues, drop = FALSE]
# aux_mat <- aux_mat[, common_tissues, drop = FALSE]

# ---------------------------
# 2) Two heatmaps (ctr and aux) with a shared scale
# ---------------------------
# Choose a common color scale range based on BOTH matrices


# Compute shared limits from BOTH matrices using robust quantiles
vals <- c(as.numeric(ctr_mat), as.numeric(aux_mat))
lims <- quantile(vals, probs = c(0.01, 0.99), na.rm = TRUE)  # try 0.02/0.98 if needed

n_breaks <- 101
breaks_shared <- seq(lims[1], lims[2], length.out = n_breaks)

pheatmap(
  ctr_mat,
  breaks = breaks_shared,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  main = "Average expression per tissue, control",
  border_color = NA,
  filename = "02-single-cell-data-analysis/out/Heatmap_CTR_shared_scale.pdf",
  width = 10,
  height = 12
)



pheatmap(
  aux_mat,
  breaks = breaks_shared,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  main = "Average expression per tissue, auxin",
  border_color = NA,
  filename = "02-single-cell-data-analysis/out/Heatmap_AUX_shared_scale.pdf",
  width = 10,
  height = 12
)



# ---------------------------
# 3) Ratio heatmap (AUX / CTR), recommended as log2 fold-change
# ---------------------------
# Add a small pseudocount to avoid division by zero
pseudocount <- 1e-6

log2_ratio <- log2((aux_mat + pseudocount) / (ctr_mat + pseudocount))

vals_r <- as.numeric(log2_ratio)
lims_r <- quantile(vals_r, c(0.01, 0.99), na.rm = TRUE)
breaks_ratio <- seq(lims_r[1], lims_r[2], length.out = 101)

pheatmap(
  log2_ratio,
  breaks = breaks_ratio,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  main = "log2(Auxin/Control)",
  border_color = NA,
  filename = "02-single-cell-data-analysis/out/Heatmap_log2_AUX_vs_CTR.pdf",
  width = 10,
  height = 12
)

# heatmaps with significance level

# Make sure row/col names match your heatmap matrices
log2_ratio_m <- as.matrix(log2_ratio)

# 1) Canonicalize names: convert '-' to '_' and remove any accidental whitespace
canon <- function(x) gsub("-", "_", trimws(x))

# Apply to heatmap matrix columns
colnames(log2_ratio_m) <- canon(colnames(log2_ratio_m))

# Apply to DE table tissue names
de_all$tissue <- canon(de_all$tissue)

colnames(log2_ratio_m) 
unique(de_all$tissue)

stopifnot(all(colnames(log2_ratio_m) %in% unique(de_all$tissue)))
stopifnot(all(rownames(log2_ratio_m) %in% unique(de_all$iSensor)))

# Pivot adjusted p-values into a matrix: rows=iSensors, cols=tissues
padj_df <- de_all %>%
  select(iSensor, tissue, p_val_adj) %>%
  distinct() %>%
  tidyr::pivot_wider(names_from = tissue, values_from = p_val_adj)

padj_mat <- as.data.frame(padj_df)
rownames(padj_mat) <- padj_mat$iSensor
padj_mat$iSensor <- NULL
padj_mat <- as.matrix(padj_mat)

# Reorder to match your heatmap matrices exactly
padj_mat <- padj_mat[rownames(log2_ratio_m), colnames(log2_ratio_m), drop = FALSE]

# Define significance (you can change the threshold)
alpha <- 0.05
sig_mat <- padj_mat < alpha
sig_mat[is.na(sig_mat)] <- FALSE


#building pheatmaps
log2_ratio_masked <- log2_ratio_m
log2_ratio_masked[!sig_mat] <- NA  # hide non-significant cells

# Robust scale for visualization (prevents one outlier from dominating)
vals_m <- as.numeric(log2_ratio_masked)
vals_m <- vals_m[is.finite(vals_m)]
lims_m <- quantile(vals_m, probs = c(0.01, 0.99), na.rm = TRUE)

breaks_masked <- seq(lims_m[1], lims_m[2], length.out = 101)

pheatmap(
  log2_ratio_masked,
  breaks = breaks_masked,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  main = "log2(AUX/CTR) per tissue (non-significant masked)",
  border_color = NA,
  na_col = "grey90",
  filename = "02-single-cell-data-analysis/out/Heatmap_log2_AUX_vs_CTR_masked.pdf",
  width = 10,
  height = 12
)


# Create star matrix aligned to log2_ratio
stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat),
                    dimnames = dimnames(padj_mat))

stars_mat[padj_mat < 0.01]  <- "*"
stars_mat[padj_mat < 0.001]  <- "**"
stars_mat[padj_mat < 0.00001] <- "***"

# Robust scale for visualization
vals_r <- as.numeric(log2_ratio_m)
lims_r <- quantile(vals_r, probs = c(0.01, 0.99), na.rm = TRUE)
breaks_ratio <- seq(lims_r[1], lims_r[2], length.out = 101)

pheatmap(
  log2_ratio_m,
  breaks = breaks_ratio,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  main = "log2(AUX/CTR) per tissue (significance annotated)",
  border_color = NA,
  display_numbers = stars_mat,
  number_color = "black",
  fontsize_number = 8,
  filename = "out/Auxin/Heatmap_log2_AUX_vs_CTR_annotated.pdf",
  width = 10,
  height = 12
)

lfc_cut <- 0.25
stars_mat2 <- stars_mat
stars_mat2[abs(log2_ratio_m) < lfc_cut] <- ""

pheatmap(
  log2_ratio_m,
  breaks = breaks_ratio,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 6,
  main = "log2(AUX/CTR) (stars only if |log2FC| >= 0.25)",
  border_color = NA,
  display_numbers = stars_mat2,
  number_color = "black",
  fontsize_number = 8,
  filename = "out/Auxin/Heatmap_log2_AUX_vs_CTR_annotated_lfcFiltered.pdf",
  width = 10,
  height = 12
)