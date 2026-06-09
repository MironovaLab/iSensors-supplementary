library(Seurat)
library(iSensors)
library(ggPlantmap)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(scales)

output_dir  <- "iSensors-supplementary/Manuscript-Figures/out"
rdylbu5     <- rev(brewer.pal(n = 5, name = "RdYlBu"))

# ── Load full Shahan root object (used for Figure 5A) ─────────────────────────
# Object already contains iSensors_mean assay — no recalculation needed.
shahan_obj <- readRDS("iSensors-supplementary/00-iSensors-objects/data/shahan-iSensors-obj.rds")
DefaultAssay(shahan_obj) <- "iSensors_mean"

# Print sensor names once to verify spellings used below
message("Available sensors in iSensors_mean:")
print(rownames(shahan_obj[["iSensors_mean"]]))


# ── Figure 5A (right panel): Bar plot of ARF iSensor per cell type ────────────
# Shows average ARF activity across root-tip cell clusters (log1p scale).
# Left panel of 5A is the realistic root layout — generated in Figure4.R (4H/4E).

arf_feat <- "ATH-aux-trans-ARF"  # adjust if printed name differs

ae <- AverageExpression(shahan_obj, assays = "iSensors_mean",
                        features = arf_feat, slot = "data",
                        verbose = FALSE)[["iSensors_mean"]]
ae <- as.matrix(ae)
val <- as.numeric(ae[1, ])
clu <- colnames(ae)

clusters_bar <- c(
  "QC", "CSC", "CSCD", "Columella", "Initials",
  "Lateral-Root-Cap",       # from "LRC"
  "Lateral-Root-Primordia", # from "LRP"
  "Protoxylem_m1", "Protoxylem_m2",
  "Procambium_m1", "Procambium_m2",
  "Pericycle_m1",  "Pericycle_m2",
  "Cortex_m1",     "Cortex_m2",
  "Endodermis_m1", "Endodermis_m2",
  "Atrichoblast_m1", "Atrichoblast_m2",
  "Trichoblast_m1",  "Trichoblast_m2"
)

cluster_rename_bar <- c(
  "LRP"  = "Lateral-Root-Primordia",
  "LRC"  = "Lateral-Root-Cap",
  "Young LRC"  = "Young-LRC",
  "Dying LRC"  = "Dying-LRC"
)

bar_df <- tibble(cluster = clu, value = val) %>%
  mutate(cluster = if_else(cluster %in% names(cluster_rename_bar),
                           cluster_rename_bar[cluster], cluster)) %>%
  filter(cluster %in% clusters_bar) %>%
  arrange(desc(value)) %>%
  mutate(cluster = factor(cluster, levels = rev(cluster)))

p_5a_bar <- ggplot(bar_df, aes(x = value, y = cluster)) +
  geom_col(width = 0.7, fill = "#fdbf6f") +
  scale_x_continuous(trans  = scales::log1p_trans(),
                     breaks = scales::pretty_breaks(n = 3)) +
  labs(x = "Average expression (log1p scale)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(panel.grid    = element_blank(),
        axis.line     = element_line(linewidth = 0.6),
        axis.ticks    = element_line(linewidth = 0.4),
        legend.position = "none")

ggsave(file.path(output_dir, "Figure5A_ARF_barplot_by_celltype.svg"),
       plot = p_5a_bar, width = 4, height = 5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure5A_ARF_barplot_by_celltype.pdf"),
       plot = p_5a_bar, width = 4, height = 5, dpi = 300, bg = "white")
message("Saved Figure 5A bar plot")


# ── Figure 5A (left): Root tip layout with ARF iSensor (ggRootCellAtlas) ──────
library(ggRootCellAtlas)

avg_all_5a <- AverageExpression(shahan_obj, assays = "iSensors_mean",
                                slot = "data", verbose = FALSE)[["iSensors_mean"]]
rownames(avg_all_5a) <- gsub("^ATH-aux-[^-]+-", "", rownames(avg_all_5a))
rownames(avg_all_5a) <- gsub("PolarAuxinTransport", "PAT", rownames(avg_all_5a))
avg_exp_5a <- list(RNA = avg_all_5a)

if ("ARF" %in% rownames(avg_exp_5a[["RNA"]])) {
  vals_5a <- as.numeric(avg_exp_5a[["RNA"]]["ARF", ])
  c2_5a   <- quantile(vals_5a[is.finite(vals_5a)], 0.98, na.rm = TRUE)

  p_5a_layout <- ggRootCellAtlas_expression(avg_exp_5a, "ARF", c1 = 0, c2 = c2_5a)
  ggsave(file.path(output_dir, "Figure5A_ARF_root_layout.pdf"),
         plot = p_5a_layout, width = 6, height = 4, dpi = 300)
  ggsave(file.path(output_dir, "Figure5A_ARF_root_layout.svg"),
         plot = p_5a_layout, width = 6, height = 4, dpi = 300)
  message("Saved Figure 5A root layout")
}


# ── Figure 5B: Epidermal UMAP + FeaturePlots ─────────────────────────────────
# Load the epidermal subset; iSensors_mean must be calculated here.
epidermis_obj <- readRDS(
  "iSensors-supplementary/00-iSensors-objects/data/Shahan_epidermis_020625.rds"
)

# Check if iSensors_mean already exists; calculate if not.
if (!"iSensors_mean" %in% Assays(epidermis_obj)) {
  message("Calculating iSensors_mean for epidermis object …")
  AuxinPanel <- LoadSensors(setName    = "Auxin",
                            species    = "ATH",
                            hormone    = "aux",
                            customPanels = FALSE,
                            randomInfo = list(n = 3, sizes = c(100, 200, 300),
                                              majortrend = TRUE))
  epidermis_obj <- CalcSensors(
    epidermis_obj,
    seurLayer = "data",
    panelSet  = AuxinPanel,
    signals   = "mean"
  )
  message("iSensors_mean calculated. Assays: ", paste(Assays(epidermis_obj), collapse = ", "))
}

DefaultAssay(epidermis_obj) <- "iSensors_mean"

# UMAP DimPlot (Figure 5B left)
cols_epi <- c("#fb9a99", "#a6cee3", "#fdbf6f")  # T, AT, mixed (if present)
p_5b_dim <- DimPlot(epidermis_obj, reduction = "umap", label = FALSE,
                    cols = cols_epi) +
  theme(axis.line  = element_blank(), axis.ticks = element_blank(),
        axis.text  = element_blank(), axis.title = element_blank(),
        plot.title = element_blank())

ggsave(file.path(output_dir, "Figure5B_epidermis_DimPlot.svg"),
       plot = p_5b_dim, width = 4, height = 3, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure5B_epidermis_DimPlot.pdf"),
       plot = p_5b_dim, width = 4, height = 3, dpi = 300, bg = "white")
message("Saved Figure 5B DimPlot")

# FeaturePlots for ARF, Synthesis, PAT (Figure 5B right panels)
# Verify these sensor names match rownames(epidermis_obj[["iSensors_mean"]])
sensors_5b <- c(
  "ATH-aux-trans-ARF",
  "ATH-aux-trans-Synthesis",
  "ATH-aux-trans-PolarAuxinTransport"  # PAT
)
names_5b   <- c("ARF", "Synthesis", "PAT")

for (i in seq_along(sensors_5b)) {
  feat     <- sensors_5b[i]
  lab      <- names_5b[i]

  if (!feat %in% rownames(epidermis_obj[["iSensors_mean"]])) {
    message("Sensor not found: ", feat, " — skipping")
    next
  }

  vals <- FetchData(epidermis_obj, vars = feat)[, 1]
  lim  <- quantile(vals, probs = 0.98, na.rm = TRUE)

  p_fp <- FeaturePlot(epidermis_obj, features = feat,
                      reduction = "umap", order = TRUE, pt.size = 0.6,
                      min.cutoff = "q05", max.cutoff = lim) +
    scale_color_gradientn(colors = rdylbu5,
                          limits = c(0, lim), oob = scales::squish) +
    guides(color = guide_colorbar(barwidth  = unit(0.4, "cm"),
                                  barheight = unit(3,   "cm"),
                                  ticks = FALSE)) +
    labs(title = lab) +
    theme(legend.position = "right",
          axis.line  = element_blank(), axis.ticks = element_blank(),
          axis.text  = element_blank(), axis.title = element_blank())

  ggsave(file.path(output_dir, paste0("Figure5B_", lab, "_FeaturePlot.svg")),
         plot = p_fp, width = 4, height = 3, dpi = 300, bg = "white")
  ggsave(file.path(output_dir, paste0("Figure5B_", lab, "_FeaturePlot.pdf")),
         plot = p_fp, width = 4, height = 3, dpi = 300, bg = "white")
  message("Saved Figure 5B FeaturePlot: ", lab)
}


# ── Figure 5C: Realistic epidermal layout ─────────────────────────────────────
# The epidermal ggPlantmap layout (AT vs T columns across developmental zones)
# was created from a custom SVG by Carmen — see CW_scripts_all/CW_AT_T_atlas.R.
# That script produces new_ggPlantmap_working.csv from EpidermisT_AT_Carmen.svg.
#
# Option 1: If the CSV has already been generated, load it here and proceed:
#
epi_map_path <- "iSensors-supplementary/Manuscript-Figures/in/new_ggPlantmap_epidermis.csv"

if (file.exists(epi_map_path)) {
  epi_map <- read.csv(epi_map_path, stringsAsFactors = FALSE)
  epi_map$x <- as.numeric(epi_map$x)
  epi_map$y <- as.numeric(epi_map$y)

  # Average iSensor signal per epidermal cluster
  avg_epi    <- AverageExpression(epidermis_obj, assays = "iSensors_mean",
                                  slot = "data", verbose = FALSE)[["iSensors_mean"]]
  avg_epi_df              <- as.data.frame(t(avg_epi))
  avg_epi_df$Cell_type    <- gsub("_", "-", rownames(avg_epi_df))
  colnames(avg_epi_df)    <- gsub("^ATH-aux-[^-]+-", "", colnames(avg_epi_df))
  colnames(avg_epi_df)    <- gsub("PolarAuxinTransport", "PAT", colnames(avg_epi_df))
  colnames(avg_epi_df)    <- gsub("-", "_", colnames(avg_epi_df))

  sensors_5c <- c("ARF", "Synthesis", "PAT")

  for (sensor in sensors_5c) {
    if (!sensor %in% colnames(avg_epi_df)) {
      message("Column not found for 5C: ", sensor)
      next
    }
    merged_epi <- ggPlantmap.merge(epi_map, avg_epi_df, "Cell_type", "Cell_type") %>%
      filter(is.finite(x), is.finite(y))

    p_5c <- ggPlantmap_heatmap(merged_epi, sensor)
    ggsave(file.path(output_dir, paste0("Figure5C_", sensor, "_epi_layout.svg")),
           plot = p_5c, width = 2.7, height = 3.6, dpi = 300, bg = "white")
    ggsave(file.path(output_dir, paste0("Figure5C_", sensor, "_epi_layout.pdf")),
           plot = p_5c, width = 2.7, height = 3.6, dpi = 300, bg = "white")
    message("Saved Figure 5C: ", sensor)
  }

} else {
  message(
    "Figure 5C skipped: epidermal ggPlantmap layout not found at ", epi_map_path, "\n",
    "Run CW_scripts_all/CW_AT_T_atlas.R to generate it from the epidermal SVG,\n",
    "then save the resulting new.ggPlantmap as:\n  ", epi_map_path
  )
}
