#!/usr/bin/env Rscript

# Create:
# 1. An all-gene Wang rice root object containing only cells with confident
#    Arabidopsis Atlas_new label-transfer assignments.
# 2. A smaller object whose RNA assay contains only genes from the Auxin OSA
#    panels, with iSensors mean signals calculated.
# 3. FeaturePlots for all OSA auxin trans panels.
#
# Run from iSensors-supplementary/:
# Rscript 00-iSensors-objects/create-confident-wang-iSensors.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(iSensors)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(scales)
})

set.seed(100)

input_rds <- paste0(
  "00-iSensors-objects/data/arabidopsis_rice_transfer/",
  "wang_rice_root_label_transferred_all_genes_umap.rds"
)
output_all_genes <- paste0(
  "00-iSensors-objects/data/",
  "wang_rice_root_confident_transferred_all_genes.rds"
)
output_isensors <- paste0(
  "00-iSensors-objects/data/",
  "wang_rice_root_confident_Auxin_OSA_iSensors.rds"
)
output_panel_summary <- paste0(
  "00-iSensors-objects/data/",
  "wang_rice_root_confident_Auxin_OSA_panel_genes.tsv"
)
local_conjugation_tsv <- paste0(
  "00-iSensors-objects/data/",
  "OSA-aux-trans-ConjugationDeconjugation_all_high_confidence.tsv"
)
plot_dir <- paste0(
  "00-iSensors-objects/data/",
  "wang_rice_root_confident_Auxin_OSA_FeaturePlots"
)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

log_step("Loading all-gene transferred Wang root object")
rice <- readRDS(input_rds)
if (!"predicted.Atlas_new" %in% colnames(rice[[]])) {
  stop("Input object lacks predicted.Atlas_new")
}
if (!"umap.allgenes" %in% Reductions(rice)) {
  stop("Input object lacks umap.allgenes")
}

confident_cells <- colnames(rice)[
  as.character(rice$predicted.Atlas_new) != "Uncertain"
]
if (!length(confident_cells)) stop("No confidently transferred cells found")

log_step(
  "Keeping ", length(confident_cells), " of ", ncol(rice),
  " cells with confident transferred labels"
)
rice_confident <- subset(rice, cells = confident_cells)
Idents(rice_confident) <- "predicted.Atlas_new"

# Keep all RNA genes, counts/data, metadata, PCA and all-gene UMAP.
rice_confident <- DietSeurat(
  rice_confident,
  assays = "RNA",
  dimreducs = c("pca", "umap.allgenes"),
  graphs = NULL,
  layers = c("counts", "data"),
  misc = TRUE
)
saveRDS(rice_confident, output_all_genes, compress = FALSE)
log_step(
  "Saved confident all-gene object: ", ncol(rice_confident),
  " cells x ", nrow(rice_confident), " genes"
)

auxin_panel <- LoadSensors(
  setName = "Auxin",
  species = "OSA",
  hormone = "aux",
  customPanels = FALSE,
  random = FALSE
)

# The released OSA registry currently lacks ConjugationDeconjugation. Build a
# local panel from high-confidence Ensembl Plants orthologs of the installed
# Arabidopsis panel, allowing one-to-one, one-to-many, and many-to-many matches.
ath_auxin_panel <- LoadSensors(
  setName = "Auxin",
  species = "ATH",
  hormone = "aux",
  customPanels = FALSE,
  random = FALSE
)
ath_conjugation_name <- "ATH-aux-trans-ConjugationDeconjugation"
ath_conjugation <- ath_auxin_panel[["panels"]][[ath_conjugation_name]]
if (is.null(ath_conjugation)) {
  stop("Installed ATH panel set lacks ", ath_conjugation_name)
}

conjugation_map <- read.delim(
  local_conjugation_tsv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_map_columns <- c(
  "arabidopsis_gene", "rice_gene", "orthology_type",
  "orthology_confidence", "rice_identity_to_ath", "ath_identity_to_rice"
)
if (!all(required_map_columns %in% colnames(conjugation_map))) {
  stop(
    "Expanded conjugation mapping lacks columns: ",
    paste(setdiff(required_map_columns, colnames(conjugation_map)), collapse = ", ")
  )
}
conjugation_map <- conjugation_map[
  conjugation_map$arabidopsis_gene %in% ath_conjugation &
    conjugation_map$orthology_confidence == 1 &
    conjugation_map$orthology_type %in% c(
      "ortholog_one2one", "ortholog_one2many", "ortholog_many2many"
    ),
  ,
  drop = FALSE
]
conjugation_map$present_in_Wang_RNA <-
  conjugation_map$rice_gene %in% rownames(rice_confident)
conjugation_map$source_panel <- ath_conjugation_name
conjugation_map$local_panel <- "OSA-aux-trans-ConjugationDeconjugation"

local_conjugation_genes <- unique(
  conjugation_map$rice_gene[conjugation_map$present_in_Wang_RNA]
)
if (!length(local_conjugation_genes)) {
  stop("No one-to-one rice orthologs for the ATH conjugation panel")
}
auxin_panel[["panels"]][["OSA-aux-trans-ConjugationDeconjugation"]] <-
  local_conjugation_genes
log_step(
  "Added expanded OSA ConjugationDeconjugation panel: ",
  length(local_conjugation_genes), " rice genes representing ",
  length(unique(conjugation_map$arabidopsis_gene[
    conjugation_map$present_in_Wang_RNA
  ])), "/", length(ath_conjugation), " ATH genes"
)

trans_panels <- auxin_panel[["panels"]]
trans_panel_names <- names(trans_panels)
panel_genes <- unique(unlist(trans_panels, use.names = FALSE))
genes_found <- intersect(panel_genes, rownames(rice_confident))

panel_summary <- do.call(rbind, lapply(trans_panel_names, function(panel) {
  genes <- unique(trans_panels[[panel]])
  data.frame(
    panel = panel,
    gene = genes,
    present_in_object = genes %in% rownames(rice_confident),
    stringsAsFactors = FALSE
  )
}))
write.table(
  panel_summary,
  output_panel_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (!length(genes_found)) stop("No OSA Auxin panel genes found in RNA assay")
missing_genes <- setdiff(panel_genes, genes_found)
log_step(
  "OSA Auxin genes found: ", length(genes_found), "/", length(panel_genes),
  if (length(missing_genes)) paste0("; missing: ", length(missing_genes)) else ""
)

# Subsetting features preserves cell metadata and dimensional reductions.
rice_panel <- subset(rice_confident, features = genes_found)
DefaultAssay(rice_panel) <- "RNA"

log_step("Calculating iSensors mean signals for seven OSA trans panels")
rice_panel <- CalcSensors(
  rice_panel,
  seurLayer = "data",
  panelSet = auxin_panel,
  signals = "mean"
)

if (!"iSensors_mean" %in% Assays(rice_panel)) {
  stop("CalcSensors did not create the iSensors_mean assay")
}

rice_panel@misc$local_OSA_ConjugationDeconjugation <- list(
  source_panel = ath_conjugation_name,
  mapping = paste(
    "Ensembl Plants high-confidence orthologs:",
    "one-to-one, one-to-many, and many-to-many"
  ),
  arabidopsis_panel_genes = length(ath_conjugation),
  represented_arabidopsis_genes = unique(
    conjugation_map$arabidopsis_gene[conjugation_map$present_in_Wang_RNA]
  ),
  mapped_rice_genes = local_conjugation_genes,
  mapping_table = local_conjugation_tsv,
  caution = paste(
    "Expanded orthology-derived panel. Paralog expansion means rice gene",
    "families may contribute multiple genes per Arabidopsis source gene."
  )
)

saveRDS(rice_panel, output_isensors, compress = FALSE)
log_step(
  "Saved panel-gene iSensors object: ", ncol(rice_panel),
  " cells, ", nrow(rice_panel[["RNA"]]), " RNA genes, ",
  nrow(rice_panel[["iSensors_mean"]]), " iSensors"
)

DefaultAssay(rice_panel) <- "iSensors_mean"
available_sensors <- rownames(rice_panel[["iSensors_mean"]])
trans_features <- intersect(trans_panel_names, available_sensors)
if (!length(trans_features)) {
  stop(
    "No expected trans-panel names found in iSensors_mean. Available: ",
    paste(available_sensors, collapse = ", ")
  )
}

feature_plot <- function(feature) {
  vals <- FetchData(rice_panel, vars = feature)[, 1L]
  upper <- unname(quantile(vals, 0.98, na.rm = TRUE))
  if (!is.finite(upper) || upper <= 0) upper <- max(vals, na.rm = TRUE)

  FeaturePlot(
    rice_panel,
    features = feature,
    reduction = "umap.allgenes",
    order = TRUE,
    pt.size = 0.28,
    min.cutoff = 0,
    max.cutoff = upper,
    raster = TRUE,
    raster.dpi = c(350, 350)
  ) +
    scale_color_gradientn(
      colors = viridis(100, option = "C"),
      limits = c(0, upper),
      oob = squish,
      na.value = "grey92"
    ) +
    labs(
      title = feature,
      color = "iSensor\nmean"
    ) +
    guides(color = guide_colorbar(
      title.position = "top",
      barheight = grid::unit(3.2, "cm")
    )) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      legend.title = element_text(face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )
}

plots <- setNames(lapply(trans_features, feature_plot), trans_features)
combined <- wrap_plots(plots, ncol = 3) +
  plot_annotation(
    title = "Auxin OSA trans-panel iSensors in confidently annotated rice root cells"
  )

ggsave(
  file.path(plot_dir, "Auxin_OSA_trans_panels_FeaturePlots.png"),
  combined,
  width = 15,
  height = 9,
  dpi = 350,
  bg = "white"
)
ggsave(
  file.path(plot_dir, "Auxin_OSA_trans_panels_FeaturePlots.pdf"),
  combined,
  width = 15,
  height = 9,
  device = cairo_pdf
)

for (feature in trans_features) {
  safe_name <- gsub("[^A-Za-z0-9_-]+", "_", feature)
  ggsave(
    file.path(plot_dir, paste0(safe_name, "_FeaturePlot.png")),
    plots[[feature]],
    width = 6,
    height = 5,
    dpi = 350,
    bg = "white"
  )
  ggsave(
    file.path(plot_dir, paste0(safe_name, "_FeaturePlot.pdf")),
    plots[[feature]],
    width = 6,
    height = 5,
    device = cairo_pdf
  )
}

log_step("Saved FeaturePlots to ", plot_dir)
