
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(ggrepel)
library(ggPlantmap)

getwd()
# output directory
output_dir <- "D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/Manuscript-Figures/out"


iSensors_obj <-readRDS(file = "D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/shahan-iSensors-obj-groundtruth.rds")

DimPlot(iSensors_obj)
iSensors_obj@active.ident
paired12 <- c(
  "#e31a1c", "#1f78b4", "#b2df8a", "#33a02c",
  "#fb9a99", "#a6cee3", "#fdbf6f", "#ff7f00",
  "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
)

if (length(levels(Idents(iSensors_obj))) > length(paired12)) {
  times <- ceiling(length(levels(Idents(iSensors_obj))) / length(paired12))
  paired12 <- rep(paired12, times = times)
}

new_order <- c("Initials", "QC", "CSC", "CSCD", "Columella", "Young LRC", "LRC","Dying LRC", 
               "Protoxylem_m1", "Protoxylem_m2", 
               "Procambium_m1", "Procambium_m2", 
               "Pericycle_m1","Pericycle_m2", 
               "Cortex_m1", "Cortex_m2", 
               "Endodermis_m1", "Endodermis_m2",
               "Atrichoblast_m1", "Atrichoblast_m2", 
               "Trichoblast_m1", "Trichoblast_m2", "LRP")

colnames(iSensors_obj@meta.data)
iSensors_obj$iSensors <- as.character(Idents(iSensors_obj))

iSensors_obj$iSensors_obj <- factor(
  iSensors_obj$iSensors,
  levels = new_order
)

p_old <- DimPlot(
  iSensors_obj,
  reduction = "umap",
  group.by = "iSensors_obj",
  cols = paired12[seq_along(levels(Idents(iSensors_obj)))],
  label = FALSE) + theme(plot.title = element_blank()) +  
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = 15,    # square
        size  = 3      # bigger size
      ),
      ncol = 1         # single column legend
    ))

print(p_old)

file_path <- file.path(output_dir, "Figure4B_dimplot.svg")
ggsave(file_path, plot = p_old,
       width = 6, height = 6, dpi = 300)


p_old_nolegend <- DimPlot(
  iSensors_obj,
  reduction = "umap",
  group.by = "iSensors_obj",
  cols = paired12[seq_along(levels(Idents(iSensors_obj)))],
  label = FALSE) + theme(plot.title = element_blank()) +  NoLegend()

print(p_old_nolegend)

file_path <- file.path(output_dir, "Figure4B_dimplot_noLegend.svg")
file_path
ggsave(file_path, plot = p_old_nolegend,
       width = 4, height = 6, dpi = 300)

file_path <- file.path(output_dir, "Figure4B_dimplot_noLegend.pdf")
file_path
ggsave(file_path, plot = p_old_nolegend,
       width = 4, height = 6, dpi = 300)

#Calculation iSensors for small shahan dataset ----

DefaultAssay(iSensors_obj)<- "iSensors_mean"

feat = "ATH-aux-trans-ARF"

vals <- FetchData(iSensors_obj, vars = feat)[, 1]
lim  <- max(abs(quantile(vals, probs = c(0.02, 0.98), na.rm = TRUE)))
lims <- c(0, lim)

rdylbu5 <- rev(brewer.pal(n = 5, name = "RdYlBu"))

p <- FeaturePlot(
  object = iSensors_obj,
  features = feat,
  combine = TRUE
)
p <- p &
  scale_color_gradientn(
    colors = rdylbu5,
    limits = lims,
    oob = scales::squish
  ) &
  guides(
    color = guide_colorbar(
      #      title = "iSensor activity",
      #      title.position = "top",
      title.hjust = 0.4,
      barwidth = unit(3, "cm"),
      barheight = unit(0.4, "cm"),
      ticks = FALSE
    )
  ) &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    axis.line  = element_blank(),
    axis.ticks = element_blank(),
    axis.text  = element_blank(),
    axis.title = element_blank()
  )


p

file_path <- file.path(output_dir, "Figure4B_featureplot_ARF.pdf")
file_path
ggsave(file_path, plot = p,
       width = 4, height = 6, dpi = 300)

file_path <- file.path(output_dir, "Figure4B_featureplot_ARF.svg")
file_path
ggsave(file_path, plot = p,
       width = 4, height = 6, dpi = 300)


# figure 4A -----
# Conceptual diagram — not generated in R.

# ── Shared setup ──────────────────────────────────────────────────────────────
library(stringr)
library(dplyr)

# Print all available sensor names so you can verify name spellings below
message("Available iSensors in iSensors_mean assay:")
print(rownames(iSensors_obj[["iSensors_mean"]]))




# Ground truth auxin rank per cell type (Petersson et al. 2009; Brunoud et al. 2012)
# QC = rank 1 (highest endogenous auxin), outer epidermis = rank 12 (lowest)
# "Dying LRC" is excluded — no clear rank in the literature.
auxin_rank <- c(
  "QC"              = 1,
  "CSC"             = 2,
  "LRP"             = 2,
  "Columella"       = 3,
  "CSCD"            = 3,
  "Young LRC"       = 3,
  "Initials"        = 4,
  "Protoxylem-m1"   = 5,
  "LRC"             = 6,
  "Protoxylem-m2"   = 6,
  "Procambium-m1"   = 7,
  "Procambium-m2"   = 7,
  "Endodermis-m1"   = 8,
  "Pericycle-m1"    = 8,
  "Endodermis-m2"   = 9,
  "Pericycle-m2"    = 9,
  "Cortex-m1"       = 10,
  "Cortex-m2"       = 10,
  "Atrichoblast-m1" = 11,
  "Trichoblast-m1"  = 11,
  "Atrichoblast-m2" = 12,
  "Trichoblast-m2"  = 12
)

# Average iSensor signal per cluster (slot = "data" matches how the object was created)

DefaultAssay(iSensors_obj) <- "iSensors_mean"
avg_exp  <- AverageExpression(iSensors_obj, assays = "iSensors_mean", slot = "data",
                              verbose = FALSE)
avg_mat  <- avg_exp[["iSensors_mean"]]

# Filter to clusters present in both the object and the ground truth rank
shared_clusters <- intersect(names(auxin_rank), colnames(avg_mat))
mat_filt        <- avg_mat[, shared_clusters, drop = FALSE]
gt_rank         <- auxin_rank[shared_clusters]


# ── Figure 4C: 2×2 Spearman scatter plots ─────────────────────────────────────
library(patchwork)

sensors_4c <- c(
  "ATH-aux-trans-ARF",
  "ATH-aux-trans-PolarAuxinTransport",
  "ATH-aux-trans-Synthesis",
  "ATH-aux-trans-IAA"
)
labels_4c <- c("ARF", "PAT", "Synthesis", "IAA")

make_scatter_4c <- function(feat, label) {
  if (!feat %in% rownames(mat_filt)) {
    message("Sensor not found: ", feat, " — check rownames printed above")
    return(NULL)
  }
  sensor_vals <- as.numeric(mat_filt[feat, ])
  sensor_rank <- rank(-sensor_vals, ties.method = "average")
  rho_test    <- cor.test(sensor_rank, gt_rank, method = "spearman", exact = FALSE)
  rho_val     <- round(as.numeric(rho_test$estimate), 2)
  p_val       <- rho_test$p.value
  p_lab       <- if (p_val < 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))

  plot_df <- data.frame(ground_truth = gt_rank, predicted = sensor_rank)

  ggplot(plot_df, aes(x = ground_truth, y = predicted)) +
    geom_point(color = "#b15928", size = 2.5) +
    geom_smooth(method = "lm", color = "grey50", se = FALSE, linewidth = 0.6) +
    scale_x_continuous(breaks = seq(2, 12, by = 2)) +
    labs(
      title    = label,
      subtitle = paste0("Rho = ", rho_val, "   ", p_lab),
      x        = "Ground truth rank",
      y        = "iSensor rank"
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 10),
      plot.subtitle = element_text(hjust = 0.5, size = 8, color = "grey40"),
      axis.line     = element_line(linewidth = 0.4),
      axis.ticks    = element_line(linewidth = 0.3),
      panel.grid    = element_blank()
    )
}

plots_4c <- mapply(make_scatter_4c, sensors_4c, labels_4c, SIMPLIFY = FALSE)
plots_4c <- Filter(Negate(is.null), plots_4c)

p_4c <- wrap_plots(plots_4c, ncol = 2)

p_4c
ggsave(file.path(output_dir, "Figure4C_spearman_2x2.svg"),
       plot = p_4c, width = 6, height = 4, dpi = 300)
ggsave(file.path(output_dir, "Figure4C_spearman_2x2.pdf"),
       plot = p_4c, width = 6, height = 4, dpi = 300)
message("Saved Figure 4C (2×2)")


# ── Figure 4D: Rho bar plot for all iSensors ─────────────────────────────────
all_sensors <- rownames(mat_filt)

spearman_list <- lapply(all_sensors, function(feat) {
  vals <- as.numeric(mat_filt[feat, ])
  rk   <- rank(-vals, ties.method = "average")
  test <- cor.test(rk, gt_rank, method = "spearman", exact = FALSE)
  data.frame(panel = feat, rho = as.numeric(test$estimate),
             p_value = test$p.value, stringsAsFactors = FALSE)
})
spearman_df <- bind_rows(spearman_list)

spearman_df <- spearman_df %>%
  mutate(
    reg_class = case_when(
      str_detect(panel, "reg") ~ "reg",
      str_detect(panel, "cis")      ~ "cis",
      str_detect(panel, "trans")    ~ "trans",
      TRUE                          ~ "other"
    ),
    clean_panel = str_replace(panel, "^ATH-aux-[^-]+-", ""),
    clean_panel = str_replace(clean_panel, "PolarAuxinTransport", "PAT"),
    significance = if_else(rho > 0.5 & p_value < 0.05,
                           "significant", "not significant")
  ) %>%
  arrange(rho) %>%
  mutate(clean_panel = factor(clean_panel, levels = unique(clean_panel)))

write.csv(spearman_df,
          file.path(output_dir, "Figure4D_spearman_sc-endo.csv"),
          row.names = FALSE)

barPlot_4d <- ggplot(spearman_df,
                     aes(x = rho, y = clean_panel, fill = reg_class,
                         alpha = significance)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(cis   = "#fdbf6f", trans = "#ff7f00",
                                reg   = "#b15928", other = "gray")) +
  scale_alpha_manual(values = c("significant" = 1, "not significant" = 0.3)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
  geom_vline(xintercept = 0.5, linewidth = 0.4, linetype = "dashed",
             color = "grey40") +
  scale_x_continuous(breaks = seq(-0.5, 1, by = 0.1)) +
  labs(x = "Spearman correlation (Rho)", y = NULL, fill = "iSensor type") +
  guides(fill = guide_legend(order = 1), alpha = guide_legend(order = 2)) +
  theme_minimal() +
  theme(panel.grid  = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.line.x = element_line())

barPlot_4d
ggsave(file.path(output_dir, "Figure4D_spearman_barplot.svg"),
       plot = barPlot_4d, width = 8, height = 10, dpi = 300)
ggsave(file.path(output_dir, "Figure4D_spearman_barplot.pdf"),
       plot = barPlot_4d, width = 8, height = 10, dpi = 300)


# ── Figure 4E & 4H: Realistic root layouts with ggRootCellAtlas ──────────────

if (!requireNamespace("ggPlantmap", quietly = TRUE)) {
  devtools::install_github("leonardojo/ggPlantmap")
} else {
  message("ggPlantmap is already installed")
}

#devtools::install_github("MariaSavina91/ggRootCellAtlas", force = TRUE)

if (!requireNamespace("ggRootCellAtlas", quietly = TRUE)) {
  devtools::install_github("MariaSavina91/ggRootCellAtlas")
} else {
  message("ggRootCellAtlas is already installed")
}
library(ggRootCellAtlas)

# This does not work yet!!! Correct data is missing...

iSensors_obj <- readRDS("D:/!GitHub/DigitalSensor-Toolbox/iSensors-supplementary/00-iSensors-objects/data/shahan-roottip-iSensors-obj.rds")
Assays(iSensors_obj)
DefaultAssay(iSensors_obj)<- "iSensors_mean"
avg_exp  <- AverageExpression(iSensors_obj, assays = "iSensors_mean", slot = "data",
                              verbose = FALSE)
avg_mat_rrca  <- avg_exp[["iSensors_mean"]]

rownames(avg_mat_rrca) <- gsub("^ATH-aux-[^-]+-", "", rownames(avg_mat_rrca))
rownames(avg_mat_rrca) <- gsub("PolarAuxinTransport", "PAT", rownames(avg_mat_rrca))
avg_exp_rrca <- list(RNA = avg_mat_rrca)

message("Sensor names available for root atlas plots:")
print(rownames(avg_exp_rrca[["RNA"]]))

# Figure 4E: five sensors shown in the manuscript
# Adjust names below if they differ from the printed list above.
sensors_4e <- c("Receptors", "PAT", "Synthesis", "IAA", "majortrend")

for (sensor in sensors_4e) {
  if (!sensor %in% rownames(avg_exp_rrca[["RNA"]])) {
    message("Sensor not found for 4E: ", sensor,
            " — available: ", paste(rownames(avg_exp_rrca[["RNA"]]), collapse = ", "))
    next
  }
  vals <- as.numeric(avg_exp_rrca[["RNA"]][sensor, ])
  c2   <- quantile(vals[is.finite(vals)], 0.98, na.rm = TRUE)

  p <- ggRootCellAtlas_expression(avg_exp_rrca, sensor)
  ggsave(file.path(output_dir, paste0("Figure4E_", sensor, "_root.pdf")),
         plot = p, width = 6, height = 4, dpi = 300)
  ggsave(file.path(output_dir, paste0("Figure4E_", sensor, "_root.svg")),
         plot = p, width = 6, height = 4, dpi = 300)
  message("Saved Figure 4E: ", sensor)
}

# Figure 4H: ARF iSensor on root layout
if ("ARF" %in% rownames(avg_exp_rrca[["RNA"]])) {
  vals_arf <- as.numeric(avg_exp_rrca[["RNA"]]["ARF", ])
  c2_arf   <- quantile(vals_arf[is.finite(vals_arf)], 0.98, na.rm = TRUE)

  p_4h <- ggRootCellAtlas_expression(avg_exp_rrca, "ARF", c1 = 0, c2 = c2_arf)
  ggsave(file.path(output_dir, "Figure4H_ARF_root.pdf"),
         plot = p_4h, width = 6, height = 4, dpi = 300)
  ggsave(file.path(output_dir, "Figure4H_ARF_root.svg"),
         plot = p_4h, width = 6, height = 4, dpi = 300)
  message("Saved Figure 4H")
}




# ── Figure 4F: Summary performance across benchmarks (arrow plot) ─────────────
# This panel requires statistics CSV files in iSensors-supplementary/Manuscript-Figures/in/:
#   - 2026-06-03_PerformanceTest1_stat.csv  (bulk Spearman — concentration)
#   - 2026-06-03_PerformanceTest3_stat.csv  (bulk Spearman — duration)
#   - Statistics-exo-scdata.csv             (exogenous sc-data; generated by Figure2.R)
# The sc-endo Spearman is taken directly from spearman_df computed above.

in_dir_stats <- "iSensors-supplementary/Manuscript-Figures/in"

polish_name <- function(x) {
  x <- gsub("^ATH-aux-cistrans-|^ATH-aux-cis-|^ATH-aux-trans-|^ATH-aux-reg-", "", x)
  x <- gsub("^AT-aux-cistrans-|^AT-aux-cis-|^AT-aux-trans-|^AT-aux-reg-",     "", x)
  gsub("PolarAuxinTransport", "PAT", x)
}

read_stat_file <- function(path) {
  df <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) {
    message("Could not read: ", path)
    return(NULL)
  }
  # Ensure first or second column is the iSensor name
  if (!any(tolower(names(df)) == "isensor")) df <- df[, -1, drop = FALSE]
  names(df)[1] <- "iSensor"
  df$iSensor <- polish_name(df$iSensor)
  df
}

bulk_conc <- read_stat_file(file.path(in_dir_stats, "2026-06-03_PerformanceTest1_stat.csv"))
bulk_dur  <- read_stat_file(file.path(in_dir_stats, "2026-06-03_PerformanceTest3_stat.csv"))
sc_exo    <- read_stat_file(file.path(in_dir_stats, "Statistics-exo-scdata.csv"))
sc_endo_raw <- read.csv(file.path(in_dir_stats, "Figure4D_spearman_sc-endo.csv"),
                       stringsAsFactors = FALSE)
sc_endo <- sc_endo_raw %>%
  transmute(iSensor  = gsub("PolarAuxinTransport", "PAT", clean_panel),
            rho_endo = rho,
            p_endo   = p_value)
sc_endo$iSensor <- as.character(sc_endo$iSensor)

# Keep only sensors with rho_endo > 0.5 (the significant ones shown in 4F)
sig_sensors <- sc_endo$iSensor[sc_endo$rho_endo > 0.5 & sc_endo$p_endo < 0.05]

if (length(sig_sensors) > 0 && !is.null(bulk_conc) && !is.null(bulk_dur) && !is.null(sc_exo)) {

  # Identify relevant columns — adjust column names if they differ in your files
  # bulk_conc: Rho column + p-value column
  # bulk_dur:  Rho column + p-value column
  # sc_exo:    effect/estimate column + adj p-value column
  # Merge all into one long table keyed on iSensor
  merged_stats <- sc_endo %>%
    filter(iSensor %in% sig_sensors)

  # Build long-format arrow table
  metric_list <- list()

  # Helper: prefer BH-adjusted p; fall back to raw p column
  pick_p_col <- function(nms) {
    adj <- grep("^BH[0-9]|p_adj|padj|\\.adj", nms, ignore.case = TRUE, value = TRUE)
    if (length(adj)) return(adj[1])
    grep("^p[0-9]|p\\.val|p_val|pval", nms, ignore.case = TRUE, value = TRUE)[1]
  }

  if (!is.null(bulk_conc)) {
    rho_col <- grep("^Rho|^rho|estimate", names(bulk_conc), value = TRUE)[1]
    p_col   <- pick_p_col(names(bulk_conc))
    if (!is.na(rho_col)) {
      tmp <- bulk_conc[bulk_conc$iSensor %in% sig_sensors, c("iSensor", rho_col, p_col)]
      names(tmp) <- c("iSensor", "value", "pval")
      tmp$metric <- "bulk_conc"
      metric_list[["bulk_conc"]] <- tmp
    }
  }
  if (!is.null(bulk_dur)) {
    rho_col <- grep("^Rho|^rho|estimate", names(bulk_dur), value = TRUE)[1]
    p_col   <- pick_p_col(names(bulk_dur))
    if (!is.na(rho_col)) {
      tmp <- bulk_dur[bulk_dur$iSensor %in% sig_sensors, c("iSensor", rho_col, p_col)]
      names(tmp) <- c("iSensor", "value", "pval")
      tmp$metric <- "bulk_dur"
      metric_list[["bulk_dur"]] <- tmp
    }
  }
  if (!is.null(sc_exo)) {
    eff_col <- grep("effect|estimate", names(sc_exo), ignore.case = TRUE, value = TRUE)[1]
    p_col   <- grep("p_adj|padj|p.adj", names(sc_exo), ignore.case = TRUE, value = TRUE)[1]
    if (!is.na(eff_col)) {
      tmp <- sc_exo[sc_exo$iSensor %in% sig_sensors, c("iSensor", eff_col, p_col)]
      names(tmp) <- c("iSensor", "value", "pval")
      tmp$metric <- "sc_exo"
      metric_list[["sc_exo"]] <- tmp
    }
  }
  sc_endo_long <- sc_endo %>%
    filter(iSensor %in% sig_sensors) %>%
    transmute(iSensor, value = rho_endo, pval = p_endo, metric = "sc_endo")
  metric_list[["sc_endo"]] <- sc_endo_long

  df_long_4f <- bind_rows(metric_list) %>%
    mutate(
      value = suppressWarnings(as.numeric(value)),
      pval  = suppressWarnings(as.numeric(pval)),
      sig   = !is.na(pval) & pval < 0.05,
      arrow = case_when(is.na(value) ~ NA_character_,
                        value > 0    ~ "↑",
                        value < 0    ~ "↓",
                        TRUE         ~ "→"),
      class = case_when(!sig      ~ "n.s.",
                        value > 0 ~ "up",
                        value < 0 ~ "down",
                        TRUE      ~ "zero"),
      metric = factor(metric,
                      levels = c("bulk_conc", "bulk_dur", "sc_exo", "sc_endo"),
                      labels = c("Bulk conc.", "Bulk dur.", "sc exo", "sc endo"))
    )

  mag_lim <- quantile(abs(df_long_4f$value[is.finite(df_long_4f$value)]), 0.98,
                      na.rm = TRUE)

  p_4f <- ggplot(df_long_4f, aes(x = metric, y = iSensor)) +
    geom_text(aes(label = arrow, size = pmin(abs(value), mag_lim), color = class),
              na.rm = TRUE) +
    scale_color_manual(values = c(up = "#D73027", down = "#4575B4",
                                   "n.s." = "grey75", zero = "grey75")) +
    scale_size(range = c(2.5, 7.5), guide = "none") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1),
          axis.text.y = element_text(size = 8),
          axis.ticks  = element_blank(),
          legend.position = "none")

  ggsave(file.path(output_dir, "Figure4F_performance_summary.svg"),
         plot = p_4f, width = 3, height = 3.5, dpi = 300)
  ggsave(file.path(output_dir, "Figure4F_performance_summary.pdf"),
         plot = p_4f, width = 3, height = 3.5, dpi = 300)
  message("Saved Figure 4F")
} else {
  message("Figure 4F skipped: run Figure2.R first to generate Statistics-exo-scdata.csv, then check all files exist in ", in_dir_stats)
}
