
library(Seurat)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(ggrepel)

getwd()
# output directory
output_dir <- "D:/!GitHub/iSensors-supplementary/Manuscript-Figures/out"


iSensors_obj <-readRDS(file = "D:/!GitHub/iSensors-supplementary/00-iSensors-objects/data/shahan-iSensors-obj.rds")

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










# comparing significance levels figure 2G ----

library(ggrepel)
# import data
exo_1 <- read.delim("D:/!GitHub/DigitalSensor-Toolbox/in/PerformanceTest1_mean_normed_corrData_2.txt", 
                    header = TRUE)

exo_3 <- read.delim("D:/!GitHub/DigitalSensor-Toolbox/in/PerformanceTest3_mean_normed_corrData_2.txt", 
                    header = TRUE)

endo <- read.csv("D:/!GitHub/DigitalSensor-Toolbox/in/v0.4_rho_values.csv",
                 header = TRUE)
endo$panel <- gsub("cist-", "cis-", endo$panel)

endo <- endo %>%
  mutate(
    cleanpanel = str_replace(panel, "AT-aux-[^-]+-", "")
  )

# make rownames into panel column
exo_1$panel <- rownames(exo_1)
exo_3$panel <- rownames(exo_3)

# Standardize panel names in exo
colnames(exo_1)[colnames(exo_1) == "Panel"] <- "panel"
exo_1$panel <- gsub("_", "-", exo_1$panel)

colnames(exo_3)[colnames(exo_3) == "Panel"] <- "panel"
exo_3$panel <- gsub("_", "-", exo_3$panel)

# Rename columns in exo
colnames(exo_1)[colnames(exo_1) == "all_corr"] <- "exo_rho"
colnames(exo_1)[colnames(exo_1) == "all_pval"] <- "exo_p"

colnames(exo_3)[colnames(exo_3) == "all_corr"] <- "exo_rho"
colnames(exo_3)[colnames(exo_3) == "all_pval"] <- "exo_p"

# Rename columns in endo
colnames(endo)[colnames(endo) == "rho.rho"] <- "endo_rho"
colnames(endo)[colnames(endo) == "p_value"] <- "endo_p"

# Merge the data frames
exo_1_endo <- merge(exo_1, endo, by = "panel", all.x = TRUE)
exo_1_endo_rel <- exo_1_endo[, c("panel", "exo_rho", "exo_p", "endo_rho", "endo_p")]

exo_3_endo <- merge(exo_3, endo, by = "panel", all.x = TRUE)
exo_3_endo_rel <- exo_3_endo[, c("panel", "exo_rho", "exo_p", "endo_rho", "endo_p")]

# as numeric
exo_1_endo_rel$exo_rho <- as.numeric(exo_1_endo_rel$exo_rho)
exo_1_endo_rel$exo_p <- as.numeric(exo_1_endo_rel$exo_p)
exo_1_endo_rel$endo_rho <- as.numeric(exo_1_endo_rel$endo_rho)
exo_1_endo_rel$endo_p <- as.numeric(exo_1_endo_rel$endo_p)

exo_3_endo_rel$exo_rho <- as.numeric(exo_3_endo_rel$exo_rho)
exo_3_endo_rel$exo_p <- as.numeric(exo_3_endo_rel$exo_p)
exo_3_endo_rel$endo_rho <- as.numeric(exo_3_endo_rel$endo_rho)
exo_3_endo_rel$endo_p <- as.numeric(exo_3_endo_rel$endo_p)

# Take log of p-values
exo_1_endo_rel$exo_logp <- -log10(exo_1_endo_rel$exo_p)
exo_1_endo_rel$endo_logp <- -log10(exo_1_endo_rel$endo_p)

exo_3_endo_rel$exo_logp <- -log10(exo_3_endo_rel$exo_p)
exo_3_endo_rel$endo_logp <- -log10(exo_3_endo_rel$endo_p)

# Assign color groups based on significance
exo_1_endo_rel$color_rho <- case_when(
  exo_1_endo_rel$exo_rho > 0.5 & exo_1_endo_rel$endo_rho > 0.5 ~ "Both significant",
  exo_1_endo_rel$exo_rho > 0.5 ~ "Significant, exogenous",
  exo_1_endo_rel$endo_rho > 0.5 ~ "Significant, endogenous",
  TRUE ~ "Not significant"
)

exo_1_endo_rel$color_p_log <- case_when(
  exo_1_endo_rel$exo_logp > -log10(0.05) & exo_1_endo_rel$endo_logp > -log10(0.05) ~ "Both significant",
  exo_1_endo_rel$exo_logp > -log10(0.05) ~ "Significant, exogenous",
  exo_1_endo_rel$endo_logp > -log10(0.05) ~ "Significant, endogenous",
  TRUE ~ "Not significant"
)

exo_3_endo_rel$color_rho <- case_when(
  exo_3_endo_rel$exo_rho > 0.5 & exo_3_endo_rel$endo_rho > 0.5 ~ "Both significant",
  exo_3_endo_rel$exo_rho > 0.5 ~ "Significant, exogenous",
  exo_3_endo_rel$endo_rho > 0.5 ~ "Significant, endogenous",
  TRUE ~ "Not significant"
)

exo_3_endo_rel$color_p_log <- case_when(
  exo_3_endo_rel$exo_logp > -log10(0.05) & exo_3_endo_rel$endo_logp > -log10(0.05) ~ "Both significant",
  exo_3_endo_rel$exo_logp > -log10(0.05) ~ "Significant, exogenous",
  exo_3_endo_rel$endo_logp > -log10(0.05) ~ "Significant, endogenous",
  TRUE ~ "Not significant"
)

# cleanpanel join
exo_1_endo_rel <- exo_1_endo_rel %>%
  left_join(endo %>% select(panel, cleanpanel), by = "panel")

exo_3_endo_rel <- exo_3_endo_rel %>%
  left_join(endo %>% select(panel, cleanpanel), by = "panel")

# Rename cleanpanel values
exo_1_endo_rel$cleanpanel <- recode(exo_1_endo_rel$cleanpanel,
                                    "R2D2" = "Response",
                                    "R2D2v2" = "Response2",
                                    "R2D2v3" = "Response3",
                                    "R2D2v4" = "Response4",
                                    "R2D2v5" = "Response5",
                                    "R2D2v6" = "Response6",
                                    "EffluxInflux" = "PAT")
exo_3_endo_rel$cleanpanel <- recode(exo_3_endo_rel$cleanpanel,
                                    "R2D2" = "Response",
                                    "R2D2v2" = "Response2",
                                    "R2D2v3" = "Response3",
                                    "R2D2v4" = "Response4",
                                    "R2D2v5" = "Response5",
                                    "R2D2v6" = "Response6",
                                    "EffluxInflux" = "PAT")

# Plot scale
x_limits <- c(0, 13)
y_limits <- c(0, 7)
x_breaks <- seq(0, 13, by = 1)
y_breaks <- seq(0, 7, by = 1)



p1 <- ggplot(
  exo_1_endo_rel,
  aes(x = exo_logp, y = endo_logp, color = color_rho, shape = color_p_log)
) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Both significant" = "#1f78b4",
    "Significant, exogenous" = "#33a02c",
    "Significant, endogenous" = "#e31a1c",
    "Not significant" = "gray70"
  )) +
  scale_shape_manual(values = c(
    "Both significant" = 15,  # square
    "Significant, exogenous" = 19,
    "Significant, endogenous" = 17,
    "Not significant" = 18
  )) +
  # LABELS driven by p-significance
  geom_text_repel(
    data = subset(exo_1_endo_rel, color_p_log == "Both significant"),
    aes(x = exo_logp, y = endo_logp, label = cleanpanel),  # <- supply x & y
    inherit.aes = FALSE,
    color = "gray40",
    size = 4,
    max.overlaps = 100
  ) +
  scale_x_continuous(breaks = x_breaks, limits = x_limits) +
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +
  # Clean background, keep axes, bigger fonts
  theme_classic(base_size = 16) +
  theme(
    axis.line   = element_line(linewidth = 0.6),  # <- linewidth, not size
    axis.ticks  = element_line(linewidth = 0.4),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  ) +
  labs(
    x = expression(-log[10]("Exogenous concentration, p-value")),
    y = expression(-log[10]("Endogenous concentration, p-value")),
    color = "rho",
    shape = "p-value"
  )

p1

file_path <- file.path(output_dir, "exo1-endo-comparison.svg")

ggsave(file_path, plot = p1, width = 7, height = 5, dpi = 300, bg = "white")


# Plot for exo_3_endo_rel
p3 <- ggplot(
  exo_3_endo_rel,
  aes(x = exo_logp, y = endo_logp, color = color_rho, shape = color_p_log)
) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Both significant" = "#1f78b4",
    "Significant, exogenous" = "#33a02c",
    "Significant, endogenous" = "#e31a1c",
    "Not significant" = "gray70"
  )) +
  scale_shape_manual(values = c(
    "Both significant" = 15,  # square
    "Significant, exogenous" = 19,
    "Significant, endogenous" = 17,
    "Not significant" = 18
  )) +
  # LABELS driven by p-significance
  geom_text_repel(
    data = subset(exo_3_endo_rel, color_p_log == "Both significant"),
    aes(x = exo_logp, y = endo_logp, label = cleanpanel),  # <- supply x & y
    inherit.aes = FALSE,
    color = "gray40",
    size = 4,
    max.overlaps = 100
  ) +
  scale_x_continuous(breaks = x_breaks, limits = x_limits) +
  scale_y_continuous(breaks = y_breaks, limits = y_limits) +
  # Clean background, keep axes, bigger fonts
  theme_classic(base_size = 16) +
  theme(
    axis.line   = element_line(linewidth = 0.6),  # <- linewidth, not size
    axis.ticks  = element_line(linewidth = 0.4),
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  ) +
  labs(
    x = expression(-log[10]("Exogenous duration, p-value")),
    y = expression(-log[10]("Endogenous concentration, p-value")),
    color = "rho",
    shape = "p-value"
  )

p3
file_path <- file.path(output_dir, "exo3-endo-comparison.svg")

ggsave(file_path, plot = p3, width = 7, height = 5, dpi = 300, bg = "white")

##############################
# Define R2D2 versions (after recoding)
R2D2_versions <- c("Response", "Response2", "Response3", "Response4", "Response5", "Response6")

# Filter R2D2 only and non-R2D2 panels for exo_1_endo_rel
exo_1_R2D2 <- exo_1_endo_rel %>% filter(cleanpanel %in% R2D2_versions)
exo_1_nonR2D2 <- exo_1_endo_rel %>% filter(!cleanpanel %in% R2D2_versions)

# Filter R2D2 only and non-R2D2 panels for exo_3_endo_rel
exo_3_R2D2 <- exo_3_endo_rel %>% filter(cleanpanel %in% R2D2_versions)
exo_3_nonR2D2 <- exo_3_endo_rel %>% filter(!cleanpanel %in% R2D2_versions)

# Function for R2D2 plots
plot_logp_R2D2 <- function(data) {
  ggplot(data, aes(x = exo_logp, y = endo_logp, color = color_rho, shape = color_p_log)) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      "Both significant" = "#1f78b4",
      "Significant for exogenous data" = "#33a02c",
      "Significant for endogenous data" = "#e31a1c",
      "Other" = "gray70"
    )) +
    scale_shape_manual(values = c(
      "Both significant" = 15,
      "Significant for exogenous data" = 19,
      "Significant for endogenous data" = 17,
      "Other" = 18
    )) +
    geom_text_repel(
      aes(label = cleanpanel),
      size = 3,
      max.overlaps = 10
    ) +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    theme_minimal() +
    labs(
      x = expression(-log[10]("Exogenous auxin p-value")),
      y = expression(-log[10]("Endogenous auxin p-value")),
      color = "rho significance",
      shape = "p significance"
    )
}

# Function for non-R2D2 plots
plot_logp_nonR2D2 <- function(data) {
  ggplot(data, aes(x = exo_logp, y = endo_logp, color = color_rho, shape = color_p_log)) +
    geom_point(size = 3) +
    scale_color_manual(values = c(
      "Both significant" = "#1f78b4",
      "Significant for exogenous data" = "#33a02c",
      "Significant for endogenous data" = "#e31a1c",
      "Other" = "gray70"
    )) +
    scale_shape_manual(values = c(
      "Both significant" = 15,
      "Significant for exogenous data" = 19,
      "Significant for endogenous data" = 17,
      "Other" = 18
    )) +
    geom_text_repel(
      data = subset(data, color_rho %in% c("Both significant", "Significant for endogenous data")),
      aes(label = cleanpanel),
      size = 3,
      max.overlaps = 10
    ) +
    scale_x_continuous(breaks = x_breaks, limits = x_limits) +
    scale_y_continuous(breaks = y_breaks, limits = y_limits) +
    theme_minimal() +
    labs(
      x = expression(-log[10]("Exogenous auxin p-value")),
      y = expression(-log[10]("Endogenous auxin p-value")),
      color = "rho significance",
      shape = "p significance"
    )
}

# Create and save the four plots
ggsave("Shahan/REDO_counts/exo_1_R2D2.png", plot = plot_logp_R2D2(exo_1_R2D2), width = 7, height = 5, dpi = 300, bg = "white")
ggsave("Shahan/REDO_counts/exo_1_nonR2D2.png", plot = plot_logp_nonR2D2(exo_1_nonR2D2), width = 7, height = 5, dpi = 300, bg = "white")
ggsave("Shahan/REDO_counts/exo_3_R2D2.png", plot = plot_logp_R2D2(exo_3_R2D2), width = 7, height = 5, dpi = 300, bg = "white")
ggsave("Shahan/REDO_counts/exo_3_nonR2D2.png", plot = plot_logp_nonR2D2(exo_3_nonR2D2), width = 7, height = 5, dpi = 300, bg = "white")


#Figure 3A epidermis ----

assay_name <- "iSensor_mean_normed"
feature <- "R2D2v2"

ae <- AverageExpression(iSensors_obj, assays = assay_name, features = feature, slot = "data", verbose = FALSE)[[assay_name]]

ae <- as.matrix(ae)                      # coerce if it’s a vector/data.frame

# Use the first (and only) row
val <- as.numeric(ae[1, ])
clu <- colnames(ae)
stopifnot(length(val) == length(clu))

df <- tibble(cluster = clu, value = val) %>%
  arrange(desc(value)) %>%                                 # largest on top
  mutate(cluster = factor(cluster, levels = rev(cluster))) # order y-axis

# Horizontal bar “heatmap”: clusters on Y, value on X (log scale)
library(scales)
p <- ggplot(df, aes(x = value, y = cluster)) +
  geom_col(width = 0.7, fill = "#fdbf6f") +        # fixed bar color
  scale_x_continuous(
    trans  = scales::log1p_trans(),                # log(1+x)
    breaks = scales::pretty_breaks(n = 3)
  ) +
  labs(x = "Average expression (log1p scale)", y = NULL) +  # no y-axis title
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.4),
    legend.position = "none"                       # hide any legend
  )

p

p <- ggplot(df, aes(x = value, y = cluster)) +
  geom_col(width = 0.7, fill = "#fdbf6f") +
  scale_x_continuous(
    trans   = scales::log1p_trans(),
    breaks  = scales::breaks_pretty(n = 2),            # fewer ticks
    guide   = guide_axis(n.dodge = 2),                 # dodge if still tight
    expand  = expansion(mult = c(0, 0.02))
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid   = element_blank(),
    axis.line    = element_line(linewidth = 0.6),
    axis.ticks   = element_line(linewidth = 0.4),
    legend.position = "none"
  )

p
file_path <- file.path(output_dir, "Response2_barplot_epidermis.svg")

ggsave(file_path, plot = p, width = 4, height = 5, dpi = 300, bg = "white")



#Figure 3B -----
in_path <- "E:/Users/Maria/Genetic sensor paper/Genetic sensor paper/"
seurat_obj <- readRDS(paste0(in_path, "ra_integrated_Shahan_20250114.rds"))

meta.data <- seurat_obj@meta.data


Idents(seurat_obj)<- factor(seurat_obj[["CellTypes", drop = TRUE]])
DimPlot(seurat_obj)
#exploring original annotation ----

clusters_to_keep <- c("Epidermis")

cells <- WhichCells(seurat_obj, idents = clusters_to_keep)
obj_sub <- subset(seurat_obj, cells = cells)

DimPlot(obj_sub)

obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub<- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:20)
obj_sub <- FindClusters(obj_sub, resolution = 1)  # adjust as needed
obj_sub <- RunUMAP(obj_sub, dims = 1:50)
DimPlot(obj_sub, reduction = "umap", label = TRUE, label.size = 5)

exclude_clusters <- c("16", "22", "25")
obj_sub <- subset(obj_sub, idents = exclude_clusters, invert = TRUE)

obj_sub <- NormalizeData(obj_sub)
obj_sub <- FindVariableFeatures(obj_sub)
obj_sub <- ScaleData(obj_sub)
obj_sub<- RunPCA(obj_sub)
obj_sub <- FindNeighbors(obj_sub, dims = 1:20)
obj_sub <- FindClusters(obj_sub, resolution = 0.5)  # adjust as needed
obj_sub <- RunUMAP(obj_sub, dims = 1:50)
DimPlot(obj_sub, reduction = "umap", label = TRUE, label.size = 5)

obj_sub <- RenameIdents(
  obj_sub,
  `0`  = "T",  `15` = "T",  `4`  = "T",  `11` = "T",  `10` = "T",  `2` = "T",  `3` = "T",
  `1`  = "AT", `14` = "AT", `13` = "AT", `6`  = "AT", `7`  = "AT", `9` = "AT", `5` = "AT",
  `8`  = "T_AT", `12` = "T_AT"
)

# (optional) set the order of the new identities
Idents(obj_sub) <- factor(Idents(obj_sub), levels = c("T", "AT", "T_AT"))

# Quick check
levels(Idents(obj_sub))
cols <- c("#fb9a99", "#a6cee3", "#fdbf6f")
p <- DimPlot(
  obj_sub,
  reduction = "umap",
  label = FALSE, label.size = 3,
  cols = cols
)+NoLegend()

p
file_path <- file.path(output_dir, "DimPlot_epidermis.svg")

ggsave(file_path, plot = p, width = 3, height = 3, dpi = 300, bg = "white")

in_dir <- "D:/!GitHub/DigitalSensor-Toolbox/in"
file_path <- file.path(in_dir, "shahan_epidermis.rds")
file_path
saveRDS(obj_sub, file_path)

obj_sub <- readRDS(file_path)

obj_sub <- iSensor_pipeline(seuratObject=obj_sub, useParallel = FALSE, signals = c("mean_normed"),
                            metaPanels = list('Response' = list('srcPanels' = c('AT_aux_trans_ARF', 'AT_aux_trans_IAA'),
                                                                'rule' = prod)))
Assays(obj_sub)


DefaultAssay(obj_sub) <- "iSensor_mean_normed"

p_transport <- FeaturePlot(
  obj_sub,
  features  = "AT-aux-trans-EffluxInflux",
  reduction = "umap",
  order     = TRUE,
  pt.size   = 0.6,
  min.cutoff = "q05",
  max.cutoff = "q95",
  cols      = c("grey90", "navy")
) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.4),
    plot.title = element_blank()
  ) 

p_transport

file_path <- file.path(output_dir, "shahan_epidermis_PAT_featureplot.svg")
file_path
ggsave(file_path, plot = p_transport,
       width = 4, height = 3, dpi = 300, bg = "white")


p_synthesis <- FeaturePlot(
  obj_sub,
  features  = "AT-aux-trans-Synthesis",
  reduction = "umap",
  order     = TRUE,
  pt.size   = 0.6,
  min.cutoff = "q05",
  max.cutoff = "q95",
  cols      = c("grey90", "navy")
) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.4),
    plot.title = element_blank()
  ) 

p_synthesis

file_path <- file.path(output_dir, "shahan_epidermis_synthesis_featureplot.svg")
file_path
ggsave(file_path, plot = p_synthesis,
       width = 4, height = 3, dpi = 300, bg = "white")


p_response <- FeaturePlot(
  obj_sub,
  features  = "Response",
  reduction = "umap",
  order     = TRUE,
  pt.size   = 0.6,
  min.cutoff = "q05",
  max.cutoff = "q95",
  cols      = c("grey90", "navy")
) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.4),
    plot.title = element_blank()
  ) 

p_response

file_path <- file.path(output_dir, "shahan_epidermis_response_featureplot.svg")
file_path
ggsave(file_path, plot = p_response,
       width = 4, height = 3, dpi = 300, bg = "white")
