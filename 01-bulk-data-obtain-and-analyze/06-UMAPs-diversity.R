library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(stringr)


install.packages("vctrs")
seurat_obj<- readRDS("D:/!GitHub/iSensors-supplementary/01-bulk-data-obtain-and-analyze/Data/seurat_auxin_microarrays_no_cb.rds")

DimPlot(seurat_obj)
meta_pseudo <- seurat_obj@meta.data


#Polish treatment
seurat_obj$Treatment <- seurat_obj$Type
seurat_obj$Type <- NULL
seurat_obj$Treatment <- dplyr::recode(
  as.character(seurat_obj$Treatment),
  "treat"   = "Auxin",
  "control" = "Control"
)


# Polish Duration 
#Convert "05h" -> 0,5, "1h" -> 1, "72h" -> 72, then re-factor in numeric order
levels(meta_pseudo$Duration)
dur_chr <- as.character(seurat_obj$Duration)

dur_hours <- dplyr::case_when(
  dur_chr == "05h" ~ 0.5,
  TRUE ~ as.numeric(sub("h$", "", dur_chr))
)

# Create ordered factor for plotting (keep original labels)
dur_levels <- unique(dur_chr)[order(dur_hours[match(unique(dur_chr), dur_chr)])]
seurat_obj$Duration <- factor(dur_chr, levels = dur_levels)


seurat_obj$Duration <- factor(seurat_obj$Duration, levels = dur_levels)

Seurat::DiscretePalette(15, palette = "polychrome")

# Polish Concentration  

seurat_obj$sample_col <- rownames(seurat_obj@meta.data)
sample_col <- "sample_col"   # <-- change if your sample id is stored elsewhere

# ---- 1) Build the mapping table from your text ----
map_tbl <- tibble::tribble(
  ~dataset_dir, ~label,
  "auxin_seedling_2h_okushima",            "Auxin (IAA 5uM)",
  "auxin_root_2h_vanneste",                "Auxin (IAA 10uM)",
  "auxin_root_6h_vanneste",                "Auxin (IAA 10uM)",
  "auxin_root_4h_stepanova",               "Auxin (IAA 1uM)",
  "auxin_seedling_05h_delker",             "Auxin (IAA 1uM)",
  "auxin_seedling_1h_delker",              "Auxin (IAA 1uM)",
  "auxin_seedling_3h_delker",              "Auxin (IAA 1uM)",
  "auxin_root_2h_bargmann",                "Auxin (IAA 5uM)",
  "auxin_root_05h_lewis",                  "Auxin (IAA 1uM)",
  "auxin_root_1h_lewis",                   "Auxin (IAA 1uM)",
  "auxin_root_2h_lewis",                   "Auxin (IAA 1uM)",
  "auxin_root_4h_lewis",                   "Auxin (IAA 1uM)",
  "auxin_root_8h_lewis",                   "Auxin (IAA 1uM)",
  "auxin_root_12h_lewis",                  "Auxin (IAA 1uM)",
  "auxin_root_24h_lewis",                  "Auxin (IAA 1uM)",
  "auxin_root_2h_derybel",                 "Auxin (NAA 10uM)",
  "auxin_root_6h_derybel",                 "Auxin (NAA 10uM)",
  "auxin_seedling_3h_shimada",             "Auxin (IAA 1uM)",
  "auxin_root-apex_6h_xuan",               "Auxin (10uM IBA)",
  "auxin_cauline-buds_18h_muller",         "Auxin (NAA 1uM)",
  "auxin_leaves_72h_romero-puertas",       "Auxin (NAA 10uM)" #it is incorrect, jst for simplisity
)

# ---- 2) Parse compound + numeric concentration (convert mM -> uM) ----
map_tbl <- map_tbl %>%
  mutate(
    compound = case_when(
      str_detect(label, "IAA")  ~ "IAA",
      str_detect(label, "NAA")  ~ "NAA",
      str_detect(label, "IBA")  ~ "IBA",
      str_detect(label, "2\\.4D") ~ "2.4D",
      TRUE ~ NA_character_
    ),
    # Extract numeric + unit (uM or mM)
    conc_value = as.numeric(str_extract(label, "\\d+(?:\\.\\d+)?")),
    conc_unit  = str_extract(label, "(uM|mM)"),
    concentration_uM = case_when(
      conc_unit == "uM" ~ conc_value,
      conc_unit == "mM" ~ conc_value * 1000,
      TRUE ~ NA_real_
    )
  )

# ---- 3) Create Concentration in metadata ----
md <- seurat_obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(sample_id = .data[[sample_col]])

# sample_id is like "control_1auxin_root_2h_lewis"? or "treat_1auxin_root_2h_lewis"?
# We'll extract:
#   group = control/treat
#   dataset_dir = everything after the first "_" + replicate prefix
md <- md %>%
  mutate(
    group = case_when(
      str_detect(sample_id, "^control_") ~ "control",
      str_detect(sample_id, "^treat_")   ~ "treat",
      TRUE ~ NA_character_
    ),
    # Remove leading "control_<rep>_" or "treat_<rep>_"
    dataset_dir = str_replace(sample_id, "^(control|treat)_\\d+_?", "")
  )

# Join to mapping table
md2 <- md %>%
  left_join(map_tbl %>% select(dataset_dir, concentration_uM, compound, label),
            by = "dataset_dir") %>%
  mutate(
    Concentration = case_when(
      group == "control" ~ 0,
      group == "treat"   ~ concentration_uM,
      TRUE               ~ NA_real_
    )
  )

# Put back into Seurat object
seurat_obj$Concentration <- md2$Concentration[match(rownames(seurat_obj@meta.data), md2$cell)]

# Optional extra columns (often helpful)
seurat_obj$Auxin_compound <- md2$compound[match(rownames(seurat_obj@meta.data), md2$cell)]
seurat_obj$Auxin_label    <- md2$label[match(rownames(seurat_obj@meta.data), md2$cell)]
seurat_obj$Dataset_dir    <- md2$dataset_dir[match(rownames(seurat_obj@meta.data), md2$cell)]

# ---- 4) Sanity checks ----
table(is.na(seurat_obj$Concentration))
sort(unique(seurat_obj$Concentration))

# show which treated samples failed to map (most important check)
unmapped_treat <- md2 %>%
  filter(group == "treat", is.na(Concentration)) %>%
  distinct(sample_id, dataset_dir)

unmapped_treat




meta_cols <- c(
  "Type",
  "Tissue",
  "Duration",
  "Concentration"
)

# Build one color map PER variable (so it stays identical each time you plot that variable)
cols_treatment  <- c("#7F7F7F", "lightgrey")
cols_tissue     <- c("#98DF8A", "#8C564B", "#BCBD22",  "#FF7F0E")
cols_duration <- Seurat::DiscretePalette(10, palette = "polychrome")
cols_concentration<- c("lightgrey", "#AEC7E8", "#17BECF",  "#1F77B4")


p1 <- DimPlot(seurat_obj, reduction="umap", group.by="Treatment", pt.size=0.5) +
  scale_color_manual(values = cols_treatment, drop = FALSE)+
    theme_classic(base_size = 7) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      legend.title = element_blank(),
      legend.text = element_text(size = 4),
      legend.key.height = unit(0.02, "cm"),
      legend.key.width  = unit(0.02, "cm"),
      legend.position = c(0.98, 0.85),      # bottom-right inside panel
      legend.justification = c(1, 0),
#      legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm"))

p1
out_pdf <- "01-bulk-data-obtain-and-analyze/Out/UMAP_treatment.pdf"

ggsave(out_pdf, plot = p1, width = 4, height = 4, units = "cm", useDingbats = FALSE, dpi = 300)

p2 <- DimPlot(seurat_obj, reduction="umap", group.by="Tissue", pt.size=0.5) +
  scale_color_manual(values = cols_tissue, drop = FALSE)+ 
  theme_classic(base_size = 7) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 4),
    legend.key.height = unit(0.02, "cm"),
    legend.key.width  = unit(0.02, "cm"),
    legend.position = c(1, 0.75),      # bottom-right inside panel
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.65), color = NA),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
  )
p2
out_pdf <- "01-bulk-data-obtain-and-analyze/Out/UMAP_tissue.pdf"

ggsave(out_pdf, plot = p2, width = 4, height = 4, units = "cm", useDingbats = FALSE, dpi = 300)


p3 <- DimPlot(seurat_obj, reduction="umap", group.by="Duration", pt.size=0.5) +
  scale_color_manual(values = cols_duration, drop = FALSE)+ 
  theme_classic(base_size = 7) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 4),
    legend.key.height = unit(0.02, "cm"),
    legend.key.width  = unit(0.02, "cm"),
    legend.position = c(1, 0.20),      # bottom-right inside panel
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.45), color = NA),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
  )

p3
out_pdf <- "01-bulk-data-obtain-and-analyze/Out/UMAP_duration.pdf"

ggsave(out_pdf, plot = p3, width = 4, height = 4, units = "cm", useDingbats = FALSE, dpi = 300)

p4 <- DimPlot(seurat_obj, reduction="umap", group.by="Concentration", pt.size=0.5) +
  scale_color_manual(values = cols_concentration, drop = FALSE)+ 
  theme_classic(base_size = 7) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(size = 8, face = "bold", hjust = 0),
    legend.title = element_blank(),
    legend.text = element_text(size = 4),
    legend.key.height = unit(0.02, "cm"),
    legend.key.width  = unit(0.02, "cm"),
    legend.position = c(1, 0.75),      # bottom-right inside panel
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = alpha("white", 0.65), color = NA),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, unit = "mm")
  )
p4
out_pdf <- "01-bulk-data-obtain-and-analyze/Out/UMAP_concentraion.pdf"

ggsave(out_pdf, plot = p4, width = 4, height = 4, units = "cm", useDingbats = FALSE, dpi = 300)


p5 <- wrap_plots(c(p1,p2,p3,p4), ncol = 1)

p5

out_pdf <- "01-bulk-data-obtain-and-analyze/Out/UMAP_4panels_vertical.pdf"
out_png <- "01-bulk-data-obtain-and-analyze/Out/UMAP_4panels_vertical.png"

# PDF (vector; best for papers)
ggsave(out_pdf, plot = p5, width = 4, height = 10, units = "cm", useDingbats = FALSE)

# PNG (raster; 300–600 dpi for slides/figures)
ggsave(out_png, plot = p5, width = 4, height = 10, units = "cm", dpi = 600, bg = "white")




#alternative way to combine plots

library(ggplot2)
library(patchwork)
library(cowplot)

# helper: extract legend grob
get_leg <- function(p) cowplot::get_legend(
  p + theme(
    legend.position = "right",
    legend.justification = "center",
    legend.text = element_text(size = 5),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width  = unit(0.25, "cm"),
    legend.background = element_rect(fill = "white", color = NA)
  )
)

# helper: bind plot + its legend in two columns
plot_with_leg <- function(p, leg_width = 0.33) {
  leg <- get_leg(p)
  p_noleg <- p + theme(legend.position = "none")
  cowplot::plot_grid(
    p_noleg, leg,
    ncol = 2,
    rel_widths = c(1, leg_width),
    align = "h",
    axis = "tb"
  )
}

# Build each row (adjust leg_width per variable; clusters may need wider)
row1 <- plot_with_leg(p1, leg_width = 0.40)  # e.g. 20 clusters
row2 <- plot_with_leg(p2, leg_width = 0.30)
row3 <- plot_with_leg(p3, leg_width = 0.30)
row4 <- plot_with_leg(p4, leg_width = 0.30)

p5 <- cowplot::plot_grid(row1, row2, row3, row4, ncol = 1, rel_heights = c(1,1,1,1))
p5

ggsave("01-bulk-data-obtain-and-analyze/Out/UMAP_4panels_vertical.pdf", p5,
       width = 4.5, height = 14, units = "cm",
       device = cairo_pdf)
