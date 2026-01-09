#This script is to generate a testing dataset for iSensors turorial. It uses single cell dataset publusehd in (Martin-Arevalillo et al., Cell, 2025)
library(Seurat)
library(iSensors)
library(tidyverse)
library(RColorBrewer)
library(scales)


library(stringr)
library(AnnotationDbi)
library(org.At.tair.db)

rdylbu5 <- rev(brewer.pal(n = 5, name = "RdYlBu"))

#The full dataset is shared via: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE241573

#Loading the dataset, subsetting and reclustering ----

load("C:/!Victoria/!Projects/DigitalSensors/Data/Martin-Arevalillo-2025/GSE241573_Seurat_Object_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13.RData")
ls()

seurat_obj <-
  UpdateSeuratObject(Seurat_Object_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13)

DimPlot(
  seurat_obj,
  reduction = "umap",
  split.by = "orig.ident2",
  combine = TRUE,label = TRUE
)
meta <- seurat_obj@meta.data

head(Idents(seurat_obj))
seu <- subset(seurat_obj, idents = "Initials_Stele_QC")
DimPlot(seu,   split.by = "orig.ident2")

Assays(seu)
Reductions(seu)
seu <- subset(seurat_obj, idents = "Initials_Stele_QC")
DefaultAssay(seu) <- "integrated"


# Optional: clear only graph/clustering results (keep reductions if you want)
seu@graphs <- list()
seu@neighbors <- list()
seu$seurat_clusters <- NULL

seu <- RunPCA(seu, npcs = 30, verbose = FALSE)

dims_use <- 1:30
seu <- FindNeighbors(seu, reduction = "pca", dims = dims_use, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = dims_use, verbose = FALSE)

# Diagnostics: should mix by treatment, separate by tissue/cell identity
DimPlot(seu, group.by = "orig.ident2") + ggtitle("Treatment mixing (orig.ident2)")
DimPlot(seu, split.by = "orig.ident2")

#Clusters annotation ----

FeaturePlot(seu, features = "AT3G11260", label = TRUE, split.by = "orig.ident2")
FeaturePlot(seu, features = "AT5G66770", label = TRUE, split.by = "orig.ident2")

# --- Marker panels (SYMBOLS; easier to read/maintain) ---
markers_sym <- list(
  QC = c("WOX5", "AGL42", "BRAVO", "PLT1", "PLT2"),
  Stele_procambium_initials = c("TMO5", "LHW", "ATHB8"),
  CEI_endodermis_lineage = c("SCR", "SHR", "CYCD6;1"),
  Endodermis = c("MYB36", "CASP1"),
  Columella_stem_cells = c("FEZ", "ACR4"),
  Rootcap_LRC_program = c("SMB", "BRN1", "BRN2"),
  Epidermis_atrichoblast = c("WER", "GL2"),
  Epidermis_trichoblast = c("CPC", "RHD6"),
  Xylem = c("VND7", "IRX3"),
  Phloem = c("APL"),
  Pericycle = c("LBD16", "GATA23")
)

# --- Map SYMBOL -> TAIR (AGI) using org.At.tair.db ---
all_sym <- unique(unlist(markers_sym))

map <- AnnotationDbi::select(
  org.At.tair.db,
  keys = all_sym,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "TAIR")
) %>%
  distinct(SYMBOL, TAIR) %>%
  filter(!is.na(TAIR))

# Helper: map symbol list -> TAIR list (drop missing mappings)
map_syms_to_tair <- function(sym_vec, map_tbl) {
  out <- map_tbl %>%
    filter(SYMBOL %in% sym_vec) %>%
    arrange(match(SYMBOL, sym_vec)) %>%
    pull(TAIR) %>%
    unique()
  out
}


# Use YOUR marker lists only (much safer than pulling every SYMBOL in the DB)
all_sym <- unique(unlist(markers_sym))

label_map <- AnnotationDbi::select(
  org.At.tair.db,
  keys = all_sym,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "TAIR")
) %>%
  distinct(SYMBOL, TAIR) %>%
  filter(!is.na(TAIR), !is.na(SYMBOL))

label_map

markers_tair <- lapply(markers_sym, map_syms_to_tair, map_tbl = map)

# Report unmapped symbols (useful in tutorial)
mapped_syms <- unique(map$SYMBOL)
unmapped <- setdiff(all_sym, mapped_syms)
if (length(unmapped) > 0) {
  message("Unmapped symbols (no TAIR ID found in org.At.tair.db): ",
          paste(unmapped, collapse = ", "))
}

# --- Handle feature naming quirks in Seurat objects ---
# Some datasets store TAIR IDs as "AT3G11260.1" etc.
# We'll create a "base" version without transcript suffix for matching.
DefaultAssay(seu) <- "RNA"
features <- rownames(seu)
features_base <- sub("\\.\\d+$", "", features)  # strip ".1" ".2" etc if present

# Function to keep only markers present (accounting for possible transcript suffixes)
keep_present_markers <- function(tair_vec) {
  # match on base IDs
  present_base <- tair_vec[tair_vec %in% features_base]
  
  # convert back to actual rownames in the object (prefer exact match; otherwise first transcript)
  present_actual <- vapply(present_base, function(x) {
    if (x %in% features) return(x)
    idx <- which(features_base == x)
    features[idx[1]]
  }, character(1))
  
  unique(present_actual)
}

markers_present <- lapply(markers_tair, keep_present_markers)

# Print missing per group (helpful for documenting snRNA limitations)
for (nm in names(markers_present)) {
  missing_n <- length(markers_tair[[nm]]) - length(markers_present[[nm]])
  message(nm, ": ", length(markers_present[[nm]]), " present; ", missing_n, " missing")
}

Idents(seu) <- "seurat_clusters"  # or "integrated_snn_res.0.6" if you prefer that clustering



p_dot <- DotPlot(
  seu,
  features = markers_present,
  assay = "RNA"
)

# Build named vector: TAIR -> SYMBOL
tair_to_symbol <- label_map$SYMBOL
names(tair_to_symbol) <- label_map$TAIR

p_dot +
  scale_x_discrete(labels = function(x) {
    # Replace TAIR with SYMBOL where possible
    out <- tair_to_symbol[x]
    out[is.na(out)] <- x[is.na(out)]
    out
  }) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 9
    )
  ) +
  labs(
    x = "Marker genes",
    y = "Cluster"
  )



anchors_sym <- c("WOX5","AGL42","BRAVO","PLT1","PLT2",
                 "SCR","SHR","CYCD6;1",
                 "TMO5","LHW","ATHB8",
                 "FEZ","SMB")

anchors_tair <- map_syms_to_tair(anchors_sym, map)
anchors_present <- keep_present_markers(anchors_tair)

pretty_label <- function(tair_ids, label_map) {
  syms <- label_map$SYMBOL[match(tair_ids, label_map$TAIR)]
  ifelse(
    is.na(syms),
    tair_ids,
    paste0(syms, " (", tair_ids, ")")
  )
}

titles <- pretty_label(anchors_present, label_map)

plots <- FeaturePlot(
  seu,
  features = anchors_present, label = TRUE,
  ncol = 4
)

for (i in seq_along(plots)) {
  plots[[i]] <- plots[[i]] + ggtitle(titles[i])
}

plots
DimPlot(seu)

AT1G73590 PIN1
AT3G62100 IAA30
AT5G57090 PIN2
AT1G70940 PIN3
AT1G70560 TAA1
AT1G07370 PCNA1
FeaturePlot(seu, features = c("AT1G07370" ,"AT1G70560", "AT1G70940", "AT5G57090", "AT1G73590", "AT3G62100"), label = TRUE)

new_ids <- c(
  "0" = "QC_CEI_S",
  "1" = "LRCEI",
  "2" = "CSC_LRCEI",
  "3" = "immature vasculature",
  "4" = "vasculature initials",
  "5" = "unknown"
)

seu <- RenameIdents(seu, new_ids)

# Save as a metadata column for plotting and downstream use
seu$cluster_annot <- Idents(seu)

# Plot using the new names
DimPlot(seu, group.by = "cluster_annot", label = TRUE, repel = TRUE) +
  ggtitle("Annotated clusters")

saveRDS(seu, file = "C:/!Victoria/GitHub/iSensors-supplementary/Tutorial-iSensors/data/stemcellniche_reannotated.rds")

