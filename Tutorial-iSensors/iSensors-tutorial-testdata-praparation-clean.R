# iSensors tutorial dataset preparation (Arabidopsis root stem cell niche)
# Source dataset: Martin-Arevalillo et al., Cell (2025); GEO: GSE241573

suppressPackageStartupMessages({
  library(Seurat)
  library(iSensors)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.At.tair.db)
  library(Matrix)
})
set.seed(100)
# ----------------------------
# 1) Load and subset to niche
# ----------------------------

load("C:/!Victoria/!Projects/DigitalSensors/Data/Martin-Arevalillo-2025/GSE241573_Seurat_Object_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13.RData")

seurat_obj <- UpdateSeuratObject(Seurat_Object_Ath_Root_Auxin_treatment_DR5_DR15_IR8_ER13)

# Inspect integration quality (optional)
DimPlot(seurat_obj, reduction = "umap", split.by = "orig.ident2", label = TRUE)

# Subset: stem cell niche / initials
Idents(seurat_obj) <- "Cell.Ident"  # adapt if needed; this is how you subsetted before
seu <- subset(seurat_obj, idents = "Initials_Stele_QC")

# ----------------------------
# 2) Recluster using integrated space (prevents treatment-driven clustering)
# ----------------------------

DefaultAssay(seu) <- "integrated"

# Clear graph-level objects to avoid confusion (keep reductions if present)
seu@graphs <- list()
seu@neighbors <- list()
seu$seurat_clusters <- NULL

# Recompute PCA/UMAP in integrated space
dims_use <- 1:30
seu <- RunPCA(seu, npcs = max(dims_use), verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "pca", dims = dims_use, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(seu, reduction = "pca", dims = dims_use, verbose = FALSE)

# Diagnostics (optional)
DimPlot(seu, group.by = "orig.ident2") + ggtitle("Treatment (orig.ident2) mixing check")
DimPlot(seu, split.by = "orig.ident2") + ggtitle("UMAP split by treatment")

# ----------------------------
# 3) Marker panels (SYMBOLS) + mapping to TAIR IDs
# ----------------------------

markers_sym <- list(
  QC = c("WOX5", "AGL42", "BRAVO", "PLT1", "PLT2"),
  Stele_procambium_initials = c("TMO5", "LHW", "ATHB8", "PIN1"),
  CEI_endodermis_lineage = c("SCR", "SHR", "CYCD6;1"),
  Columella_stem_cells = c("FEZ"),
  Rootcap_LRC_program = c("SMB"),
  Epidermis = c("WER", "CPC", "PIN2"),
  Xylem = c("VND7", "IRX3", "APL"),
  Pericycle_LRP = c("LBD16", "GATA23")
)

all_sym <- unique(unlist(markers_sym))

# Robust mapping: query only needed symbols (avoids occasional DB issues)
tair_vec <- AnnotationDbi::mapIds(
  org.At.tair.db,
  keys = all_sym,
  keytype = "SYMBOL",
  column = "TAIR",
  multiVals = "first"   # if multiple TAIR IDs exist, take the first
)

label_map <- tibble::tibble(
  SYMBOL = names(tair_vec),
  TAIR   = unname(tair_vec)
) %>%
  filter(!is.na(TAIR), !is.na(SYMBOL)) %>%
  distinct(SYMBOL, TAIR)

# Sanity check (print only a small part)
print(head(label_map, 10))


# Helper: map symbol vector -> TAIR vector
map_syms_to_tair <- function(sym_vec, label_map) {
  label_map %>%
    filter(SYMBOL %in% sym_vec) %>%
    arrange(match(SYMBOL, sym_vec)) %>%
    pull(TAIR) %>%
    unique()
}

# Convert marker panels to TAIR IDs
markers_tair <- lapply(markers_sym, map_syms_to_tair, label_map = label_map)

# Report any unmapped symbols
unmapped <- setdiff(all_sym, unique(label_map$SYMBOL))
if (length(unmapped) > 0) {
  message("Unmapped symbols (no TAIR mapping found): ", paste(unmapped, collapse = ", "))
}

# ----------------------------
# 4) Handle transcript suffixes and keep only markers present
# ----------------------------

DefaultAssay(seu) <- "RNA"  # for plotting gene expression
features <- rownames(seu)
features_base <- sub("\\.\\d+$", "", features)  # handle ATxGxxxxx.1

keep_present_markers <- function(tair_vec) {
  present_base <- tair_vec[tair_vec %in% features_base]
  present_actual <- vapply(present_base, function(x) {
    if (x %in% features) return(x)
    idx <- which(features_base == x)
    features[idx[1]]
  }, character(1))
  unique(present_actual)
}

markers_present <- lapply(markers_tair, keep_present_markers)

# ----------------------------
# 5) DotPlot with SYMBOL labels (vertical)
# ----------------------------

Idents(seu) <- "seurat_clusters"

tair_to_symbol <- setNames(label_map$SYMBOL, label_map$TAIR)

p_dot <- DotPlot(seu, features = markers_present, assay = "RNA") +
  scale_x_discrete(labels = function(x) {
    out <- tair_to_symbol[x]
    out[is.na(out)] <- x[is.na(out)]
    out
  }) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9)) +
  labs(x = "Marker genes", y = "Cluster")

p_dot

# ----------------------------
# 6) FeaturePlots with titles "SYMBOL (TAIR)"
# ----------------------------

anchors_sym <- c("WOX5","AGL42","BRAVO","PLT1","PLT2",
                 "SCR","SHR","CYCD6;1",
                 "TMO5","LHW","ATHB8",
                 "FEZ","SMB", "PIN1", "PIN2")

anchors_tair <- map_syms_to_tair(anchors_sym, label_map)
anchors_present <- keep_present_markers(anchors_tair)

pretty_label <- function(tair_ids) {
  syms <- label_map$SYMBOL[match(sub("\\.\\d+$", "", tair_ids), label_map$TAIR)]
  ifelse(is.na(syms), tair_ids, paste0(syms, " (", tair_ids, ")"))
}

titles <- pretty_label(anchors_present)

plots <- FeaturePlot(seu, features = anchors_present, ncol = 4, label = TRUE, order = TRUE)
for (i in seq_along(plots)) plots[[i]] <- plots[[i]] + ggtitle(titles[i])
plots


# ----------------------------
# 7) Rename clusters (store annotation in metadata)
# ----------------------------

# Save numeric clusters before renaming (important for later subsetting)
seu$cluster_id <- as.character(Idents(seu))

new_ids <- c(
  "0" = "Daughters_CEI",
  "1" = "Daughters_LRCEI",
  "2" = "EarlyRootCap",
  "3" = "Daughters_Vasc",
  "4" = "Initials",
  "5" = "TransitInitials"
)

seu <- RenameIdents(seu, new_ids)
seu$cluster_annot <- as.character(Idents(seu))

DimPlot(seu, group.by = "cluster_annot", label = TRUE, repel = TRUE) +
  ggtitle("Annotated clusters")

# ----------------------------
# 8) Export a LIGHT tutorial dataset restricted to clusters 0, 1, 4
#    (counts + umap + meta only; GitHub-friendly)
# ----------------------------

# Subset by the original numeric cluster IDs (stored in seu$cluster_id)
keep_clusters <- c("0", "1", "4", "5")
seu_small <- subset(seu, subset = cluster_id %in% keep_clusters)

# Extract original UMAP
umap <- Embeddings(seu_small, reduction = "umap")

# (optional) Rotate 180°: (x, y) -> (-x, -y) for visualization purposes
umap_rot <- cbind(
  UMAP_1 = -umap[, 1],
  UMAP_2 = -umap[, 2]
)
rownames(umap_rot) <- rownames(umap)

# Replace the existing UMAP reduction
seu_small[["umap"]] <- CreateDimReducObject(
  embeddings = umap_rot,
  key = "UMAP_",
  assay = DefaultAssay(seu_small)
)


p <- DimPlot(
  seu_small,
  reduction = "umap",
  group.by = "cluster_annot",
  repel = TRUE,
  split.by = "orig.ident2"
) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.line        = element_blank(),
    axis.ticks       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    plot.title       = element_blank() 
  ) +
  coord_equal()

p

ggsave(
  "Tutorial-iSensors/out/DimPlot_stemcellniche_auxin.png",
  plot = p,
  width = 10,
  height = 4,
  dpi = 300
)



# Export minimal dataset
counts <- GetAssayData(seu_small, assay = "RNA", layer = "counts")
counts <- as(counts, "dgCMatrix")

umap <- Embeddings(seu_small, "umap")  # rownames are cells
meta <- seu_small@meta.data            # rownames are cells

minimal <- list(
  counts = counts,
  umap = umap,
  meta = meta,
  features = rownames(counts),
  cells = colnames(counts)
)

out_file <- "C:/!Victoria/GitHub/iSensors-supplementary/Tutorial-iSensors/data/stemcellniche_minimal_github.rds"
saveRDS(minimal, file = out_file, compress = "xz")

