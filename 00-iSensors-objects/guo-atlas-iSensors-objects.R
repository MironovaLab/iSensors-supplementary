# Generate iSensors objects and AverageExpression tables for the Guo et al. 2025 atlas
# Run with working directory = iSensors-supplementary/
#
# Input:  00-iSensors-objects/data/Guo/*_rename.slim.RDS  (18 datasets, 16 organs)
# Output per dataset:
#   00-iSensors-objects/data/Guo/<name>_iSensors.rds     (Seurat object with iSensors_mean)
# Combined output:
#   00-iSensors-objects/data/Guo/guo_avg_exp_iSensors_mean.rds
#     — named list; each element is the AverageExpression matrix (sensors × clusters)
#       for that dataset, e.g. guo_avg_exp[["D0_silique"]]

library(iSensors)
library(Seurat)
library(tidyverse)

data_dir <- "00-iSensors-objects/data/Guo"

# ── Load Auxin panel (once) ───────────────────────────────────────────────────
AuxinPanel <- LoadSensors(
  setName      = "Auxin",
  species      = "ATH",
  hormone      = "aux",
  customPanels = FALSE,
  randomInfo   = list(n = 3, sizes = c(100, 200, 300), majortrend = TRUE)
)

# ── Discover all slim RDS files ───────────────────────────────────────────────
rds_files <- list.files(data_dir, pattern = "_rename\\.slim\\.RDS$",
                        full.names = TRUE, ignore.case = TRUE)
message("Found ", length(rds_files), " datasets:")
message(paste(" -", basename(rds_files), collapse = "\n"))

# ── Process each dataset ──────────────────────────────────────────────────────
guo_avg_exp <- list()

for (f in rds_files) {
  dataset_name <- sub("_rename\\.slim\\.RDS$", "", basename(f), ignore.case = TRUE)
  message("\n── Processing: ", dataset_name, " ──")

  obj <- readRDS(f)
  DefaultAssay(obj) <- "RNA"

  # Calculate iSensors_mean if not already present
  if (!"iSensors_mean" %in% Assays(obj)) {

    # Option A (default): use full object — requires ~20+ GB RAM per large dataset
    # obj_for_calc <- obj

    # Option B (low-memory): subset RNA to panel genes before CalcSensors.
    # Avoids dense coercion of the full genome matrix; use if Option A runs out of memory.
    panel_genes  <- unique(unlist(AuxinPanel))
    genes_keep   <- intersect(panel_genes, rownames(obj))
    message("  Subsetting RNA to ", length(genes_keep), " panel genes (of ",
            nrow(obj), " total)")
    obj_for_calc <- subset(obj, features = genes_keep)
    gc()

    message("  Running CalcSensors...")
    obj_for_calc <- CalcSensors(obj_for_calc, seurLayer = "data",
                                panelSet = AuxinPanel, signals = "mean")

    # Transfer iSensors_mean assay back to the full object
    obj[["iSensors_mean"]] <- obj_for_calc[["iSensors_mean"]]
    rm(obj_for_calc); gc()

  } else {
    message("  iSensors_mean already present — skipping CalcSensors")
  }

  # Save updated object
  out_rds <- file.path(data_dir, paste0(dataset_name, "_iSensors.rds"))
  saveRDS(obj, out_rds)
  message("  Saved: ", basename(out_rds))

  # AverageExpression per cluster
  DefaultAssay(obj) <- "iSensors_mean"
  avg <- AverageExpression(obj, assays = "iSensors_mean", slot = "data",
                           verbose = FALSE)[["iSensors_mean"]]
  guo_avg_exp[[dataset_name]] <- avg
  message("  AverageExpression: ", nrow(avg), " sensors × ", ncol(avg), " clusters")
}

# ── Save combined list ────────────────────────────────────────────────────────
combined_path <- file.path(data_dir, "guo_avg_exp_iSensors_mean.rds")
saveRDS(guo_avg_exp, combined_path)
message("\nSaved combined avg_exp list: ", combined_path)
message("Datasets in list: ", paste(names(guo_avg_exp), collapse = ", "))
