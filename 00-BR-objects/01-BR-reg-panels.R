# 01-BR-reg-panels.R
# Creates 3 BR iSensors reg-panels from bulk RNA-Seq DEGs (GSE147589):
#   BL (brassinolide) vs mock at 1hr, 2hr, and 4hr
#
#   ATH-br-reg-1hr-up  / ATH-br-reg-1hr-down
#   ATH-br-reg-2hr-up  / ATH-br-reg-2hr-down
#   ATH-br-reg-4hr-up  / ATH-br-reg-4hr-down
#
# DEG method: limma lmFit on batch-corrected (by replicate) log1p-expression,
# design ~0 + treatment, contrasts: BL - WT per timepoint.
# Significance threshold: adj.P.Val < 0.05, |logFC| > 0.3 (log1p scale: ~1.35-fold)

suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
  library(dplyr)
  library(iSensors)
})

setwd("D:/!GitHub/DigitalSensor-Toolbox/BRiSensors")

# ── 1. Load Seurat object with proper batch correction ───────────────────────
# Stefan's object has treatment/timepoint/replicate columns and a batchcorrected
# assay where batch = replicate (correct approach for DEG analysis).
cat("Loading Seurat object...\n")
obj <- readRDS("StefanReintjes_internship/in/GSE147589_BR_treatment_Seuratobj.rds")

cat("Samples:", ncol(obj), "\n")
cat("Genes:  ", nrow(obj), "\n")
cat("Assays:", paste(Assays(obj), collapse=", "), "\n")

cat("\nTreatment:\n"); print(table(obj$treatment))
cat("\nTimepoint:\n"); print(table(obj$timepoint))
cat("\nReplicate:\n"); print(table(obj$replicate))

# ── 2. Extract batch-corrected expression matrix ─────────────────────────────
# batchcorrected assay: limma::removeBatchEffect with batch = replicate,
# preserving treatment and timepoint signals. Use directly with lmFit.
expr_mat <- GetAssayData(obj, assay = "batchcorrected", layer = "data")
cat("\nUsing batchcorrected assay (batch = replicate correctly removed)\n")
cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")
cat("Value range:", round(range(as.matrix(expr_mat[1:50, ])), 2), "\n")

# ── 3. DEG analysis with limma per timepoint ─────────────────────────────────
# Design: ~0 + treatment  (replicate batch already removed from expression)
timepoints <- c("1hr", "2hr", "4hr")
deg_results <- list()

for (tp in timepoints) {

  cat("\n── Timepoint:", tp, "──────────────────────────────────\n")

  # Subset to this timepoint
  keep <- obj$timepoint == tp
  if (sum(keep) == 0) {
    cat("  No samples found for timepoint:", tp, "- skipping\n")
    next
  }

  expr_tp <- as.matrix(expr_mat[, keep])
  meta_tp  <- obj@meta.data[keep, ]
  cat("  Samples:", ncol(expr_tp), "\n")
  print(table(meta_tp$treatment))

  # Require both BL and WT
  if (!all(c("BL", "WT") %in% meta_tp$treatment)) {
    cat("  Missing treatment group - skipping\n")
    next
  }

  # Design: treatment only (batch-corrected data)
  treatment <- factor(meta_tp$treatment, levels = c("WT", "BL"))
  design    <- model.matrix(~0 + treatment)
  colnames(design) <- c("WT", "BL")
  cat("  Design columns:", paste(colnames(design), collapse = ", "), "\n")

  # Fit linear model
  fit  <- lmFit(expr_tp, design)

  # Contrast: BL - WT
  cont_mat <- makeContrasts(BL - WT, levels = design)
  fit2 <- contrasts.fit(fit, cont_mat)
  fit2 <- eBayes(fit2)

  # All genes, sorted by p-value
  all_degs <- topTable(fit2, coef = 1, n = Inf, sort.by = "p")
  n_sig <- sum(all_degs$adj.P.Val < 0.05 & abs(all_degs$logFC) > 0.3, na.rm = TRUE)
  cat("  Significant DEGs (adj.P<0.05, |logFC|>0.3):", n_sig, "\n")

  sig_degs   <- all_degs[all_degs$adj.P.Val < 0.05 & abs(all_degs$logFC) > 0.3, ]
  up_genes   <- rownames(sig_degs[sig_degs$logFC > 0, ])
  down_genes <- rownames(sig_degs[sig_degs$logFC < 0, ])
  cat("  Up-regulated:", length(up_genes), "\n")
  cat("  Down-regulated:", length(down_genes), "\n")

  deg_results[[tp]] <- list(all = all_degs, up = up_genes, down = down_genes)

  dir.create("out", showWarnings = FALSE)
  write.csv(all_degs,
            file = paste0("out/BR_DEGs_", tp, "_BLvsMock.csv"),
            row.names = TRUE)
}

# ── 4. Report summary ─────────────────────────────────────────────────────────
cat("\n── DEG summary ──────────────────────────────────────────────────────────\n")
for (tp in names(deg_results)) {
  cat(tp, "- up:", length(deg_results[[tp]]$up),
      "  down:", length(deg_results[[tp]]$down), "\n")
}

# ── 5. Create reg-panels ─────────────────────────────────────────────────────
# Each timepoint gets two panels: -up (BL-activated) and -down (BL-repressed)
cat("\n── Creating reg-panels ──────────────────────────────────────────────────\n")

for (tp in names(deg_results)) {
  res <- deg_results[[tp]]

  # Up panel (genes activated by BR at this timepoint)
  if (length(res$up) >= 5) {
    panel_name_up <- paste0("ATH-br-reg-", tp, "-up")
    cat("Creating:", panel_name_up, "(", length(res$up), "genes)\n")
    iSensorsTransPanelCreate(
      panel_name        = panel_name_up,
      species           = "Arabidopsis Thaliana",
      gene_list         = res$up,
      panel_description = paste0("Genes upregulated by BL vs mock at ", tp,
                                 " (adj.P<0.05, logFC>0.3, limma)")
    )
  } else {
    cat("Skipping", paste0("ATH-br-reg-", tp, "-up"),
        ": too few significant genes (", length(res$up), ")\n")
  }

  # Down panel (genes repressed by BR at this timepoint)
  if (length(res$down) >= 5) {
    panel_name_dn <- paste0("ATH-br-reg-", tp, "-down")
    cat("Creating:", panel_name_dn, "(", length(res$down), "genes)\n")
    iSensorsTransPanelCreate(
      panel_name        = panel_name_dn,
      species           = "Arabidopsis Thaliana",
      gene_list         = res$down,
      panel_description = paste0("Genes downregulated by BL vs mock at ", tp,
                                 " (adj.P<0.05, logFC<-0.3, limma)")
    )
  } else {
    cat("Skipping", paste0("ATH-br-reg-", tp, "-down"),
        ": too few significant genes (", length(res$down), ")\n")
  }
}

cat("\nReg-panels created (saved to ./iSensors/). To load all BR panels:\n")
cat("  BRpanels <- LoadSensors(setName='BR', defaultPanels=FALSE, customPanels=TRUE)\n")

# ── 6. Save DEG results for downstream use ───────────────────────────────────
saveRDS(deg_results, "out/BR_DEG_results_limma.rds")
cat("\nDEG results saved to: out/BR_DEG_results_limma.rds\n")
