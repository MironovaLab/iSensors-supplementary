# FigureS_gradient_correlation.R
# Auxin gradient Spearman rank correlation -- benchmark iSensors vs UCell / AUCell
#
# Design replicates 02-single-cell-data-analysis/02-endo-sc-root-spearman.R:
#   1. Per-cell-type average sensor score
#   2. Rank those averages: rank(-means, ties="average")  [rank 1 = highest]
#   3. Spearman rho between predicted rank and literature auxin rank
#   4. Significant: rho > 0.5 AND p < 0.05  (same threshold as Figure 4D)
#
# Ground truth:
#   Literature-based auxin rank for 22 cell types from the GROUNDTRUTH
#   Shahan object (shahan-iSensors-obj-groundtruth.rds), which is the SAME
#   object used in Figure 4D.  Cell type names use dashes as in that object.
#   Rank 1 = highest auxin (QC), rank 12 = lowest (Atrichoblast-m2).
#   Source: Petersson 2009; Brunoud 2012.
#
# Input:
#   00-iSensors-objects/data/shahan-iSensors-obj-groundtruth.rds
#     -- the groundtruth Seurat object (already has iSensors_mean assay)
#
# Outputs:
#   out/scoring_comparison/FigureS_C_gradient_spearman.pdf/.svg
#   out/scoring_comparison/FigureS_C_gradient_extended.pdf/.svg
#   out/scoring_comparison/gradient_spearman.csv
#
# Run from: iSensors-supplementary/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(iSensors)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

for (pkg in c("UCell", "AUCell", "BiocParallel")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({ library(UCell); library(AUCell) })

output_dir <- "Manuscript-Figures/out/scoring_comparison"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Ground-truth auxin rank (Figure 4D / Petersson 2009 / Brunoud 2012)
# 22 cell types; names use underscores to match the groundtruth Seurat object's
# Idents() format.  Rank 1 = QC (highest auxin), rank 12 = Atrichoblast_m2.
# ==============================================================================
auxin_rank <- c(
  "QC"              = 1,
  "CSC"             = 2,
  "LRP"             = 2,
  "Columella"       = 3,
  "CSCD"            = 3,
  "Young LRC"       = 3,
  "Initials"        = 4,
  "Protoxylem_m1"   = 5,
  "LRC"             = 6,
  "Protoxylem_m2"   = 6,
  "Procambium_m1"   = 7,
  "Procambium_m2"   = 7,
  "Endodermis_m1"   = 8,
  "Pericycle_m1"    = 8,
  "Endodermis_m2"   = 9,
  "Pericycle_m2"    = 9,
  "Cortex_m1"       = 10,
  "Cortex_m2"       = 10,
  "Atrichoblast_m1" = 11,
  "Trichoblast_m1"  = 11,
  "Atrichoblast_m2" = 12,
  "Trichoblast_m2"  = 12
)
gt_cell_types <- names(auxin_rank)
message("Ground truth: ", length(gt_cell_types),
        " cell types, ranks 1-", max(auxin_rank))

# ==============================================================================
# Panel setup
# ==============================================================================
message("Loading ATH Auxin panel...")
AuxinPanel <- LoadSensors(
  setName = "Auxin", species = "ATH", hormone = "aux",
  customPanels = FALSE,
  randomInfo   = list(n = 3, sizes = c(100, 200, 300), majortrend = TRUE)
)
panels_list     <- AuxinPanel[["panels"]]
additional_list <- AuxinPanel[["additional"]]
all_panels_raw  <- c(panels_list, additional_list)

safe_key <- function(nm) gsub("[^A-Za-z0-9]", "_", nm)

display_label <- function(nm) {
  nm <- gsub("^ATH-aux-trans-", "", nm)
  nm <- gsub("^ATH-aux-reg-DR5-", "DR5.", nm)
  nm <- gsub("^ATH-aux-reg-IR8-", "IR8.", nm)
  nm <- gsub("^ATH-aux-cis-DR5-", "cis(DR5).", nm)
  nm <- gsub("^ATH-aux-cis-IR8-", "cis(IR8).", nm)
  nm <- gsub("PolarAuxinTransport", "PAT", nm)
  nm <- gsub("ConjugationDeconjugation", "Conjugation", nm)
  nm <- gsub("-up$", "(up)", nm)
  nm <- gsub("-down$", "(dn)", nm)
  nm <- gsub("^random$", "Random", nm)
  nm
}

get_up <- function(el) {
  if (is.list(el) && "up" %in% names(el)) el$up
  else if (is.character(el)) el
  else character(0)
}

rand_el   <- all_panels_raw[["random"]]
rand_flat <- if (is.list(rand_el) && !("up" %in% names(rand_el))) rand_el[[1]] else get_up(rand_el)
if (!is.character(rand_flat)) rand_flat <- character(0)
all_panels_flat <- all_panels_raw
all_panels_flat[["random"]] <- rand_flat

gene_sets_all_raw <- Filter(function(x) length(x) >= 5,
                            lapply(all_panels_flat, get_up))
panel_key_map <- setNames(names(gene_sets_all_raw),
                          vapply(names(gene_sets_all_raw), safe_key, character(1)))
gene_sets_all <- setNames(gene_sets_all_raw, names(panel_key_map))

sel_panel_names <- c(
  ARF         = "ATH-aux-trans-ARF",
  Synthesis   = "ATH-aux-trans-Synthesis",
  Conjugation = "ATH-aux-trans-ConjugationDeconjugation",
  PAT         = "ATH-aux-trans-PolarAuxinTransport",
  IAA         = "ATH-aux-trans-IAA",
  AARF        = "ATH-aux-trans-A-ARF",
  DR5_ARF1_up = "ATH-aux-reg-DR5-ARF1-up",
  DR5_ARF8_dn = "ATH-aux-reg-DR5-ARF8-down",
  DR5_ARF2_dn = "ATH-aux-reg-DR5-ARF2-down",
  Random      = "random"
)
gene_sets_sel <- Filter(function(x) length(x) >= 5,
                        setNames(lapply(sel_panel_names,
                                        function(nm) get_up(all_panels_flat[[nm]])),
                                 names(sel_panel_names)))

sel_display <- setNames(vapply(sel_panel_names, display_label, character(1)),
                        names(sel_panel_names))

method_colors <- c(iSensors = "#2166ac", UCell = "#4dac26", AUCell = "#f4a582")
method_shapes <- c(iSensors = 16, UCell = 15, AUCell = 18)

# ==============================================================================
# Load the groundtruth Seurat object (same as Figure 4D)
# ==============================================================================
gt_rds <- "00-iSensors-objects/data/shahan-iSensors-obj-groundtruth.rds"
if (!file.exists(gt_rds))
  stop("Groundtruth object not found: ", gt_rds)

message("Loading ", gt_rds, " ...")
gt_obj <- readRDS(gt_rds)
message("Loaded: ", ncol(gt_obj), " cells; assays: ",
        paste(Assays(gt_obj), collapse = ", "))
message("Cell types in object: ",
        paste(sort(unique(as.character(Idents(gt_obj)))), collapse = ", "))

# Confirm ground-truth cell types are present
found_ct <- intersect(gt_cell_types, unique(as.character(Idents(gt_obj))))
missing_ct <- setdiff(gt_cell_types, found_ct)
message("Ground-truth cell types matched: ", length(found_ct), " / ", length(gt_cell_types))
if (length(missing_ct) > 0)
  message("  Missing: ", paste(missing_ct, collapse = ", "))

# ==============================================================================
# iSensors: AverageExpression from existing iSensors_mean assay
# Exactly as Figure 4D does it.
# ==============================================================================
DefaultAssay(gt_obj) <- "iSensors_mean"
avg_isensors <- AverageExpression(gt_obj, assays = "iSensors_mean",
                                  slot = "data", verbose = FALSE)[["iSensors_mean"]]
# Seurat v5 converts underscores to dashes in AverageExpression column names; revert.
colnames(avg_isensors) <- gsub("-", "_", colnames(avg_isensors))
message("iSensors avg matrix: ", nrow(avg_isensors), " sensors x ",
        ncol(avg_isensors), " cell types")

# ==============================================================================
# UCell + AUCell scoring on the groundtruth object
# ==============================================================================
DefaultAssay(gt_obj) <- "RNA"

# Selected sensors
n_genes    <- nrow(gt_obj[["RNA"]])
gs_in_sel  <- lapply(gene_sets_sel, function(g) intersect(g, rownames(gt_obj)))
gs_ok_sel  <- gs_in_sel[lengths(gs_in_sel) >= 5]
max_rank_sel <- as.integer(min(n_genes, max(max(lengths(gs_ok_sel)), 1500L)))
message("Scoring selected sensors -- UCell maxRank: ", max_rank_sel)

gt_obj <- suppressWarnings(
  AddModuleScore_UCell(gt_obj, features = gs_ok_sel,
                       maxRank = max_rank_sel,
                       BPPARAM = BiocParallel::SerialParam())
)
for (nm in names(gs_ok_sel)) {
  old <- paste0(nm, "_UCell")
  if (old %in% names(gt_obj@meta.data)) {
    gt_obj@meta.data[[paste0("UCell_", nm)]] <- gt_obj@meta.data[[old]]
    gt_obj@meta.data[[old]] <- NULL
  }
}

expr_counts   <- GetAssayData(gt_obj, assay = "RNA", layer = "counts")
rankings_sel  <- AUCell_buildRankings(expr_counts, nCores = 1,
                                      plotStats = FALSE, verbose = FALSE)
auc_sel       <- AUCell_calcAUC(gs_ok_sel, rankings_sel, verbose = FALSE)
auc_sel_df    <- as.data.frame(t(getAUC(auc_sel)))
names(auc_sel_df) <- paste0("AUCell_", names(auc_sel_df))
gt_obj@meta.data <- cbind(gt_obj@meta.data, auc_sel_df)

# All sensors
gs_in_all   <- lapply(gene_sets_all, function(g) intersect(g, rownames(gt_obj)))
gs_ok_all   <- gs_in_all[lengths(gs_in_all) >= 5]
max_rank_all <- as.integer(min(n_genes, max(max(lengths(gs_ok_all)), 1500L)))
message("Scoring all sensors -- UCell maxRank: ", max_rank_all)

gt_obj <- suppressWarnings(
  AddModuleScore_UCell(gt_obj, features = gs_ok_all,
                       maxRank = max_rank_all,
                       BPPARAM = BiocParallel::SerialParam())
)
for (nm in names(gs_ok_all)) {
  old <- paste0(nm, "_UCell")
  if (old %in% names(gt_obj@meta.data)) {
    gt_obj@meta.data[[paste0("all_UCell_", nm)]] <- gt_obj@meta.data[[old]]
    gt_obj@meta.data[[old]] <- NULL
  }
}

auc_all     <- AUCell_calcAUC(gs_ok_all, rankings_sel, verbose = FALSE)
auc_all_df  <- as.data.frame(t(getAUC(auc_all)))
names(auc_all_df) <- paste0("all_AUCell_", names(auc_all_df))
gt_obj@meta.data <- cbind(gt_obj@meta.data, auc_all_df)

# ==============================================================================
# AverageExpression for UCell / AUCell scores per cell type
# ==============================================================================
ucell_cols  <- grep("^(UCell_|all_UCell_)",   names(gt_obj@meta.data), value = TRUE)
aucell_cols <- grep("^(AUCell_|all_AUCell_)", names(gt_obj@meta.data), value = TRUE)
all_score_cols <- c(ucell_cols, aucell_cols)

ct_labels <- as.character(Idents(gt_obj))
ct_means_ucell <- do.call(rbind, lapply(found_ct, function(ct) {
  idx  <- which(ct_labels == ct)
  vals <- colMeans(gt_obj@meta.data[idx, all_score_cols, drop = FALSE],
                   na.rm = TRUE)
  c(cell_type = ct, vals)
}))
ct_means_ucell <- as.data.frame(ct_means_ucell, stringsAsFactors = FALSE)
for (col in all_score_cols)
  ct_means_ucell[[col]] <- as.numeric(ct_means_ucell[[col]])

rm(gt_obj); gc()

# ==============================================================================
# Spearman rho computation
# ==============================================================================
rank_vec <- auxin_rank[found_ct]

compute_spearman_vec <- function(vals, rk) {
  ok <- !is.na(vals) & !is.na(rk)
  if (sum(ok) < 5) return(c(rho = NA_real_, pval = NA_real_))
  vals_ranked <- rank(-vals[ok], ties.method = "average")
  tst <- tryCatch(
    cor.test(vals_ranked, rk[ok], method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(tst)) return(c(rho = NA_real_, pval = NA_real_))
  c(rho = unname(tst$estimate), pval = tst$p.value)
}

# iSensors: from avg_isensors matrix (already averaged per cell type)
isensors_cols <- rownames(avg_isensors)
avg_isensors_gt <- avg_isensors[, found_ct, drop = FALSE]

spear_isensors <- do.call(rbind, lapply(isensors_cols, function(feat) {
  vals <- as.numeric(avg_isensors_gt[feat, ])
  res  <- compute_spearman_vec(vals, rank_vec)
  data.frame(col = feat, rho = res["rho"], pval = res["pval"],
             Method = "iSensors", SensorKey = feat,
             Pool = if (feat %in% sel_panel_names) "sel" else "all",
             stringsAsFactors = FALSE)
}))

# UCell and AUCell: from ct_means_ucell
parse_ucell_col <- function(col) {
  if      (grepl("^all_UCell_",   col)) list(method = "UCell",   sk = sub("^all_UCell_",   "", col), pool = "all")
  else if (grepl("^all_AUCell_",  col)) list(method = "AUCell",  sk = sub("^all_AUCell_",  "", col), pool = "all")
  else if (grepl("^UCell_",       col)) list(method = "UCell",   sk = sub("^UCell_",       "", col), pool = "sel")
  else if (grepl("^AUCell_",      col)) list(method = "AUCell",  sk = sub("^AUCell_",      "", col), pool = "sel")
  else                                  list(method = NA, sk = col, pool = NA)
}

spear_uca <- do.call(rbind, lapply(all_score_cols, function(col) {
  p    <- parse_ucell_col(col)
  vals <- as.numeric(ct_means_ucell[[col]])
  res  <- compute_spearman_vec(vals, rank_vec)
  data.frame(col = col, rho = res["rho"], pval = res["pval"],
             Method = p$method, SensorKey = p$sk, Pool = p$pool,
             stringsAsFactors = FALSE)
}))

# Combine; rho sign is already correct (rank 1 = highest both for sensor and auxin)
spear_all_raw <- bind_rows(spear_isensors, spear_uca) %>%
  filter(!is.na(Method)) %>%
  mutate(significant = !is.na(rho) & rho > 0.5 & pval < 0.05)

# ==============================================================================
# Selected-sensor table
# ==============================================================================
# iSensors: SensorKey = full panel name e.g. "ATH-aux-trans-ARF"
# UCell/AUCell: SensorKey = short label e.g. "ARF" (from gene_sets_sel names)
spear_sel <- spear_all_raw %>%
  filter(Pool == "sel") %>%
  mutate(
    ShortKey = case_when(
      Method == "iSensors" ~ names(sel_panel_names)[match(SensorKey, sel_panel_names)],
      SensorKey %in% names(sel_panel_names) ~ SensorKey,
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(ShortKey)) %>%
  mutate(
    Method       = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
    DisplayLabel = sel_display[ShortKey],
    DisplayLabel = factor(DisplayLabel,
                          levels = rev(sel_display[names(sel_panel_names)]))
  )

# All-sensor table for the extended figure.
# iSensors: include ALL sensors (sel + all pools) so selected sensors (ARF, PAT…)
#   appear alongside UCell/AUCell.
# UCell/AUCell: use only "all" pool (which covers every panel including selected ones).
spear_all_clean <- bind_rows(
  spear_all_raw %>% filter(Method == "iSensors"),
  spear_all_raw %>% filter(Method != "iSensors", Pool == "all")
) %>%
  mutate(
    Method = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
    OrigKey = case_when(
      Method == "iSensors" ~ SensorKey,
      !is.na(panel_key_map[SensorKey]) ~ panel_key_map[SensorKey],
      TRUE ~ SensorKey
    ),
    DisplayLabel = display_label(OrigKey)
  ) %>%
  group_by(DisplayLabel) %>%
  mutate(isens_rho = { r <- rho[Method == "iSensors"]; if (length(r) > 0) r[1] else NA_real_ }) %>%
  ungroup() %>%
  arrange(isens_rho) %>%
  mutate(DisplayLabel = factor(DisplayLabel, levels = unique(DisplayLabel)))

# Print
message("\nSpearman rho for selected sensors (rho > 0 = tracks auxin gradient):")
wide <- spear_sel %>%
  select(ShortKey, Method, rho, significant) %>%
  pivot_wider(names_from = Method, values_from = c(rho, significant))
print(as.data.frame(wide), digits = 3)

# Save
write.csv(
  bind_rows(
    spear_sel       %>% select(Method, ShortKey, DisplayLabel,
                               Spearman_rho = rho, pval, significant),
    spear_all_clean %>% select(Method, SensorKey, DisplayLabel,
                               Spearman_rho = rho, pval, significant)
  ),
  file.path(output_dir, "gradient_spearman.csv"), row.names = FALSE
)

# ==============================================================================
# Figures
# ==============================================================================
make_spear_plot <- function(df, rho_col, title, x_lim = c(-0.8, 1.05)) {
  df[[rho_col]] <- as.numeric(df[[rho_col]])
  ggplot(df, aes(x = .data[[rho_col]], y = DisplayLabel,
                 color = Method, shape = Method)) +
    geom_vline(xintercept = 0.5,
               linetype = "dashed", color = "grey30", linewidth = 0.8) +
    geom_point(size = 2.8, alpha = 0.9,
               position = position_dodge(width = 0.6)) +
    scale_color_manual(values = method_colors) +
    scale_shape_manual(values = method_shapes) +
    scale_x_continuous(limits = x_lim, breaks = seq(-1, 1, 0.25)) +
    labs(x = paste0("Spearman rho (ground-truth auxin rank, ",
                    length(found_ct), " cell types)"),
         y = NULL, color = NULL, shape = NULL, title = title) +
    theme_classic(base_size = 10) +
    theme(plot.title         = element_text(face = "bold", size = 10),
          legend.position    = "bottom",
          legend.key.size    = unit(0.4, "cm"),
          axis.line          = element_line(linewidth = 0.4),
          panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88"))
}

# Scatter: ARF sensor score rank vs auxin rank, one facet per method
arf_isensor_vals <- as.numeric(avg_isensors_gt["ATH-aux-trans-ARF", found_ct])
ucell_arf_vals   <- as.numeric(ct_means_ucell[["UCell_ARF"]])
aucell_arf_vals  <- as.numeric(ct_means_ucell[["AUCell_ARF"]])

scatter_df <- data.frame(
  cell_type  = rep(found_ct, 3),
  score      = c(arf_isensor_vals, ucell_arf_vals, aucell_arf_vals),
  Method     = factor(rep(c("iSensors", "UCell", "AUCell"), each = length(found_ct)),
                      levels = c("iSensors", "UCell", "AUCell")),
  auxin_rank = rep(rank_vec, 3)
)

arf_rhos <- spear_sel %>%
  filter(ShortKey == "ARF") %>%
  mutate(label = paste0("rho = ", round(rho, 2),
                        ifelse(significant, "*", "")))

p_scatter <- ggplot(scatter_df,
                    aes(x = auxin_rank, y = score, color = Method)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text(data = arf_rhos,
            aes(label = label),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4,
            size = 3, inherit.aes = FALSE, color = "grey30") +
  scale_color_manual(values = method_colors) +
  scale_x_reverse(breaks = c(1, 3, 5, 7, 9, 11),
                  labels = c("1\n(QC)", "3", "5", "7", "9", "11\n(Epidermis)")) +
  facet_wrap(~ Method, nrow = 1, scales = "free_y") +
  labs(x = "Ground-truth auxin rank (1 = highest)",
       y = "ARF sensor score\n(per cell type mean)",
       title = paste0("ARF sensor vs ground-truth auxin rank (",
                      length(found_ct), " cell types)")) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 9),
        legend.position  = "none",
        axis.line        = element_line(linewidth = 0.4))

p_C <- make_spear_plot(
  spear_sel, "rho",
  "Auxin gradient Spearman rho -- selected sensors\n(dotted line = significance threshold rho > 0.5; same as Figure 4D)"
)

fig_C <- (p_scatter / p_C) +
  plot_layout(heights = c(1, 2.4)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 13))
  ) & theme(plot.margin = margin(4, 8, 4, 8, "pt"))

ggsave(file.path(output_dir, "FigureS_C_gradient_spearman.pdf"),
       plot = fig_C, width = 9, height = 9, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "FigureS_C_gradient_spearman.svg"),
       plot = fig_C, width = 9, height = 9, dpi = 300, bg = "white")
message("Saved FigureS_C_gradient_spearman")

# Extended: all sensors
p_C_ext <- make_spear_plot(
  spear_all_clean, "rho",
  "Auxin gradient Spearman rho -- all sensors (groundtruth object)"
)
h_ext <- max(6, 0.22 * length(unique(spear_all_clean$DisplayLabel)))
ggsave(file.path(output_dir, "FigureS_C_gradient_extended.pdf"),
       plot = p_C_ext, width = 7, height = h_ext, dpi = 300, bg = "white",
       limitsize = FALSE)
ggsave(file.path(output_dir, "FigureS_C_gradient_extended.svg"),
       plot = p_C_ext, width = 7, height = h_ext, dpi = 300, bg = "white",
       limitsize = FALSE)
message("Saved FigureS_C_gradient_extended (", round(h_ext, 1), " in tall)")
message("Done.")
