# FigureS_scoring_comparison.R
# -----------------------------------------------------------------------------
# Experimental design:
#   Goal: Position iSensors relative to standard gene-set scoring methods
#         (UCell and AUCell) using identical gene sets applied to the same datasets.
#
#   Evaluation datasets:
#     1. Martin-Arevalillo 2025 (GSE241573) -- exogenous auxin treatment
#        Ground truth: orig.ident2 = "Auxin" vs "Control"
#     2. Shahan 2022 root tip -- endogenous auxin gradient
#        Ground truth: QC / Columella / CSC / CSCD / Young LRC / LRCEI / SI / CEI
#        = high-auxin zone (DR5 reporter; Benkova 2003; Blilou 2005)
#
#   Methods: iSensors | UCell | AUCell
#   Gene sets: ATH Auxin panel (up-genes for UCell/AUCell; full up+down for iSensors)
#   Metric: AUROC per (method, sensor, dataset)
#
#   Outputs:
#     FigureS_scoring_comparison.pdf/.svg      -- Panel A (violins) + Panel B (AUROC, selected sensors)
#     FigureS_B_extended_D1.pdf/.svg           -- AUROC for ALL sensors, Dataset 1
#     FigureS_B_extended_D2.pdf/.svg           -- AUROC for ALL sensors, Dataset 2
#
# Run with working directory = iSensors-supplementary/
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
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " ...")
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
suppressPackageStartupMessages({ library(UCell); library(AUCell) })

output_dir <- "Manuscript-Figures/out/scoring_comparison"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -- AUROC utility -------------------------------------------------------------
auroc <- function(scores, is_positive) {
  n1 <- sum(is_positive, na.rm = TRUE)
  n2 <- sum(!is_positive, na.rm = TRUE)
  if (n1 < 3 || n2 < 3) return(NA_real_)
  valid <- !is.na(scores) & !is.na(is_positive)
  w <- suppressWarnings(
    wilcox.test(scores[valid & is_positive], scores[valid & !is_positive],
                exact = FALSE, correct = FALSE)$statistic
  )
  unname(w / (n1 * n2))
}


# ==============================================================================
# Section 1 -- ATH Auxin panel: build gene set collections
# ==============================================================================
message("Loading ATH Auxin panel...")
AuxinPanel <- LoadSensors(
  setName      = "Auxin", species = "ATH", hormone = "aux",
  customPanels = FALSE,
  randomInfo   = list(n = 3, sizes = c(100, 200, 300), majortrend = TRUE)
)

panels_list    <- AuxinPanel[["panels"]]
additional_list <- AuxinPanel[["additional"]]
all_panels_raw <- c(panels_list, additional_list)
message("Total panel entries: ", length(all_panels_raw))

get_up <- function(el) {
  if (is.list(el) && "up" %in% names(el)) el$up
  else if (is.character(el)) el
  else character(0)
}

# Safe column-name key: replaces all non-alphanumeric chars with "_"
safe_key <- function(nm) gsub("[^A-Za-z0-9]", "_", nm)

# Display label for plots: strip ATH-aux-*- prefix, prettify up/down
display_label <- function(nm) {
  nm <- gsub("^ATH-aux-trans-", "", nm)
  nm <- gsub("^ATH-aux-reg-DR5-", "DR5.", nm)  # .
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

# -- Build FULL gene set collection (all panels) -------------------------------
# random in additional_list is nested; flatten it
rand_el <- all_panels_raw[["random"]]
rand_flat <- if (is.list(rand_el) && !("up" %in% names(rand_el))) {
  rand_el[[1]]
} else {
  get_up(rand_el)
}
if (!is.character(rand_flat)) rand_flat <- character(0)

# Replace the nested random entry with the flattened version
all_panels_flat <- all_panels_raw
all_panels_flat[["random"]] <- rand_flat

gene_sets_all_raw <- lapply(all_panels_flat, get_up)
# Keep sensors with >= 5 up-genes (filtered per-dataset later)
gene_sets_all_raw <- Filter(function(x) length(x) >= 5, gene_sets_all_raw)
message("Panels with >= 5 up-genes: ", length(gene_sets_all_raw))

# Rename keys to safe R names; keep mapping to original panel names
panel_key_map <- setNames(names(gene_sets_all_raw),
                          vapply(names(gene_sets_all_raw), safe_key, character(1)))
gene_sets_all <- setNames(gene_sets_all_raw, names(panel_key_map))  # safe names

# -- Build SELECTED gene set collection (for main Figure B) -------------------
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
# Display labels for selected sensors
sel_display <- vapply(sel_panel_names, display_label, character(1))
names(sel_display) <- names(sel_panel_names)
message("Selected sensors: ", paste(names(sel_panel_names), collapse = ", "))

gene_sets_sel <- Filter(function(x) length(x) >= 5,
                        setNames(lapply(sel_panel_names,
                                        function(nm) get_up(all_panels_flat[[nm]])),
                                 names(sel_panel_names)))

# iSensors assay row names (dashes, as stored)
isensor_row_names <- sel_panel_names  # panel key IS the assay row name


# ==============================================================================
# Section 2 -- Helper functions
# ==============================================================================

# Extract ALL iSensor scores from iSensors_mean assay into metadata
add_all_isensor_meta <- function(obj) {
  if (!"iSensors_mean" %in% Assays(obj)) {
    warning("No iSensors_mean assay in object -- iSensors scores will be NA")
    return(obj)
  }
  # Single GetAssayData call -- avoids per-feature FetchData overhead and
  # handles dashes in feature names without backtick workarounds
  ad    <- GetAssayData(obj, assay = "iSensors_mean", layer = "data")
  avail <- rownames(ad)
  message("    iSensors_mean has ", length(avail), " sensors")

  # Selected sensors
  for (label in names(isensor_row_names)) {
    row_nm <- isensor_row_names[[label]]
    obj@meta.data[[paste0("iSensors_", label)]] <-
      if (row_nm %in% avail) as.numeric(ad[row_nm, ]) else NA_real_
  }
  # All panel sensors (safe key names as column suffix)
  for (sk in names(panel_key_map)) {
    row_nm <- panel_key_map[[sk]]
    obj@meta.data[[paste0("iSensors_all_", sk)]] <-
      if (row_nm %in% avail) as.numeric(ad[row_nm, ]) else NA_real_
  }
  DefaultAssay(obj) <- "RNA"
  obj
}

# Run UCell + AUCell on a combined gene set list; prefix result columns
score_ucell_aucell <- function(obj, gs_named, col_prefix = "", label = "") {
  message("\n  Scoring ", label, ": UCell + AUCell on ", length(gs_named),
          " gene sets...")
  DefaultAssay(obj) <- "RNA"
  gs_in <- lapply(gs_named, function(g) intersect(g, rownames(obj)))
  gs_ok <- gs_in[lengths(gs_in) >= 5]
  message("    Gene sets passing filter: ", length(gs_ok))

  # UCell -- maxRank must be >= largest gene set AND >= 1500 (standard minimum)
  max_rank <- as.integer(min(nrow(obj[["RNA"]]),
                             max(max(lengths(gs_ok)), 1500L)))
  message("    UCell maxRank: ", max_rank)
  obj <- suppressWarnings(
    AddModuleScore_UCell(obj, features = gs_ok, maxRank = max_rank,
                         BPPARAM = BiocParallel::SerialParam())
  )
  for (nm in names(gs_ok)) {
    old <- paste0(nm, "_UCell")
    if (old %in% names(obj@meta.data)) {
      obj@meta.data[[paste0(col_prefix, "UCell_", nm)]] <- obj@meta.data[[old]]
      obj@meta.data[[old]] <- NULL
    }
  }

  # AUCell
  expr_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  rankings    <- AUCell_buildRankings(expr_counts, nCores = 1,
                                      plotStats = FALSE, verbose = FALSE)
  auc_res  <- AUCell_calcAUC(gs_ok, rankings, verbose = FALSE)
  auc_vals <- as.data.frame(t(getAUC(auc_res)))
  names(auc_vals) <- paste0(col_prefix, "AUCell_", names(auc_vals))
  obj@meta.data <- cbind(obj@meta.data, auc_vals)

  obj
}

# AUROC table for one dataset
make_auroc_table <- function(meta, is_positive, sensor_keys, col_prefix = "") {
  methods <- c("iSensors", "UCell", "AUCell")
  results <- expand.grid(SensorKey = sensor_keys, Method = methods,
                         stringsAsFactors = FALSE)
  results$AUROC <- mapply(function(sk, method) {
    col <- if (method == "iSensors") paste0("iSensors_", sk)
           else paste0(col_prefix, method, "_", sk)
    if (!col %in% names(meta)) return(NA_real_)
    auroc(meta[[col]], is_positive)
  }, results$SensorKey, results$Method)
  results
}

# AUROC table for ALL panels
make_auroc_table_all <- function(meta, is_positive) {
  methods <- c("iSensors", "UCell", "AUCell")
  all_keys <- names(panel_key_map)  # safe key names
  results  <- expand.grid(SensorKey = all_keys, Method = methods,
                          stringsAsFactors = FALSE)
  results$AUROC <- mapply(function(sk, method) {
    col <- if (method == "iSensors") paste0("iSensors_all_", sk)
           else paste0("all_", method, "_", sk)
    if (!col %in% names(meta)) return(NA_real_)
    auroc(meta[[col]], is_positive)
  }, results$SensorKey, results$Method)
  results$DisplayLabel <- display_label(panel_key_map[results$SensorKey])
  results
}


# ==============================================================================
# Section 3 -- Dataset 1: Martin-Arevalillo 2025 (exogenous auxin)
# ==============================================================================
message("\n== Dataset 1: Martin-Arevalillo auxin treatment ==")

auxin_obj <- readRDS(
  "00-iSensors-objects/data/iSensors-Martin-Arevalillo-auxin-root2025_mean.rds"
)
message("Cells: ", ncol(auxin_obj), "  treatment: ",
        paste(names(table(auxin_obj$orig.ident2)),
              table(auxin_obj$orig.ident2), sep = "=", collapse = ", "))

auxin_obj <- add_all_isensor_meta(auxin_obj)

# Score with SELECTED gene sets
auxin_obj <- score_ucell_aucell(auxin_obj, gene_sets_sel,
                                col_prefix = "", label = "selected sensors")
# Score with ALL gene sets (col_prefix "all_" to avoid collision)
auxin_obj <- score_ucell_aucell(auxin_obj, gene_sets_all,
                                col_prefix = "all_", label = "all sensors")

is_auxin <- auxin_obj$orig.ident2 == "Auxin"

auroc_sel_d1 <- make_auroc_table(auxin_obj@meta.data, is_auxin,
                                 names(gene_sets_sel))
auroc_sel_d1$Dataset <- "Exogenous auxin\n(Martin-Arevalillo 2025)"

auroc_all_d1 <- make_auroc_table_all(auxin_obj@meta.data, is_auxin)
auroc_all_d1$Dataset <- "Exogenous auxin\n(Martin-Arevalillo 2025)"

# Violin data for Panel A (ARF sensor, 3 methods)
violin_meta <- auxin_obj@meta.data %>%
  select(Treatment = orig.ident2,
         iSensors_ARF,
         UCell_ARF,
         AUCell_ARF)

saveRDS(auxin_obj@meta.data, file.path(output_dir, "auxin_meta_scores.rds"))
rm(auxin_obj); gc()


# ==============================================================================
# Section 4 -- Dataset 2: Shahan root tip (endogenous gradient)
# ==============================================================================
message("\n== Dataset 2: Shahan root tip gradient ==")

shahan_obj <- readRDS("00-iSensors-objects/data/shahan-roottip-iSensors-obj.rds")
if ("root_tip" %in% names(shahan_obj@meta.data))
  Idents(shahan_obj) <- shahan_obj$root_tip

ct <- as.character(Idents(shahan_obj))
message("Cell types: ", paste(sort(unique(ct)), collapse = ", "))

# High-auxin zone includes initials by their original labels (LRCEI, SI, CEI)
high_auxin_types <- c("QC", "CSC", "Columella", "CSCD", "Young LRC",
                      "Initials", "LRP", "LRCEI", "SI", "CEI")
is_high_auxin <- ct %in% high_auxin_types
message("High-auxin cells: ", sum(is_high_auxin), " / ", length(is_high_auxin),
        "  (", paste(sort(unique(ct[is_high_auxin])), collapse = ", "), ")")

shahan_obj <- add_all_isensor_meta(shahan_obj)
shahan_obj <- score_ucell_aucell(shahan_obj, gene_sets_sel,
                                 col_prefix = "", label = "selected sensors")
shahan_obj <- score_ucell_aucell(shahan_obj, gene_sets_all,
                                 col_prefix = "all_", label = "all sensors")

auroc_sel_d2 <- make_auroc_table(shahan_obj@meta.data, is_high_auxin,
                                 names(gene_sets_sel))
auroc_sel_d2$Dataset <- "Endogenous gradient\n(Shahan root tip)"

auroc_all_d2 <- make_auroc_table_all(shahan_obj@meta.data, is_high_auxin)
auroc_all_d2$Dataset <- "Endogenous gradient\n(Shahan root tip)"

shahan_obj$cell_type_gradient <- ct  # preserve for gradient analysis
saveRDS(shahan_obj@meta.data, file.path(output_dir, "shahan_meta_scores.rds"))
rm(shahan_obj); gc()


# ==============================================================================
# Section 5 -- Figures
# ==============================================================================

method_colors <- c(iSensors = "#2166ac", UCell = "#4dac26", AUCell = "#f4a582")
method_shapes <- c(iSensors = 16,        UCell = 15,        AUCell = 18)

# -- AUROC dotplot helper ------------------------------------------------------
make_auroc_plot <- function(df, title = "Discrimination performance (AUROC)",
                            x_lim = c(0.25, 1.01)) {
  ggplot(df, aes(x = AUROC, y = DisplayLabel, color = Method, shape = Method)) +
    geom_vline(xintercept = 0.5, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_point(size = 2.8, alpha = 0.9,
               position = position_dodge(width = 0.6)) +
    scale_color_manual(values = method_colors) +
    scale_shape_manual(values = method_shapes) +
    scale_x_continuous(limits = x_lim,
                       breaks = seq(0.3, 1.0, 0.1)) +
    facet_wrap(~ Dataset, ncol = 1) +
    labs(x = "AUROC", y = NULL, color = NULL, shape = NULL, title = title) +
    theme_classic(base_size = 10) +
    theme(strip.background   = element_blank(),
          strip.text         = element_text(face = "bold", size = 9),
          plot.title         = element_text(face = "bold", size = 10),
          legend.position    = "bottom",
          legend.key.size    = unit(0.4, "cm"),
          axis.line          = element_line(linewidth = 0.4),
          panel.grid.major.x = element_line(linewidth = 0.2, color = "grey88"))
}

# -- Panel A: Violins -- ARF sensor, exogenous auxin dataset -------------------
violin_df <- violin_meta %>%
  pivot_longer(-Treatment, names_to = "raw_col", values_to = "Score") %>%
  mutate(
    Method    = recode(raw_col,
                       iSensors_ARF = "iSensors",
                       UCell_ARF    = "UCell",
                       AUCell_ARF   = "AUCell"),
    Method    = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
    Treatment = factor(Treatment, levels = c("Control", "Auxin"))
  ) %>%
  filter(!is.na(Score), !is.na(Treatment))

arf_auroc_annot <- auroc_sel_d1 %>%
  filter(SensorKey == "ARF") %>%
  mutate(Method = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
         label  = paste0("AUROC=", round(AUROC, 2)))

p_A <- ggplot(violin_df, aes(x = Treatment, y = Score, fill = Treatment)) +
  geom_violin(trim = TRUE, scale = "width", linewidth = 0.3, alpha = 0.75) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white",
               linewidth = 0.35, coef = 0) +
  geom_text(data = arf_auroc_annot,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.4, size = 2.8, colour = "grey30") +
  scale_fill_manual(values = c(Control = "#b2b2b2", Auxin = "#e34a33"),
                    guide = "none") +
  facet_wrap(~ Method, nrow = 1, scales = "free_y") +
  labs(x = NULL, y = "ARF sensor score",
       title = "ARF sensor -- exogenous auxin response") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 9),
        plot.title       = element_text(face = "bold", size = 10),
        axis.line        = element_line(linewidth = 0.4),
        axis.ticks       = element_line(linewidth = 0.3),
        axis.text.x      = element_text(angle = 30, hjust = 1))

# -- Panel B: AUROC dotplot -- SELECTED sensors x both datasets ----------------
auroc_sel_both <- bind_rows(auroc_sel_d1, auroc_sel_d2) %>%
  mutate(
    Method = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
    DisplayLabel = sel_display[SensorKey],
    DisplayLabel = factor(DisplayLabel,
                          levels = rev(sel_display[names(gene_sets_sel)])),
    Dataset = factor(Dataset,
                     levels = c("Exogenous auxin\n(Martin-Arevalillo 2025)",
                                "Endogenous gradient\n(Shahan root tip)"))
  )

write.csv(auroc_sel_both, file.path(output_dir, "auroc_selected_comparison.csv"),
          row.names = FALSE)

p_B <- make_auroc_plot(auroc_sel_both,
                       title = "Discrimination performance (AUROC) -- selected sensors")

# -- Combined main figure ------------------------------------------------------
fig_combined <- (p_A / p_B) +
  plot_layout(heights = c(1, 2.5)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 13))
  ) & theme(plot.margin = margin(4, 8, 4, 8, "pt"))

ggsave(file.path(output_dir, "FigureS_scoring_comparison.pdf"),
       plot = fig_combined, width = 9, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "FigureS_scoring_comparison.svg"),
       plot = fig_combined, width = 9, height = 11, dpi = 300, bg = "white")
message("Saved FigureS_scoring_comparison")

# -- Extended Figure B -- ALL sensors, Dataset 1 -------------------------------
auroc_all_d1_plot <- auroc_all_d1 %>%
  filter(!is.na(AUROC)) %>%
  mutate(Method = factor(Method, levels = c("iSensors", "UCell", "AUCell"))) %>%
  group_by(SensorKey) %>%
  mutate(isens_auroc = AUROC[Method == "iSensors"][1]) %>%
  ungroup() %>%
  arrange(isens_auroc) %>%
  mutate(DisplayLabel = factor(DisplayLabel, levels = unique(DisplayLabel)))

p_B_ext_D1 <- make_auroc_plot(
  auroc_all_d1_plot,
  title = "All sensors -- Exogenous auxin (Martin-Arevalillo 2025)",
  x_lim = c(0.2, 1.01)
)

h_d1 <- max(5, 0.22 * length(unique(auroc_all_d1_plot$DisplayLabel)))
ggsave(file.path(output_dir, "FigureS_B_extended_D1.pdf"),
       plot = p_B_ext_D1, width = 7, height = h_d1, dpi = 300, bg = "white",
       limitsize = FALSE)
ggsave(file.path(output_dir, "FigureS_B_extended_D1.svg"),
       plot = p_B_ext_D1, width = 7, height = h_d1, dpi = 300, bg = "white",
       limitsize = FALSE)
message("Saved FigureS_B_extended_D1 (", round(h_d1, 1), " in tall)")

# -- Extended Figure B -- ALL sensors, Dataset 2 -------------------------------
auroc_all_d2_plot <- auroc_all_d2 %>%
  filter(!is.na(AUROC)) %>%
  mutate(Method = factor(Method, levels = c("iSensors", "UCell", "AUCell"))) %>%
  group_by(SensorKey) %>%
  mutate(isens_auroc = AUROC[Method == "iSensors"][1]) %>%
  ungroup() %>%
  arrange(isens_auroc) %>%
  mutate(DisplayLabel = factor(DisplayLabel, levels = unique(DisplayLabel)))

p_B_ext_D2 <- make_auroc_plot(
  auroc_all_d2_plot,
  title = "All sensors -- Endogenous gradient (Shahan root tip)",
  x_lim = c(0.2, 1.01)
)

h_d2 <- max(5, 0.22 * length(unique(auroc_all_d2_plot$DisplayLabel)))
ggsave(file.path(output_dir, "FigureS_B_extended_D2.pdf"),
       plot = p_B_ext_D2, width = 7, height = h_d2, dpi = 300, bg = "white",
       limitsize = FALSE)
ggsave(file.path(output_dir, "FigureS_B_extended_D2.svg"),
       plot = p_B_ext_D2, width = 7, height = h_d2, dpi = 300, bg = "white",
       limitsize = FALSE)
message("Saved FigureS_B_extended_D2 (", round(h_d2, 1), " in tall)")

write.csv(bind_rows(auroc_all_d1, auroc_all_d2),
          file.path(output_dir, "auroc_all_panels.csv"), row.names = FALSE)


# ==============================================================================
# Section 6 -- Auxin gradient correlation (Shahan root tip, cell-type level)
#
# Rationale: AUROC tests binary discrimination (high-auxin zone vs rest),
# which can be inflated because QC/columella cells are transcriptionally
# distinct for many reasons beyond auxin. Spearman rho across all cell types
# is a more stringent test: the sensor must track the SMOOTH gradient from
# very high auxin (QC, columella) through intermediate zones to low-auxin
# differentiated types.
#
# Reference gradient: iSensors DR5-ARF1-up score, averaged per cell type.
#   DR5 promoter gene set = canonical synthetic auxin reporter (independent
#   of the sensors being benchmarked: ARF, IAA, PAT, Synthesis, A-ARF).
#
# Outputs:
#   FigureS_C_gradient_spearman.pdf/.svg    -- selected sensors
#   FigureS_C_gradient_extended.pdf/.svg   -- all 57 sensors
#   gradient_spearman.csv
# ==============================================================================
message("\n== Section 6: gradient correlation (cell-type level) ==")

shahan_meta <- readRDS(file.path(output_dir, "shahan_meta_scores.rds"))

# Detect cell type column
ct_col <- intersect(c("cell_type_gradient", "root_tip", "seurat_clusters"),
                    names(shahan_meta))
if (length(ct_col) == 0) stop("No cell type column found in shahan metadata")
ct_col <- ct_col[1]
message("Cell type column: ", ct_col)

# All score columns (iSensors, UCell, AUCell for selected + all panels)
score_cols <- names(shahan_meta)[
  grepl("^(iSensors|UCell|AUCell|all_UCell|all_AUCell|iSensors_all)_",
        names(shahan_meta))
]
message("Score columns: ", length(score_cols))

# Compute mean score per cell type
ct_means <- shahan_meta %>%
  rename(cell_type = all_of(ct_col)) %>%
  group_by(cell_type) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE),
            n_cells = n(), .groups = "drop") %>%
  filter(n_cells >= 5)
message("Cell types retained (>= 5 cells): ", nrow(ct_means))

# Auxin gradient reference: iSensors DR5-ARF1-up per cell type
ref_col <- "iSensors_DR5_ARF1_up"
if (!ref_col %in% names(ct_means)) {
  stop("Reference column ", ref_col, " not found -- check sensor was computed")
}
auxin_ref <- ct_means[[ref_col]]
message("Reference range: ", round(min(auxin_ref, na.rm=TRUE), 3), " - ",
        round(max(auxin_ref, na.rm=TRUE), 3))

# Compute Spearman rho for every score column vs the reference
spear_all <- data.frame(
  col = score_cols,
  rho = vapply(score_cols, function(col) {
    cor(auxin_ref, ct_means[[col]], method = "spearman", use = "complete.obs")
  }, numeric(1)),
  stringsAsFactors = FALSE
)

# Parse method and sensor key from column name
# Selected sensors: "iSensors_ARF", "UCell_ARF", "AUCell_ARF"
# All sensors:      "iSensors_all_ATH_aux_trans_ARF", "all_UCell_ATH_aux_trans_ARF"
parse_col <- function(col) {
  if (grepl("^iSensors_all_", col)) {
    list(method = "iSensors", sk = sub("^iSensors_all_", "", col), pool = "all")
  } else if (grepl("^all_UCell_", col)) {
    list(method = "UCell", sk = sub("^all_UCell_", "", col), pool = "all")
  } else if (grepl("^all_AUCell_", col)) {
    list(method = "AUCell", sk = sub("^all_AUCell_", "", col), pool = "all")
  } else if (grepl("^iSensors_", col)) {
    list(method = "iSensors", sk = sub("^iSensors_", "", col), pool = "sel")
  } else if (grepl("^UCell_", col)) {
    list(method = "UCell", sk = sub("^UCell_", "", col), pool = "sel")
  } else if (grepl("^AUCell_", col)) {
    list(method = "AUCell", sk = sub("^AUCell_", "", col), pool = "sel")
  } else {
    list(method = NA, sk = col, pool = NA)
  }
}

parsed <- lapply(spear_all$col, parse_col)
spear_all$Method    <- vapply(parsed, function(x) x$method, character(1))
spear_all$SensorKey <- vapply(parsed, function(x) x$sk,     character(1))
spear_all$Pool      <- vapply(parsed, function(x) x$pool,   character(1))

# Selected-sensor table
spear_sel <- spear_all %>%
  filter(Pool == "sel", SensorKey %in% names(gene_sets_sel),
         !is.na(Method)) %>%
  mutate(
    Method = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
    DisplayLabel = sel_display[SensorKey],
    DisplayLabel = factor(DisplayLabel,
                          levels = rev(sel_display[names(gene_sets_sel)]))
  )

# All-sensor table (use "all" pool to avoid double-counting sel sensors)
spear_all_clean <- spear_all %>%
  filter(Pool == "all", !is.na(Method)) %>%
  mutate(
    Method = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
    DisplayLabel = display_label(panel_key_map[SensorKey])
  ) %>%
  group_by(SensorKey) %>%
  mutate(isens_rho = rho[Method == "iSensors"][1]) %>%
  ungroup() %>%
  arrange(isens_rho) %>%
  mutate(DisplayLabel = factor(DisplayLabel, levels = unique(DisplayLabel)))

# Save table
write.csv(bind_rows(
  spear_sel %>% select(Method, SensorKey, DisplayLabel, Spearman_rho = rho),
  spear_all_clean %>% select(Method, SensorKey, DisplayLabel, Spearman_rho = rho)
), file.path(output_dir, "gradient_spearman.csv"), row.names = FALSE)

# Helper: gradient dotplot
make_spear_plot <- function(df, title,
                            x_lim = c(-1.01, 1.01)) {
  ggplot(df, aes(x = rho, y = DisplayLabel, color = Method, shape = Method)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.5) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted",
               color = "grey70", linewidth = 0.3) +
    geom_point(size = 2.8, alpha = 0.9,
               position = position_dodge(width = 0.6)) +
    scale_color_manual(values = method_colors) +
    scale_shape_manual(values = method_shapes) +
    scale_x_continuous(limits = x_lim, breaks = seq(-1, 1, 0.25)) +
    labs(x = "Spearman rho (vs DR5-ARF1-up gradient)", y = NULL,
         color = NULL, shape = NULL, title = title) +
    theme_classic(base_size = 10) +
    theme(strip.background    = element_blank(),
          strip.text          = element_text(face = "bold", size = 9),
          plot.title          = element_text(face = "bold", size = 10),
          legend.position     = "bottom",
          legend.key.size     = unit(0.4, "cm"),
          axis.line           = element_line(linewidth = 0.4),
          panel.grid.major.x  = element_line(linewidth = 0.2, color = "grey88"))
}

# Panel C: selected sensors gradient dotplot
p_C <- make_spear_plot(spear_sel,
                       "Auxin gradient correlation -- selected sensors\n(Spearman rho per cell type, Shahan root tip)")

# Inset scatter: cell-type mean ARF sensor score vs DR5-ARF1-up reference
scatter_df <- ct_means %>%
  select(cell_type,
         DR5_ref     = all_of(ref_col),
         iSensors_ARF = iSensors_ARF,
         UCell_ARF   = UCell_ARF,
         AUCell_ARF  = AUCell_ARF) %>%
  pivot_longer(c(iSensors_ARF, UCell_ARF, AUCell_ARF),
               names_to = "Method", values_to = "ARF_score") %>%
  mutate(Method = sub("_ARF$", "", Method),
         Method = factor(Method, levels = c("iSensors", "UCell", "AUCell")))

rho_labs <- spear_sel %>%
  filter(SensorKey == "ARF") %>%
  mutate(label = paste0(Method, " rho=", round(rho, 2)))

p_scatter <- ggplot(scatter_df,
                    aes(x = DR5_ref, y = ARF_score, color = Method)) +
  geom_point(size = 1.6, alpha = 0.75) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  scale_color_manual(values = method_colors) +
  facet_wrap(~ Method, nrow = 1, scales = "free_y") +
  labs(x = "DR5-ARF1-up iSensor (auxin reference)",
       y = "ARF sensor score",
       title = "ARF sensor vs auxin gradient per cell type") +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 9),
        legend.position  = "none",
        axis.line        = element_line(linewidth = 0.4))

fig_C <- (p_scatter / p_C) +
  plot_layout(heights = c(1, 2.2)) +
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
  spear_all_clean,
  "Auxin gradient correlation -- all sensors (Shahan root tip)"
)
h_ext <- max(5, 0.22 * length(unique(spear_all_clean$DisplayLabel)))
ggsave(file.path(output_dir, "FigureS_C_gradient_extended.pdf"),
       plot = p_C_ext, width = 7, height = h_ext, dpi = 300, bg = "white",
       limitsize = FALSE)
ggsave(file.path(output_dir, "FigureS_C_gradient_extended.svg"),
       plot = p_C_ext, width = 7, height = h_ext, dpi = 300, bg = "white",
       limitsize = FALSE)
message("Saved FigureS_C_gradient_extended (", round(h_ext, 1), " in tall)")

message("Done.")

