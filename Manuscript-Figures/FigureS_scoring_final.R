# FigureS_scoring_final.R
# Assembles the final supplementary figure from pre-computed results.
#
# Panel A  -- ARF sensor score distributions by treatment (violin/box)
#             Exogenous auxin dataset, Martin-Arevalillo 2025
#             Fix: sqrt y-axis so UCell/AUCell signal is visible above the
#             large zero mass.
#
# Panel B  -- AUROC across all ATH auxin gene sets, Dataset 1
#             (formerly FigureS_B_extended_D1)
#
# Panel C  -- Spearman rank correlation with auxin cell-type gradient
#             across all ATH auxin gene sets (formerly FigureS_C_gradient_extended)
#
# Inputs (all in Manuscript-Figures/out/scoring_comparison/):
#   auxin_meta_scores.rds    -- per-cell scores, Dataset 1
#   auroc_all_panels.csv     -- AUROC for all panels, both datasets
#   gradient_spearman.csv    -- Spearman rho for all panels (groundtruth object)
#
# Output:
#   FigureS_scoring_final.pdf/.svg
#
# Run from: iSensors-supplementary/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
})

output_dir <- "Manuscript-Figures/out/scoring_comparison"

method_colors <- c(iSensors = "#2166ac", UCell = "#4dac26", AUCell = "#f4a582")
method_shapes <- c(iSensors = 16, UCell = 15, AUCell = 18)

# ==============================================================================
# Panel A -- violin: ARF sensor, exogenous auxin treatment
# ==============================================================================
message("Building Panel A ...")
auxin_meta <- readRDS(file.path(output_dir, "auxin_meta_scores.rds"))

violin_df <- auxin_meta %>%
  select(Treatment = orig.ident2,
         iSensors  = iSensors_ARF,
         UCell     = UCell_ARF,
         AUCell    = AUCell_ARF) %>%
  pivot_longer(-Treatment, names_to = "Method", values_to = "Score") %>%
  mutate(
    Method    = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
    Treatment = factor(Treatment, levels = c("Control", "Auxin"))
  ) %>%
  filter(!is.na(Score), !is.na(Treatment))

# AUROC labels -- read from the CSV so we don't re-compute
auroc_csv <- read.csv(file.path(output_dir, "auroc_selected_comparison.csv"),
                      stringsAsFactors = FALSE)
arf_auroc <- auroc_csv %>%
  filter(SensorKey == "ARF",
         grepl("Martin", Dataset)) %>%
  mutate(Method = factor(Method, levels = c("iSensors", "UCell", "AUCell")),
         label  = paste0("AUROC = ", round(AUROC, 2)))

p_A <- ggplot(violin_df,
              aes(x = Treatment, y = Score, fill = Treatment)) +
  geom_violin(trim = TRUE, scale = "width",
              linewidth = 0.3, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white",
               linewidth = 0.35, coef = 0) +
  geom_text(data = arf_auroc,
            aes(x = 1.5, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 2.8, colour = "grey25") +
  scale_fill_manual(values = c(Control = "#b2b2b2", Auxin = "#e34a33"),
                    guide = "none") +
  scale_y_sqrt(expand = expansion(c(0, 0.08))) +
  facet_wrap(~ Method, nrow = 1, scales = "free_y") +
  labs(x = NULL,
       y = "ARF sensor score  (sqrt scale)",
       title = "ARF gene set scoring by method -- exogenous auxin treatment",
       subtitle = "Martin-Arevalillo et al. 2025 (GSE241573); Control vs. +IAA") +
  theme_classic(base_size = 10) +
  theme(strip.background  = element_blank(),
        strip.text         = element_text(face = "bold", size = 9),
        plot.title         = element_text(face = "bold", size = 10),
        plot.subtitle      = element_text(size = 8, colour = "grey40"),
        axis.line          = element_line(linewidth = 0.4),
        axis.ticks         = element_line(linewidth = 0.3),
        axis.text.x        = element_text(angle = 30, hjust = 1))

# ==============================================================================
# Panel B -- AUROC dotplot, ALL sensors, Dataset 1
# ==============================================================================
message("Building Panel B ...")
auroc_all <- read.csv(file.path(output_dir, "auroc_all_panels.csv"),
                      stringsAsFactors = FALSE)

auroc_d1 <- auroc_all %>%
  filter(grepl("Martin", Dataset), !is.na(AUROC)) %>%
  mutate(Method = factor(Method, levels = c("iSensors", "UCell", "AUCell"))) %>%
  group_by(SensorKey) %>%
  mutate(isens_auroc = AUROC[Method == "iSensors"][1]) %>%
  ungroup() %>%
  arrange(isens_auroc) %>%
  mutate(DisplayLabel = factor(DisplayLabel, levels = unique(DisplayLabel)))

p_B <- ggplot(auroc_d1,
              aes(x = AUROC, y = DisplayLabel,
                  color = Method, shape = Method)) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  geom_point(size = 2.2, alpha = 0.9,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_x_continuous(limits = c(0.2, 1.01),
                     breaks = seq(0.2, 1.0, 0.2)) +
  labs(x = "AUROC",
       y = NULL,
       color = NULL, shape = NULL,
       title = "Gene set scoring performance across all ATH auxin panels",
       subtitle = "AUROC: +IAA vs. control (Martin-Arevalillo et al. 2025)") +
  theme_classic(base_size = 9) +
  theme(plot.title          = element_text(face = "bold", size = 9),
        plot.subtitle       = element_text(size = 7.5, colour = "grey40"),
        legend.position     = "bottom",
        legend.key.size     = unit(0.35, "cm"),
        axis.line           = element_line(linewidth = 0.4),
        axis.text.y         = element_text(size = 7),
        panel.grid.major.x  = element_line(linewidth = 0.2, color = "grey88"))

# ==============================================================================
# Panel C -- Spearman gradient dotplot, ALL sensors
# ==============================================================================
message("Building Panel C ...")
spear_csv <- read.csv(file.path(output_dir, "gradient_spearman.csv"),
                      stringsAsFactors = FALSE)

# Rows from spear_all_clean have SensorKey (not ShortKey); use those for extended
spear_ext <- spear_csv %>%
  filter(!is.na(SensorKey)) %>%
  mutate(Method = factor(Method, levels = c("iSensors", "UCell", "AUCell"))) %>%
  group_by(DisplayLabel) %>%
  mutate(isens_rho = {
    r <- Spearman_rho[Method == "iSensors"]
    if (length(r) > 0) r[1] else NA_real_
  }) %>%
  ungroup() %>%
  arrange(isens_rho) %>%
  mutate(DisplayLabel = factor(DisplayLabel, levels = unique(DisplayLabel)))

p_C <- ggplot(spear_ext,
              aes(x = Spearman_rho, y = DisplayLabel,
                  color = Method, shape = Method)) +
  geom_vline(xintercept = 0.5, linetype = "dashed",
             color = "grey30", linewidth = 0.7) +
  geom_point(size = 2.2, alpha = 0.9,
             position = position_dodge(width = 0.6)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_x_continuous(limits = c(-0.85, 1.05),
                     breaks = seq(-0.75, 1.0, 0.25)) +
  labs(x = "Spearman rho (sensor rank vs. auxin gradient rank)",
       y = NULL,
       color = NULL, shape = NULL,
       title = "Auxin gradient tracking across all ATH auxin panels",
       subtitle = "Spearman rho vs. literature-based cell-type auxin rank (Shahan et al. 2022);\ndashed line: rho = 0.5 significance threshold") +
  theme_classic(base_size = 9) +
  theme(plot.title          = element_text(face = "bold", size = 9),
        plot.subtitle       = element_text(size = 7.5, colour = "grey40"),
        legend.position     = "bottom",
        legend.key.size     = unit(0.35, "cm"),
        axis.line           = element_line(linewidth = 0.4),
        axis.text.y         = element_text(size = 7),
        panel.grid.major.x  = element_line(linewidth = 0.2, color = "grey88"))

# ==============================================================================
# Assemble
# ==============================================================================
message("Assembling final figure ...")

n_B <- length(unique(auroc_d1$DisplayLabel))
n_C <- length(unique(spear_ext$DisplayLabel))
row_h <- 0.19   # inches per sensor row

h_BC <- max(n_B, n_C) * row_h
h_A  <- 3.5

fig_final <- (p_A) / (p_B | p_C) +
  plot_layout(heights = c(h_A, h_BC)) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 13))
  ) & theme(plot.margin = margin(4, 8, 4, 8, "pt"))

total_h <- h_A + h_BC + 0.5   # small buffer

message(sprintf("Figure size: 14 x %.1f inches", total_h))

ggsave(file.path(output_dir, "FigureS_scoring_final.pdf"),
       plot = fig_final, width = 14, height = total_h,
       dpi = 300, bg = "white", limitsize = FALSE)
ggsave(file.path(output_dir, "FigureS_scoring_final.svg"),
       plot = fig_final, width = 14, height = total_h,
       dpi = 300, bg = "white", limitsize = FALSE)

message("Saved FigureS_scoring_final.pdf/.svg")
message("Done.")
