# Figure 7 root diagnostics: decompose the ARF prediction into sensor contributions
# and test whether a root-specific offset improves the global rice model.
#
# Run with working directory = iSensors-supplementary/
# Input:  Manuscript-Figures/in/Figure7_rice_model_data.csv
# Output: Manuscript-Figures/out/Figure7_root_model_contributions.*
#         Manuscript-Figures/out/Figure7_root_model_coefficient_sensitivity.*
#         Manuscript-Figures/in/Figure7_root_model_diagnostics.csv
#         Manuscript-Figures/in/Figure7_root_model_tests.txt

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
})

input_file <- "Manuscript-Figures/in/Figure7_rice_model_data.csv"
output_dir <- "Manuscript-Figures/out"
table_dir <- "Manuscript-Figures/in"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  select(Organ, Tissue_Cleaned, ARF, Synthesis, Conjugation, PAT) %>%
  drop_na()

# Refit rather than relying on rounded/exported predictions.
m_global <- lm(ARF ~ 0 + Synthesis + Conjugation + PAT, data = d)
b <- coef(m_global)

d <- d %>%
  mutate(
    is_root = as.integer(Organ == "Root"),
    is_cortex = as.integer(str_detect(Tissue_Cleaned, regex("cortex", ignore_case = TRUE))),
    Synthesis_contribution = b[["Synthesis"]] * Synthesis,
    Conjugation_contribution = b[["Conjugation"]] * Conjugation,
    PAT_contribution = b[["PAT"]] * PAT,
    predicted = Synthesis_contribution + Conjugation_contribution + PAT_contribution,
    residual = ARF - predicted
  )

# A root-specific additive shift is estimable. A root-by-cortex interaction is
# not tested because root cortex is represented by only one aggregated point.
m_root_shift <- lm(
  ARF ~ 0 + Synthesis + Conjugation + PAT + is_root,
  data = d
)
m_root_cortex_shift <- lm(
  ARF ~ 0 + Synthesis + Conjugation + PAT + is_root + is_cortex,
  data = d
)

root_summary <- d %>%
  filter(Organ == "Root") %>%
  summarise(
    n = n(),
    observed_mean = mean(ARF),
    predicted_mean = mean(predicted),
    mean_residual = mean(residual),
    se_residual = sd(residual) / sqrt(n())
  )

cortex_summary <- d %>%
  filter(is_cortex == 1) %>%
  group_by(Organ, Tissue_Cleaned) %>%
  summarise(
    observed = mean(ARF),
    predicted = mean(predicted),
    residual = mean(residual),
    n = n(),
    .groups = "drop"
  )

root_cortex <- d %>%
  filter(Organ == "Root", is_cortex == 1)

diagnostic_table <- d %>%
  filter(Organ == "Root") %>%
  arrange(desc(residual)) %>%
  select(
    Organ, Tissue_Cleaned, ARF, predicted, residual,
    Synthesis_contribution, Conjugation_contribution, PAT_contribution
  )

write.csv(
  diagnostic_table,
  file.path(table_dir, "Figure7_root_model_diagnostics.csv"),
  row.names = FALSE
)

capture.output(
  {
    cat("GLOBAL MODEL\n")
    print(summary(m_global))
    cat("\nROOT-SHIFT MODEL\n")
    print(summary(m_root_shift))
    cat("\nGLOBAL VS ROOT-SHIFT MODEL\n")
    print(anova(m_global, m_root_shift))
    cat("\nROOT + CORTEX ADDITIVE-SHIFT MODEL\n")
    print(summary(m_root_cortex_shift))
    cat("\nGLOBAL VS ROOT + CORTEX ADDITIVE-SHIFT MODEL\n")
    print(anova(m_global, m_root_cortex_shift))
    cat("\nROOT SUMMARY\n")
    print(root_summary)
    cat("\nCORTEX OBSERVATIONS\n")
    print(cortex_summary)
    cat("\nROOT CORTEX\n")
    print(root_cortex)
    cat("\nNOTE\n")
    cat("Residual = observed ARF - predicted ARF. Negative values denote overprediction.\n")
    cat("No root-by-cortex interaction test was performed: root cortex has one aggregated observation.\n")
  },
  file = file.path(table_dir, "Figure7_root_model_tests.txt")
)

# Contribution decomposition for each root tissue.
root_long <- d %>%
  filter(Organ == "Root") %>%
  select(
    Tissue_Cleaned, ARF, predicted, residual,
    Synthesis = Synthesis_contribution,
    Conjugation = Conjugation_contribution,
    PAT = PAT_contribution
  ) %>%
  arrange(ARF) %>%
  mutate(Tissue_Cleaned = factor(Tissue_Cleaned, levels = Tissue_Cleaned)) %>%
  pivot_longer(
    cols = c(Synthesis, Conjugation, PAT),
    names_to = "Component",
    values_to = "Contribution"
  )

sensor_cols <- c(
  Synthesis = "#2166ac",
  Conjugation = "#4dac26",
  PAT = "#d01c8b"
)

p_contributions <- ggplot(
  root_long,
  aes(x = Tissue_Cleaned, y = Contribution, fill = Component)
) +
  geom_col(width = 0.72) +
  geom_point(
    data = distinct(root_long, Tissue_Cleaned, ARF),
    aes(x = Tissue_Cleaned, y = ARF),
    inherit.aes = FALSE,
    shape = 21, fill = "white", colour = "black", size = 2.6, stroke = 0.7
  ) +
  scale_fill_manual(values = sensor_cols) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Predicted ARF contribution",
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  )

# A separate residual panel prevents stacked contributions from obscuring the
# direction of the discrepancy. Negative = prediction exceeds observed ARF.
root_residuals <- d %>%
  filter(Organ == "Root") %>%
  mutate(
    Tissue_Cleaned = factor(
      Tissue_Cleaned,
      levels = levels(root_long$Tissue_Cleaned)
    ),
    Highlight = if_else(is_cortex == 1, "Cortex", "Other root tissue")
  )

p_residuals <- ggplot(root_residuals, aes(Tissue_Cleaned, residual)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_col(aes(fill = Highlight), width = 0.62) +
  scale_fill_manual(values = c("Cortex" = "#d73027", "Other root tissue" = "grey65")) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Observed - predicted ARF",
    fill = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank()
  )

p_combined <- p_contributions + p_residuals +
  plot_layout(widths = c(1.6, 1), guides = "collect") &
  theme(legend.position = "top")

ggsave(
  file.path(output_dir, "Figure7_root_model_contributions.pdf"),
  p_combined, width = 8.2, height = 4.8, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "Figure7_root_model_contributions.png"),
  p_combined, width = 8.2, height = 4.8, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "Figure7_root_model_contributions.svg"),
  p_combined, width = 8.2, height = 4.8, dpi = 300, bg = "white"
)

# Check whether roots materially alter the three global sensor coefficients.
coefficient_frame <- function(model, fit_name) {
  tab <- summary(model)$coefficients
  data.frame(
    Term = rownames(tab),
    Estimate = tab[, "Estimate"],
    SE = tab[, "Std. Error"],
    Fit = fit_name,
    row.names = NULL
  ) %>%
    filter(Term %in% c("Synthesis", "Conjugation", "PAT"))
}

m_no_root <- lm(
  ARF ~ 0 + Synthesis + Conjugation + PAT,
  data = filter(d, Organ != "Root")
)

coef_compare <- bind_rows(
  coefficient_frame(m_global, "All tissues"),
  coefficient_frame(m_no_root, "Roots excluded")
) %>%
  mutate(
    Term = factor(Term, levels = c("PAT", "Conjugation", "Synthesis")),
    Fit = factor(Fit, levels = c("All tissues", "Roots excluded"))
  )

write.csv(
  coef_compare,
  file.path(table_dir, "Figure7_root_coefficient_sensitivity.csv"),
  row.names = FALSE
)

p_coef <- ggplot(coef_compare, aes(Term, Estimate, colour = Fit)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbar(
    aes(ymin = Estimate - 1.96 * SE, ymax = Estimate + 1.96 * SE),
    width = 0.12,
    position = position_dodge(width = 0.35)
  ) +
  geom_point(size = 2.6, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = c("All tissues" = "black", "Roots excluded" = "#d95f02")) +
  coord_flip() +
  labs(x = NULL, y = "Scaling coefficient (95% CI)", colour = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", panel.grid = element_blank())

ggsave(
  file.path(output_dir, "Figure7_root_model_coefficient_sensitivity.pdf"),
  p_coef, width = 4.8, height = 3.3, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "Figure7_root_model_coefficient_sensitivity.png"),
  p_coef, width = 4.8, height = 3.3, dpi = 300, bg = "white"
)
ggsave(
  file.path(output_dir, "Figure7_root_model_coefficient_sensitivity.svg"),
  p_coef, width = 4.8, height = 3.3, dpi = 300, bg = "white"
)

message("Saved Figure 7 root model diagnostics.")
message(
  "Root cortex: observed = ", round(root_cortex$ARF, 4),
  ", predicted = ", round(root_cortex$predicted, 4),
  ", residual = ", round(root_cortex$residual, 4),
  " (negative means overprediction)."
)
