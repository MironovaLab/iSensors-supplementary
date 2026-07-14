# Figure 6 (no-root, DR5-ARF3-up variant): iSensor activity across the Guo
# et al. 2025 multi-organ atlas. Identical to Figure6_noroot.R but the linear
# model uses the DR5-ARF3-up reg-panel (34 genes, DEG-derived, induced under
# the ARF3-driven DR5 synthetic reporter context) as the response variable
# instead of the full ARF trans-panel:
#   DR5-ARF3-up = beta1*Synthesis + beta2*Conjugation + beta3*PAT
# All outputs are saved under distinct "DR5-ARF3-up" filenames so the
# original ARF-based figure/CSV outputs are not overwritten.
#
# Implementation note: "DR5-ARF3-up" is not a syntactically valid R name
# (the hyphens would be parsed as subtraction in formulas/`$` access), so
# the response column is renamed to a generic "Response" immediately after
# pivot_wider; RESPONSE_SENSOR (the raw string) is only used for filtering
# and display labels, where hyphens are fine.
#
# Run with working directory = iSensors-supplementary/

library(Seurat)
library(iSensors)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)

RESPONSE_SENSOR <- "DR5-ARF3-up"

output_dir <- "Manuscript-Figures/out"

# ── Load combined AverageExpression list ──────────────────────────────────────
guo_avg_exp <- readRDS("00-iSensors-objects/data/Guo/guo_avg_exp_iSensors_mean.rds")
message("Loaded ", length(guo_avg_exp), " organs: ", paste(names(guo_avg_exp), collapse = ", "))

# ── Helpers ───────────────────────────────────────────────────────────────────
polish_name <- function(x) {
  x <- gsub("^ATH-aux-cistrans-|^ATH-aux-cis-|^ATH-aux-trans-|^ATH-aux-reg-", "", x)
  x <- gsub("^AT-aux-cistrans-|^AT-aux-cis-|^AT-aux-trans-|^AT-aux-reg-",     "", x)
  gsub("PolarAuxinTransport", "PAT", x)
}

organ_rename_map <- c(
  cauline      = "Cauline 42 DAG",
  D0_silique   = "Silique 0 DPA",
  D2_silique   = "Silique 2 DPA",
  D3_silique   = "Silique 3 DPA",
  D4_silique   = "Silique 4 DPA",
  D5_silique   = "Silique 5 DPA",
  D6_root      = "Root 6 DAG",
  D11_root     = "Root 11 DAG",
  Early_flower = "Early flower",
  Middle_flower= "Middle flower",
  Late_flower  = "Late flower",
  S1_leaf      = "Rosette 14 DAG",
  S2_leaf      = "Rosette 21 DAG",
  S3_leaf      = "Rosette 28 DAG",
  S4_leaf      = "Rosette 35 DAG",
  S5_leaf      = "Rosette 42 DAG",
  S6_leaf      = "Rosette 49 DAG",
  stem         = "Stem 42 DAG"
)

organ_order <- c(
  "Rosette 14 DAG", "Rosette 21 DAG", "Rosette 28 DAG",
  "Rosette 35 DAG", "Rosette 42 DAG", "Rosette 49 DAG",
  "Early flower", "Middle flower", "Late flower",
  "Silique 0 DPA", "Silique 2 DPA", "Silique 3 DPA",
  "Silique 4 DPA", "Silique 5 DPA",
  "Cauline 42 DAG", "Stem 42 DAG"
)

standardize_tissue <- function(x) {
  x <- str_trim(str_to_lower(x))
  x <- str_remove(x, "-\\d+$")
  x <- if_else(str_detect(x, "^(phloem|phloem[- ]parenchyma.*|protophloem)$"), "phloem", x)
  x <- if_else(str_detect(x, "^(xylem|protoxylem.*|xylem/procambium)$"),        "xylem", x)
  x <- if_else(str_detect(x, "atrichoblast|trichoblast|epidermis"),              "epidermis", x)
  x <- if_else(str_detect(x, "seed coat.*"),                                     "seed coat", x)
  x <- if_else(str_detect(x, "^(pollen|sperm|tapetum|microsporocyte)$"),         "male reproductive", x)
  x <- if_else(str_detect(x, "^(integument|stigma|nectary)$"),                  "female reproductive", x)
  x <- if_else(str_detect(x, "embryo|endosperm"),                               "endosperm/embryo", x)
  x <- if_else(x %in% c("vascular", "vascular cells", "interfascicular fibers"),"vascular (unspecified)", x)
  x <- if_else(x %in% c("g2m stage", "s stage"),                                "dividing cell", x)
  x <- if_else(x %in% c("companion cell"),                                       "phloem", x)
  x <- if_else(x %in% c("guard cells"),                                          "guard cell", x)
  x <- if_else(x %in% c("bundle sheaths"),                                       "bundle sheath", x)
  x
}

# ── Build long-format data, excluding root samples ────────────────────────────
all_data <- lapply(names(guo_avg_exp), function(org_key) {
  mat     <- guo_avg_exp[[org_key]]
  org_lab <- organ_rename_map[org_key]
  if (is.na(org_lab)) org_lab <- org_key
  rownames(mat) <- polish_name(rownames(mat))
  as.data.frame(mat) %>%
    tibble::rownames_to_column("sensor") %>%
    pivot_longer(-sensor, names_to = "tissue", values_to = "value") %>%
    mutate(organ = org_lab, tissue_simple = standardize_tissue(tissue))
}) %>%
  bind_rows() %>%
  filter(!organ %in% c("Root 6 DAG", "Root 11 DAG")) %>%
  mutate(organ = factor(organ, levels = organ_order))

all_data_agg <- all_data %>%
  group_by(organ, tissue_simple, sensor) %>%
  summarise(value = max(value, na.rm = TRUE), .groups = "drop")

message("Available sensors after polish_name:")
print(sort(unique(all_data_agg$sensor)))
stopifnot(RESPONSE_SENSOR %in% all_data_agg$sensor)


# ── Figure 6A (no-root): 1×4 stacked barplots ────────────────────────────────
sensors_6a <- c(RESPONSE_SENSOR, "Synthesis", "ConjugationDeconjugation", "PAT")
sensors_6a  <- intersect(sensors_6a, unique(all_data_agg$sensor))

tissue_order_6a <- all_data_agg %>%
  filter(sensor %in% sensors_6a) %>%
  distinct(organ, tissue_simple) %>%
  count(tissue_simple, name = "n_organs") %>%
  arrange(n_organs) %>%
  pull(tissue_simple)

tissue_cols <- setNames(
  colorRampPalette(c(
    brewer.pal(12, "Paired"),
    brewer.pal(8,  "Set2"),
    brewer.pal(8,  "Dark2")
  ))(length(tissue_order_6a)),
  tissue_order_6a
)

bar_ylabs <- setNames(
  c(RESPONSE_SENSOR, "Synthesis", "Conjugation", "PAT"),
  c(RESPONSE_SENSOR, "Synthesis", "ConjugationDeconjugation", "PAT")
)

make_bar_6a <- function(sensor_name, show_x = FALSE) {
  df <- all_data_agg %>%
    filter(sensor == sensor_name) %>%
    mutate(
      organ = factor(organ, levels = organ_order),
      tissue_simple = factor(tissue_simple, levels = tissue_order_6a)
    )

  ggplot(df, aes(x = organ, y = value, fill = tissue_simple)) +
    geom_col(position = "stack", width = 0.85) +
    scale_fill_manual(values = tissue_cols, name = "Cell type", drop = FALSE) +
    scale_x_discrete(limits = organ_order[organ_order %in% levels(df$organ)]) +
    labs(x = NULL, y = bar_ylabs[[sensor_name]]) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = if (show_x) {
        element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 7)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line(linewidth = 0.3)
      } else {
        element_blank()
      },
      axis.line = element_line(linewidth = 0.4),
      axis.ticks.y = element_line(linewidth = 0.3),
      panel.grid = element_blank(),
      legend.text = element_text(size = 5),
      legend.key.size = unit(0.2, "cm"),
      legend.title = element_text(size = 6)
    )
}

plots_6a <- lapply(seq_along(sensors_6a), function(i) {
  make_bar_6a(
    sensors_6a[i],
    show_x = i == length(sensors_6a)
  )
})

p_6a <- wrap_plots(plots_6a, nrow = 4, guides = "collect") &
  theme(legend.position = "bottom")

p_6a

ggsave(file.path(output_dir, sprintf("Figure6A_noroot_%s_stacked_barplot_4sensors.pdf", RESPONSE_SENSOR)),
       plot = p_6a, width = 4, height = 10, dpi = 300, bg = "white")
ggsave(file.path(output_dir, sprintf("Figure6A_noroot_%s_stacked_barplot_4sensors.svg", RESPONSE_SENSOR)),
       plot = p_6a, width = 4, height = 10, dpi = 300, bg = "white")
message("Saved Figure 6A no-root (1x4)")

fig6a_csv_path <- sprintf("Manuscript-Figures/in/Figure6A_noroot_%s_avg_exp_per_organ_tissue.csv", RESPONSE_SENSOR)
write.csv(
  filter(all_data_agg, sensor %in% sensors_6a),
  fig6a_csv_path,
  row.names = FALSE
)
message("Exported Figure 6A no-root data table")


# ══════════════════════════════════════════════════════════════════════════════
# Figures 6C-E (no-root): Linear model - Response ~ Synthesis + Conjugation + PAT
# ══════════════════════════════════════════════════════════════════════════════

fig6a_csv <- read.csv(fig6a_csv_path, stringsAsFactors = FALSE)

df_wide <- fig6a_csv %>%
  filter(sensor %in% c(RESPONSE_SENSOR, "Synthesis", "ConjugationDeconjugation", "PAT")) %>%
  pivot_wider(id_cols     = c(organ, tissue_simple),
              names_from  = sensor,
              values_from = value) %>%
  rename(Organ = organ, Tissue_Cleaned = tissue_simple,
         Conjugation = ConjugationDeconjugation,
         Response = all_of(RESPONSE_SENSOR)) %>%
  drop_na()

m2 <- lm(Response ~ 0 + Synthesis + Conjugation + PAT, data = df_wide)

df_wide$Response_hat <- predict(m2)
df_wide$resid        <- resid(m2)
rmse <- sqrt(mean(df_wide$resid^2))
rse  <- summary(m2)$sigma

s      <- summary(m2)
r2_adj <- round(s$adj.r.squared, 2)
fstat  <- s$fstatistic
pval   <- pf(fstat[1], fstat[2], fstat[3], lower.tail = FALSE)
plab   <- if (pval < 2.2e-16) "p < 2.2e-16" else paste0("p = ", signif(pval, 2))


# ── Figure 6C (no-root) ───────────────────────────────────────────────────────
response_range <- range(df_wide$Response, na.rm = TRUE)
band_6c <- tibble(x = seq(response_range[1], response_range[2], length.out = 200)) %>%
  mutate(lo = x - rmse, hi = x + rmse)

p_6c <- ggplot(df_wide, aes(x = Response, y = Response_hat)) +
  geom_ribbon(data = band_6c, aes(x = x, ymin = lo, ymax = hi),
              inherit.aes = FALSE, fill = "grey80", alpha = 0.45) +
  geom_point(size = 2, alpha = 0.8, color = "#4d4d4d") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.8) +
  coord_equal(xlim = response_range, ylim = response_range) +
  annotate("text", x = response_range[1], y = response_range[2],
           label = paste0("adj. R² = ", r2_adj, "\n", plab),
           hjust = 0, vjust = 1, size = 3.2, color = "grey30") +
  labs(x = paste0("Observed ", RESPONSE_SENSOR, " iSensor"), y = "Predicted auxin level") +
  theme_classic(base_size = 11) +
  theme(axis.line  = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid = element_blank())
p_6c
ggsave(file.path(output_dir, sprintf("Figure6C_noroot_%s_scatter.pdf", RESPONSE_SENSOR)),
       plot = p_6c, width = 4, height = 4, dpi = 300, bg = "white")
ggsave(file.path(output_dir, sprintf("Figure6C_noroot_%s_scatter.svg", RESPONSE_SENSOR)),
       plot = p_6c, width = 4, height = 4, dpi = 300, bg = "white")
message("Saved Figure 6C no-root")


# ── Figure 6D (no-root) ───────────────────────────────────────────────────────
coef_tbl <- summary(m2)$coefficients
coef_df <- tibble(
  Term     = rownames(coef_tbl),
  Estimate = coef_tbl[, "Estimate"],
  SE       = coef_tbl[, "Std. Error"],
  p_value  = coef_tbl[, "Pr(>|t|)"]
) %>%
  mutate(
    signif  = case_when(p_value < 0.001 ~ "***",
                        p_value < 0.01  ~ "**",
                        p_value < 0.05  ~ "*",
                        TRUE            ~ "ns"),
    label_y = if_else(Estimate >= 0,
                      Estimate + SE + 0.05 * max(abs(Estimate + SE)),
                      Estimate - SE - 0.05 * max(abs(Estimate - SE))),
    Term    = factor(Term, levels = rev(c("Synthesis", "Conjugation", "PAT")))
  )

p_6d <- ggplot(coef_df, aes(x = Term, y = Estimate)) +
  geom_col(fill = "grey65", width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                width = 0.2, linewidth = 0.5) +
  geom_text(aes(y = label_y, label = signif), size = 4, vjust = 0.5) +
  labs(x = NULL, y = "Scaling coefficient") +
  theme_classic(base_size = 11) +
  theme(axis.line  = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid = element_blank())
p_6d
ggsave(file.path(output_dir, sprintf("Figure6D_noroot_%s_model_coefficients.pdf", RESPONSE_SENSOR)),
       plot = p_6d, width = 3.5, height = 3, dpi = 300, bg = "white")
ggsave(file.path(output_dir, sprintf("Figure6D_noroot_%s_model_coefficients.svg", RESPONSE_SENSOR)),
       plot = p_6d, width = 3.5, height = 3, dpi = 300, bg = "white")
message("Saved Figure 6D no-root")


# ── Figure 6E (no-root) ───────────────────────────────────────────────────────
make_forest <- function(summary_df, x_col, title_lab) {
  df <- summary_df %>% arrange(mean_resid) %>%
    mutate(!!x_col := factor(.data[[x_col]], levels = .data[[x_col]]))

  ggplot(df, aes(x = .data[[x_col]], y = mean_resid)) +
    geom_rect(ymin = -rse,       ymax =  rse,
              xmin = -Inf, xmax = Inf, fill = "grey90", alpha = 0.6) +
    geom_rect(ymin = -0.3 * rse, ymax =  0.3 * rse,
              xmin = -Inf, xmax = Inf, fill = "grey75", alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5,
               color = "grey40") +
    geom_errorbar(aes(ymin = mean_resid - se_resid,
                      ymax = mean_resid + se_resid),
                  width = 0.25, linewidth = 0.5) +
    geom_point(size = 2.5) +
    coord_flip() +
    labs(x = NULL, y = "Mean residual\n(observed - predicted)", title = title_lab) +
    theme_classic(base_size = 10) +
    theme(plot.title  = element_text(face = "bold", size = 10),
          axis.line   = element_line(linewidth = 0.4),
          axis.ticks  = element_line(linewidth = 0.3),
          axis.text.y = element_text(size = 8),
          panel.grid  = element_blank())
}

organ_summary <- df_wide %>%
  group_by(Organ) %>%
  summarise(mean_resid = mean(resid),
            se_resid   = sd(resid) / sqrt(n()),
            n          = n(), .groups = "drop")

tissue_summary <- df_wide %>%
  group_by(Tissue_Cleaned) %>%
  summarise(mean_resid = mean(resid),
            se_resid   = sd(resid) / sqrt(n()),
            n          = n(), .groups = "drop") %>%
  filter(n >= 3)

p_6e <- make_forest(organ_summary,  "Organ",         "Organ level") +
        make_forest(tissue_summary, "Tissue_Cleaned","Tissue level")

p_6e
ggsave(file.path(output_dir, sprintf("Figure6E_noroot_%s_residuals_forest.pdf", RESPONSE_SENSOR)),
       plot = p_6e, width = 4.5, height = 6, dpi = 300, bg = "white")
ggsave(file.path(output_dir, sprintf("Figure6E_noroot_%s_residuals_forest.svg", RESPONSE_SENSOR)),
       plot = p_6e, width = 4.5, height = 6, dpi = 300, bg = "white")
message("Saved Figure 6E no-root")


# ── Exploratory: forest plot of model residuals per reproductive cell type ─────
repro_wide_raw <- all_data %>%
  filter(tissue_simple %in% c("male reproductive", "female reproductive"),
         sensor %in% c(RESPONSE_SENSOR, "Synthesis", "ConjugationDeconjugation", "PAT")) %>%
  group_by(organ, tissue, tissue_simple, sensor) %>%
  summarise(value = max(value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = sensor, values_from = value) %>%
  rename(Conjugation = ConjugationDeconjugation,
         Response = all_of(RESPONSE_SENSOR)) %>%
  drop_na(Response, Synthesis, Conjugation, PAT) %>%
  mutate(cell_type = str_to_title(
    str_replace_all(str_remove(tissue, "-\\d+$"), "[_-]", " ")
  ))

repro_wide_raw$Response_hat <- predict(m2, newdata = repro_wide_raw)
repro_wide_raw$resid        <- repro_wide_raw$Response - repro_wide_raw$Response_hat

repro_ct_summary <- repro_wide_raw %>%
  group_by(cell_type, tissue_simple) %>%
  summarise(mean_resid = mean(resid),
            se_resid   = sd(resid) / sqrt(n()),
            n          = n(), .groups = "drop") %>%
  arrange(mean_resid) %>%
  mutate(cell_type = factor(cell_type, levels = cell_type))

p_repro_forest <- ggplot(repro_ct_summary,
                          aes(x = cell_type, y = mean_resid)) +
  geom_rect(ymin = -rse,       ymax =  rse,
            xmin = -Inf, xmax = Inf,
            fill = "grey90", alpha = 0.6, inherit.aes = FALSE) +
  geom_rect(ymin = -0.3 * rse, ymax =  0.3 * rse,
            xmin = -Inf, xmax = Inf,
            fill = "grey75", alpha = 0.8, inherit.aes = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5,
             colour = "grey40") +
  geom_errorbar(aes(ymin = mean_resid - se_resid,
                    ymax = mean_resid + se_resid),
                width = 0.25, linewidth = 0.5) +
  geom_point(size = 3) +
  facet_wrap(~ tissue_simple, scales = "free_y", ncol = 2) +
  coord_flip() +
  labs(x = NULL,
       y = paste0("Mean residual (observed ", RESPONSE_SENSOR, " - predicted)")) +
  theme_classic(base_size = 10) +
  theme(axis.line   = element_line(linewidth = 0.4),
        axis.ticks  = element_line(linewidth = 0.3),
        axis.text.y = element_text(size = 9),
        panel.grid  = element_blank(),
        strip.text  = element_text(face = "bold", size = 10))
p_repro_forest

ggsave(file.path(output_dir, sprintf("Figure6_noroot_%s_reproductive_residuals_forest.pdf", RESPONSE_SENSOR)),
       plot = p_repro_forest, width = 4.5, height = 1.5, dpi = 300, bg = "white")
ggsave(file.path(output_dir, sprintf("Figure6_noroot_%s_reproductive_residuals_forest.svg", RESPONSE_SENSOR)),
       plot = p_repro_forest, width = 4.5, height = 1.5, dpi = 300, bg = "white")
message("Saved reproductive tissue residuals forest plot (no root)")


# ═══════════════════════════════════════════════════════════════════════════════
# Figure 6 no-root combined: publication-ready layout
#
# Left  30%: A - 4 stacked barplots (nrow = 4), legend at bottom (4 cols)
# Right 70%: Row1: B (half) | C (quarter) | D (quarter)
#            Row2: E (4/5) | F (1/5)
#               E = organ forest | tissue forest
#               F = female repro (top) / male repro (bottom)
# ═══════════════════════════════════════════════════════════════════════════════

# ── Panel A: 4 stacked barplots stacked vertically, legend at bottom ──────────
plots_6a_comb <- lapply(seq_along(sensors_6a), function(i)
  make_bar_6a(sensors_6a[i], show_x = (i == length(sensors_6a)))
)
p_A <- wrap_plots(plots_6a_comb, nrow = 4, guides = "collect") &
  theme(legend.position  = "bottom",
        legend.direction  = "horizontal",
        legend.key.size   = unit(0.22, "cm"),
        legend.text       = element_text(size = 5),
        legend.title      = element_text(size = 6)) &
  guides(fill = guide_legend(ncol = 4, byrow = TRUE))

# ── Panel B: biological model -> statistical formula ──────────────────────────
p_B <- ggplot() +
  annotate("text", x = 0.5, y = 0.82,
           label = "Auxin = Synthesis + Conjugation + Transport",
           hjust = 0.5, vjust = 0.5, size = 4.2, colour = "black") +
  annotate("segment", x = 0.5, xend = 0.5, y = 0.64, yend = 0.52,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
           colour = "black", linewidth = 0.8) +
  annotate("text", x = 0.5, y = 0.38,
           label = as.expression(bquote(.(RESPONSE_SENSOR) * " = " * beta[1] * "·Synthesis + " *
                              beta[2] * "·Conjugation + " * beta[3] * "·PAT")),
           hjust = 0.5, vjust = 0.5, size = 4.0, colour = "#2166ac") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()

# ── Split reproductive forest: female (top) / male (bottom) ──────────────────
make_repro_forest_sex <- function(sex_label) {
  df <- repro_ct_summary %>%
    filter(tissue_simple == sex_label) %>%
    arrange(mean_resid) %>%
    mutate(cell_type = factor(cell_type, levels = cell_type))
  ggplot(df, aes(x = cell_type, y = mean_resid)) +
    geom_rect(ymin = -rse,       ymax =  rse,
              xmin = -Inf, xmax = Inf,
              fill = "grey90", alpha = 0.6, inherit.aes = FALSE) +
    geom_rect(ymin = -0.3 * rse, ymax =  0.3 * rse,
              xmin = -Inf, xmax = Inf,
              fill = "grey75", alpha = 0.8, inherit.aes = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5,
               colour = "grey40") +
    geom_errorbar(aes(ymin = mean_resid - se_resid,
                      ymax = mean_resid + se_resid),
                  width = 0.25, linewidth = 0.5) +
    geom_point(size = 2.5) +
    coord_flip() +
    labs(x = NULL,
         y = "Mean residual (obs - pred)",
         title = str_to_title(sex_label)) +
    theme_classic(base_size = 10) +
    theme(plot.title  = element_text(face = "bold", size = 9),
          axis.line   = element_line(linewidth = 0.4),
          axis.ticks  = element_line(linewidth = 0.3),
          axis.text.y = element_text(size = 8),
          panel.grid  = element_blank())
}

p_female <- make_repro_forest_sex("female reproductive")
p_male   <- make_repro_forest_sex("male reproductive")
p_F_comb <- p_female / p_male

# ── Panel E: organ + tissue forest (re-created as separate objects) ───────────
p_6e_organ  <- make_forest(organ_summary,  "Organ",          "Organ level")
p_6e_tissue <- make_forest(tissue_summary, "Tissue_Cleaned", "Tissue level")
p_E_comb    <- p_6e_organ + p_6e_tissue + plot_layout(ncol = 2)

# ── Panel D (colored): coefficient barplot ────────────────────────────────────
p_6d_col <- ggplot(coef_df, aes(x = Term, y = Estimate, fill = Term)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                width = 0.2, linewidth = 0.5) +
  geom_text(aes(y = label_y, label = signif), size = 4, vjust = 0.5) +
  scale_fill_manual(
    values = c(Synthesis = "#2166ac", Conjugation = "#4dac26", PAT = "#d01c8b"),
    guide  = "none"
  ) +
  labs(x = NULL, y = "Scaling coefficient") +
  theme_classic(base_size = 11) +
  theme(axis.line  = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.3),
        panel.grid = element_blank())

# ── Row 1: B (half width) | C (quarter) | D (quarter) ────────────────────────
row_BCD <- (p_B | p_6c | p_6d_col) + plot_layout(widths = c(2, 1, 1))

# ── Row 2: E (4/5 width) | F (1/5 width, stacked female/male) ────────────────
row_EF <- (p_E_comb | p_F_comb) + plot_layout(widths = c(4, 1))

# ── Right column: BCD row (narrow) / EF row (tall) ───────────────────────────
right_col <- (row_BCD / row_EF) + plot_layout(heights = c(1, 4))

# ── Full figure (widths 3:7) with panel labels ────────────────────────────────
# patchwork depth-first leaf order:
#   A(4 barplots), B, C, D, E(organ+tissue=2), F(female+male=2)  -> 11 leaves
fig6_combined <- ((p_A | right_col) +
  plot_layout(widths = c(3, 7)) +
  plot_annotation(
    tag_levels = list(c("A", "", "", "", "B", "C", "D", "E", "", "F", "")),
    theme = theme(plot.tag = element_text(face = "bold", size = 13))
  )) & theme(plot.margin = margin(4, 8, 4, 8, "pt"))

combined_pdf <- file.path(output_dir, sprintf("Figure6_noroot_%s_combined.pdf", RESPONSE_SENSOR))
combined_svg <- file.path(output_dir, sprintf("Figure6_noroot_%s_combined.svg", RESPONSE_SENSOR))
combined_png <- file.path(output_dir, sprintf("Figure6_noroot_%s_combined.png", RESPONSE_SENSOR))
ggsave(combined_pdf, plot = fig6_combined, width = 20, height = 22, dpi = 300, bg = "white")
ggsave(combined_svg, plot = fig6_combined, width = 20, height = 22, dpi = 300, bg = "white")
ggsave(combined_png, plot = fig6_combined, width = 20, height = 22, dpi = 150, bg = "white")
message("Saved Figure 6 no-root combined: ", combined_pdf, " / .svg / .png")
