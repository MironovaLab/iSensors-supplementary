# Figure8_BR_concentration_duration.R
# Manuscript figure: BR iSensor sensitivity to exogenous BR concentration and
# treatment duration, validated on GSE212230 (Nolan lab BL time-course /
# BRZ dataset). 1 column, 2 rows:
#   Row 1 - BR concentration response (BRZ=low, Control=intermediate, BL 8h=high)
#   Row 2 - BR duration response (BL 0.5-8h time course, single 100nM dose)
#
# Statistics: per-CELL Spearman rank correlation between iSensor score and the
# ordinal exposure variable (BR level or hours), following the Auxin manuscript
# methodology. Significance bar is derived EMPIRICALLY as the max |rho| among
# 3 negative-control panels (random1, random2, majortrend; generic iSensors
# panels already scored on this object) rather than using a fixed threshold -
# the fixed Auxin |rho| > 0.45 bar was too strict for BR single-cell data (see
# BRiSensors/20-...-spearman.R); the empirical bar is derived per-axis in
# BRiSensors/21-BR-concentration-duration-empirical-threshold.R.
#
# Panel naming: "4hr" -> "late", "2hr" -> "early" (e.g. reg-4hr-down -> reg-late-down)
#
# Run with working directory = iSensors-supplementary/
# Input:  ../BRiSensors/out/GSE212230_iSensors_obj.rds
# Output: Manuscript-Figures/out/Figure8_BR_concentration_duration.pdf / .png / .svg
#         Manuscript-Figures/out/Figure8_BR_concentration_stats.csv
#         Manuscript-Figures/out/Figure8_BR_duration_stats.csv

suppressPackageStartupMessages({
  library(Seurat); library(iSensors)
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(patchwork)
})

output_dir <- "Manuscript-Figures/out"
dir.create(output_dir, showWarnings = FALSE)

BR_SENSORS <- c(
  "ATH-br-reg-4hr-down", "ATH-br-reg-2hr-down", "ATH-br-reg-4hr-up",
  "ATH-br-trans-NegativeSignaling", "ATH-br-reg-2hr-up", "ATH-br-trans-TF-induced",
  "ATH-br-trans-Homeostasis", "ATH-br-trans-TF-repressed", "ATH-br-trans-TF",
  "ATH-br-trans-PositiveSignaling", "ATH-br-trans-Biosynthesis"
)
NEG_CONTROLS <- c("random1", "random2", "majortrend")
SENSORS <- c(BR_SENSORS, NEG_CONTROLS)

# "4hr" -> "late", "2hr" -> "early"
SENSOR_LABELS <- c(
  "ATH-br-reg-4hr-down"             = "reg-late-down",
  "ATH-br-reg-2hr-down"             = "reg-early-down",
  "ATH-br-reg-4hr-up"               = "reg-late-up",
  "ATH-br-trans-NegativeSignaling"  = "trans-NegSignaling",
  "ATH-br-reg-2hr-up"               = "reg-early-up",
  "ATH-br-trans-TF-induced"         = "trans-TF-induced",
  "ATH-br-trans-Homeostasis"        = "trans-Homeostasis",
  "ATH-br-trans-TF-repressed"       = "trans-TF-repressed",
  "ATH-br-trans-TF"                 = "trans-TF",
  "ATH-br-trans-PositiveSignaling"  = "trans-PosSignaling",
  "ATH-br-trans-Biosynthesis"       = "trans-Biosynthesis",
  "random1"                         = "neg-random1",
  "random2"                         = "neg-random2",
  "majortrend"                      = "neg-majortrend"
)
LABEL_ORDER <- SENSOR_LABELS[SENSORS]  # BR sensors first (best -> worst), neg controls last

cat("Loading GSE212230 iSensors object...\n")
obj <- readRDS("00-BR-objects/out/GSE212230_iSensors_obj.rds")
DefaultAssay(obj) <- "iSensors_mean"

meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, sample, cell_type, condition)

sens <- FetchData(obj, vars = SENSORS) %>% tibble::rownames_to_column("cell")

df <- meta %>% inner_join(sens, by = "cell")
rm(obj); invisible(gc())

# ══════════════════════════════════════════════════════════════════════════
# Helper: per-cell Spearman correlation, Fisher-z 95% CI, empirical threshold
# ══════════════════════════════════════════════════════════════════════════
spearman_per_sensor <- function(data, predictor_col, sensor_list) {
  purrr::map_dfr(sensor_list, function(s) {
    x <- data[[predictor_col]]
    y <- data[[s]]
    ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
    n  <- length(x)
    rho <- unname(ct$estimate)
    z  <- atanh(pmin(pmax(rho, -0.9999), 0.9999))
    se <- 1 / sqrt(n - 3)
    tibble(sensor = s, rho = rho, p.value = ct$p.value, n = n,
           lo = tanh(z - 1.96 * se), hi = tanh(z + 1.96 * se))
  })
}

classify_empirical <- function(res) {
  thr <- max(abs(res$rho[res$sensor %in% NEG_CONTROLS]))
  res <- res %>%
    mutate(
      p_adj = p.adjust(p.value, method = "BH"),
      class = case_when(
        sensor %in% NEG_CONTROLS ~ "Negative control",
        abs(rho) > thr & p.value < 0.05 & rho > 0 ~ "Up",
        abs(rho) > thr & p.value < 0.05 & rho < 0 ~ "Down",
        TRUE ~ "n.s."
      ),
      sensor_label = factor(SENSOR_LABELS[sensor], levels = rev(LABEL_ORDER))
    )
  attr(res, "threshold") <- thr
  res
}

make_forest <- function(res, thr, xlab, title) {
  rng <- range(c(res$lo, res$hi), na.rm = TRUE)
  pad <- diff(rng) * 0.08
  xlim <- c(rng[1] - pad, rng[2] + pad)

  ggplot(res, aes(x = sensor_label, y = rho)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = c(-thr, thr), linetype = "dotted", color = "grey40", linewidth = 0.5) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.4, color = "grey50") +
    geom_point(aes(color = class), size = 2) +
    scale_color_manual(values = c(Up = "#d73027", Down = "#4575b4", "n.s." = "grey75",
                                   "Negative control" = "black"), guide = "none") +
    coord_flip(ylim = xlim) +
    labs(x = NULL, y = xlab, title = title) +
    theme_classic(base_size = 11) +
    theme(axis.text.y = element_text(size = 7.5),
          axis.text.x = element_text(size = 7.5),
          axis.title.x = element_text(size = 8.5, lineheight = 0.9),
          plot.title = element_text(size = 9, face = "bold", lineheight = 1.0),
          plot.margin = margin(4, 6, 4, 2))
}

# ══════════════════════════════════════════════════════════════════════════
# Row 1 - BR concentration: BRZ (low) < Control (intermediate) < BL 8h (high)
# ══════════════════════════════════════════════════════════════════════════
conc_cells <- df %>%
  filter(condition %in% c("BRZ", "Control", "BL_8h")) %>%
  mutate(level = case_when(condition == "BRZ" ~ 1, condition == "Control" ~ 2, condition == "BL_8h" ~ 3))

res_conc <- classify_empirical(spearman_per_sensor(conc_cells, "level", SENSORS))
thr_conc <- attr(res_conc, "threshold")

p_conc <- make_forest(res_conc, thr_conc,
                       "Spearman rho\n(score vs BR level)",
                       "BR concentration response\n(low → intermediate → high)")

# ══════════════════════════════════════════════════════════════════════════
# Row 2 - BR duration: BL time course (0.5-8h), single dose
# ══════════════════════════════════════════════════════════════════════════
dur_cells <- df %>%
  filter(condition %in% c("BL_0.5h", "BL_1h", "BL_2h", "BL_4h", "BL_8h")) %>%
  mutate(hours = case_when(
    condition == "BL_0.5h" ~ 0.5, condition == "BL_1h" ~ 1, condition == "BL_2h" ~ 2,
    condition == "BL_4h" ~ 4, condition == "BL_8h" ~ 8))

res_dur <- classify_empirical(spearman_per_sensor(dur_cells, "hours", SENSORS))
thr_dur <- attr(res_dur, "threshold")

p_duration <- make_forest(res_dur, thr_dur,
                           "Spearman rho\n(score vs BL exposure duration)",
                           "BR duration response\n(0.5, 1, 2, 4, 8 hrs)")

# ══════════════════════════════════════════════════════════════════════════
# Combine (1 column, 2 rows) and save
# ══════════════════════════════════════════════════════════════════════════
p_final <- p_conc / p_duration + plot_layout(heights = c(1, 1))

ggsave(file.path(output_dir, "Figure8_BR_concentration_duration.pdf"), plot = p_final, width = 3.5, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure8_BR_concentration_duration.png"), plot = p_final, width = 3.5, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure8_BR_concentration_duration.svg"), plot = p_final, width = 3.5, height = 11, dpi = 300, bg = "white")
cat("Saved Figure8_BR_concentration_duration.pdf/.png/.svg\n")

write.csv(res_conc %>% select(sensor, sensor_label, rho, lo, hi, n, p.value, p_adj, class),
          file.path(output_dir, "Figure8_BR_concentration_stats.csv"), row.names = FALSE)
write.csv(res_dur %>% select(sensor, sensor_label, rho, lo, hi, n, p.value, p_adj, class),
          file.path(output_dir, "Figure8_BR_duration_stats.csv"), row.names = FALSE)
cat("Saved stats CSVs.\n")

cat("\nEmpirical threshold (concentration axis):", round(thr_conc, 4), "\n")
cat("Empirical threshold (duration axis):     ", round(thr_dur, 4), "\n")
