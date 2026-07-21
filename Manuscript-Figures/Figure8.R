# Figure8.R
# Manuscript Figure 8 - BR iSensors pipeline and validation.
#
# Row 1: A | B | C
#   Panel A   Pipeline schematic (illustrative pictograms), VERTICAL
#             orientation (3 cards stacked top-to-bottom, no captions):
#               1) BR iSensors: 7 trans + 4 reg panels
#               2) Validation: concentration, duration, cell-type
#               3) Evaluate BR and Auxin iSensors co-expression
#   Panel B   Per-cell Spearman rank correlation forest plots
#             (concentration | duration), empirical negative-control
#             threshold, GSE212230.
#   Panel C   BR response (reg-late-up) vs Auxin Response (ARF),
#             longitudinal root layout pair (Shahan atlas).
#             = original "Panel A" from
#             BRiSensors/14-BR-ARF-spatial-anticorrelation.R /
#             out/BR_ARF_anticorrelation.pdf
#
# Row 2: D | E | F
#   Panel D   Quadrant scatter (BR vs ARF, threshold lines, Spearman +
#             Fisher stats). = original "Panel C".
#   Panel E   Tissue composition (mutual exclusion vs co-localization,
#             stacked bar). = original "Panel D".
#   Panel F   BR response (reg-late-up) vs Auxin Response (ARF),
#             epidermis layout pair. = original "Panel B".
#
# Run with working directory = iSensors-supplementary/
# Input:  ../BRiSensors/out/GSE212230_iSensors_obj.rds
#         ../BRiSensors/out/shahan_BR_Aux_joined.csv
#         ../BRiSensors/out/shahan_root_BR_scores.csv
#         ../BRiSensors/out/shahan_auxin_scores.csv
#         ../BRiSensors/out/shahan_epidermis_BR_scores.csv
#         ../RealisticLayouts/out/new_ggPlantmap_epidermis.csv
# Output: Manuscript-Figures/out/Figure8.pdf / .png / .svg
#         Manuscript-Figures/out/Figure8_concentration_stats.csv
#         Manuscript-Figures/out/Figure8_duration_stats.csv

suppressPackageStartupMessages({
  library(Seurat); library(iSensors)
  library(dplyr); library(tidyr); library(tibble); library(purrr)
  library(ggplot2); library(patchwork); library(scales)
  library(ggRootCellAtlas); library(ggPlantmap)
})

output_dir <- "Manuscript-Figures/out"
dir.create(output_dir, showWarnings = FALSE)
BR_DIR <- "../BRiSensors"

COL_REG   <- "#b15928"  # matches reg color convention from the Auxin reference figure
COL_TRANS <- "#ff7f00"  # matches trans color convention from the Auxin reference figure
COL_BR    <- "#ff7f00"
COL_AUX   <- "#33a02c"

circle_df <- function(cx, cy, r, n = 80) {
  th <- seq(0, 2 * pi, length.out = n)
  tibble(x = cx + r * cos(th), y = cy + r * sin(th))
}
hex_df <- function(cx, cy, r) {
  th <- seq(0, 2 * pi, length.out = 7)
  tibble(x = cx + r * cos(th + pi / 6), y = cy + r * sin(th + pi / 6))
}

# ══════════════════════════════════════════════════════════════════════════
# Panel A - pipeline schematic, VERTICAL orientation (3 cards stacked)
# ══════════════════════════════════════════════════════════════════════════
cx <- 2.0
sh_top <- 10; sh_mid <- 5; sh_bot <- 0   # vertical shifts (bottom card unshifted)

card_box <- function(shift) tibble(xmin = 0.1, xmax = 3.9, ymin = 0.1 + shift, ymax = 3.9 + shift)

p_schematic <- ggplot() +
  geom_rect(data = bind_rows(card_box(sh_top), card_box(sh_mid), card_box(sh_bot)),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "grey97", color = "grey70", linewidth = 0.4) +
  # vertical arrows, pointing down
  geom_segment(aes(x = cx, xend = cx, y = 0.1 + sh_top - 0.05, yend = 3.9 + sh_mid + 0.1),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"), color = "grey30", linewidth = 0.6) +
  geom_segment(aes(x = cx, xend = cx, y = 0.1 + sh_mid - 0.05, yend = 3.9 + sh_bot + 0.1),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"), color = "grey30", linewidth = 0.6) +

  # ── Card 1 (top): BR iSensors panel counts ────────────────────────────
  annotate("text", x = cx, y = 3.55 + sh_top, label = "BR iSensors", fontface = "bold", size = 3.4) +
  geom_rect(aes(xmin = cx - 0.75, xmax = cx - 0.35, ymin = 1.1 + sh_top, ymax = 1.1 + sh_top + 4 * 0.26), fill = COL_REG) +
  geom_rect(aes(xmin = cx + 0.05, xmax = cx + 0.45, ymin = 1.1 + sh_top, ymax = 1.1 + sh_top + 7 * 0.26), fill = COL_TRANS) +
  annotate("text", x = cx - 0.55, y = 1.1 + sh_top + 4 * 0.26 + 0.25, label = "4", size = 3.2) +
  annotate("text", x = cx + 0.25, y = 1.1 + sh_top + 7 * 0.26 + 0.25, label = "7", size = 3.2) +
  annotate("text", x = cx - 0.55, y = 0.75 + sh_top, label = "reg", size = 3.2, color = COL_REG) +
  annotate("text", x = cx + 0.25, y = 0.75 + sh_top, label = "trans", size = 3.2, color = COL_TRANS) +

  # ── Card 2 (middle): Validation (concentration, cell-type / duration) ──
  # 2-row layout: row 1 = concentration | cell type; row 2 = duration (centered)
  annotate("text", x = cx, y = 3.55 + sh_mid, label = "Validation", fontface = "bold", size = 3.4) +
  geom_point(aes(x = cx - 1.0, y = 2.75 + sh_mid), size = 1.5, color = "grey30") +
  geom_point(aes(x = cx - 0.8, y = 2.75 + sh_mid), size = 2.8, color = "grey30") +
  geom_point(aes(x = cx - 0.6, y = 2.75 + sh_mid), size = 4.2, color = "grey30") +
  annotate("text", x = cx - 0.8, y = 2.15 + sh_mid, label = "concentration", size = 3.0) +
  geom_polygon(data = hex_df(cx + 0.85, 2.75 + sh_mid, 0.42), aes(x, y), fill = "grey90", color = "grey30", linewidth = 0.5) +
  annotate("text", x = cx + 0.85, y = 2.15 + sh_mid, label = "cell type", size = 3.0) +
  geom_path(data = circle_df(cx, 1.35 + sh_mid, 0.42), aes(x, y), color = "grey30", linewidth = 0.5) +
  geom_segment(aes(x = cx, xend = cx, y = 1.35 + sh_mid, yend = 1.65 + sh_mid), color = "grey30", linewidth = 0.6) +
  geom_segment(aes(x = cx, xend = cx + 0.28, y = 1.35 + sh_mid, yend = 1.35 + sh_mid), color = "grey30", linewidth = 0.6) +
  annotate("text", x = cx, y = 0.75 + sh_mid, label = "duration", size = 3.0) +

  # ── Card 3 (bottom): BR-Auxin co-expression ────────────────────────────
  annotate("text", x = cx, y = 3.55 + sh_bot, label = "Co-expression", fontface = "bold", size = 3.4) +
  geom_polygon(data = circle_df(cx - 0.45, 2.1 + sh_bot, 0.65), aes(x, y), fill = COL_BR, alpha = 0.55, color = NA) +
  geom_polygon(data = circle_df(cx + 0.45, 2.1 + sh_bot, 0.65), aes(x, y), fill = COL_AUX, alpha = 0.55, color = NA) +
  annotate("text", x = cx - 0.55, y = 2.1 + sh_bot, label = "BR", size = 3.6, fontface = "bold", color = "white") +
  annotate("text", x = cx + 0.55, y = 2.1 + sh_bot, label = "Aux", size = 3.6, fontface = "bold", color = "white") +

  coord_fixed(xlim = c(0, 4), ylim = c(0, 14), expand = FALSE, clip = "off") +
  theme_void()

# ══════════════════════════════════════════════════════════════════════════
# Panel B - per-cell Spearman rho (concentration | duration), empirical bar
# ══════════════════════════════════════════════════════════════════════════
BR_SENSORS <- c(
  "ATH-br-reg-4hr-down", "ATH-br-reg-2hr-down", "ATH-br-reg-4hr-up",
  "ATH-br-trans-NegativeSignaling", "ATH-br-reg-2hr-up", "ATH-br-trans-TF-induced",
  "ATH-br-trans-Homeostasis", "ATH-br-trans-TF-repressed", "ATH-br-trans-TF",
  "ATH-br-trans-PositiveSignaling", "ATH-br-trans-Biosynthesis"
)
NEG_CONTROLS <- c("random1", "random2", "majortrend")
SENSORS <- c(BR_SENSORS, NEG_CONTROLS)

# Class prefix dropped from the label text (to fit the bigger axis font) -
# indicated instead via the colored strip in make_forest() / p_class_strip
SENSOR_LABELS <- c(
  "ATH-br-reg-4hr-down"             = "late-down",
  "ATH-br-reg-2hr-down"             = "early-down",
  "ATH-br-reg-4hr-up"               = "late-up",
  "ATH-br-trans-NegativeSignaling"  = "NegSignaling",
  "ATH-br-reg-2hr-up"               = "early-up",
  "ATH-br-trans-TF-induced"         = "TF-induced",
  "ATH-br-trans-Homeostasis"        = "Homeostasis",
  "ATH-br-trans-TF-repressed"       = "TF-repressed",
  "ATH-br-trans-TF"                 = "TF",
  "ATH-br-trans-PositiveSignaling"  = "PosSignaling",
  "ATH-br-trans-Biosynthesis"       = "Biosynthesis",
  "random1"                         = "random1",
  "random2"                         = "random2",
  "majortrend"                      = "majortrend"
)
SENSOR_CLASS <- c(
  "ATH-br-reg-4hr-down" = "reg", "ATH-br-reg-2hr-down" = "reg", "ATH-br-reg-4hr-up" = "reg",
  "ATH-br-trans-NegativeSignaling" = "trans", "ATH-br-reg-2hr-up" = "reg",
  "ATH-br-trans-TF-induced" = "trans", "ATH-br-trans-Homeostasis" = "trans",
  "ATH-br-trans-TF-repressed" = "trans", "ATH-br-trans-TF" = "trans",
  "ATH-br-trans-PositiveSignaling" = "trans", "ATH-br-trans-Biosynthesis" = "trans",
  "random1" = "neg", "random2" = "neg", "majortrend" = "neg"
)
CLASS_COLORS <- c(reg = COL_REG, trans = COL_TRANS, neg = "grey40")
LABEL_ORDER <- SENSOR_LABELS[SENSORS]

cat("Loading GSE212230 iSensors object...\n")
obj <- readRDS(file.path(BR_DIR, "out/GSE212230_iSensors_obj.rds"))
DefaultAssay(obj) <- "iSensors_mean"

meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  select(cell, sample, cell_type, condition)

sens <- FetchData(obj, vars = SENSORS) %>% tibble::rownames_to_column("cell")

df <- meta %>% inner_join(sens, by = "cell")
rm(obj); invisible(gc())

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
      sensor_label = factor(SENSOR_LABELS[sensor], levels = rev(LABEL_ORDER)),
      sensor_class = SENSOR_CLASS[sensor]
    )
  attr(res, "threshold") <- thr
  res
}

AXIS_TEXT_SIZE <- 9  # matches Panel E's y-axis (tissue) label size

make_forest <- function(res, thr, xlab, title, show_ylabels = TRUE) {
  rng <- range(c(res$lo, res$hi), na.rm = TRUE)
  pad <- diff(rng) * 0.08
  xlim <- c(rng[1] - pad, rng[2] + pad)

  p <- ggplot(res, aes(x = sensor_label, y = rho)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = c(-thr, thr), linetype = "dotted", color = "grey40", linewidth = 0.5) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.4, color = "grey50") +
    geom_point(aes(color = class), size = 2) +
    scale_color_manual(values = c(Up = "#d73027", Down = "#4575b4", "n.s." = "grey75",
                                   "Negative control" = "black"), guide = "none") +
    coord_flip(ylim = xlim) +
    labs(x = NULL, y = xlab, title = title) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(size = AXIS_TEXT_SIZE),
          axis.title.x = element_text(size = 9.5, lineheight = 0.9),
          plot.title = element_text(size = 9, face = "bold", lineheight = 1.0),
          plot.margin = margin(4, 6, 4, 2))

  if (show_ylabels) {
    p <- p + theme(axis.text.y = element_text(size = AXIS_TEXT_SIZE))
  } else {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  p
}

# Narrow colored strip (reg/trans/neg) placed left of the y-axis labels
make_class_strip <- function(res) {
  strip_df <- res %>% distinct(sensor_label, sensor_class)
  ggplot(strip_df, aes(x = 1, y = sensor_label, fill = sensor_class)) +
    geom_tile() +
    scale_fill_manual(values = CLASS_COLORS, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(4, 1, 4, 2))
}

conc_cells <- df %>%
  filter(condition %in% c("BRZ", "Control", "BL_8h")) %>%
  mutate(level = case_when(condition == "BRZ" ~ 1, condition == "Control" ~ 2, condition == "BL_8h" ~ 3))
res_conc <- classify_empirical(spearman_per_sensor(conc_cells, "level", SENSORS))
thr_conc <- attr(res_conc, "threshold")

p_conc <- make_forest(res_conc, thr_conc,
                       "Spearman rho\n(score vs BR level)",
                       "BR concentration response\n(low → intermediate → high)",
                       show_ylabels = TRUE)

dur_cells <- df %>%
  filter(condition %in% c("BL_0.5h", "BL_1h", "BL_2h", "BL_4h", "BL_8h")) %>%
  mutate(hours = case_when(
    condition == "BL_0.5h" ~ 0.5, condition == "BL_1h" ~ 1, condition == "BL_2h" ~ 2,
    condition == "BL_4h" ~ 4, condition == "BL_8h" ~ 8))
res_dur <- classify_empirical(spearman_per_sensor(dur_cells, "hours", SENSORS))
thr_dur <- attr(res_dur, "threshold")

p_duration <- make_forest(res_dur, thr_dur,
                           "Spearman rho\n(score vs BL exposure duration)",
                           "BR duration response\n(0.5, 1, 2, 4, 8 hrs)",
                           show_ylabels = FALSE)

p_strip <- make_class_strip(res_conc)
p_panelB <- p_strip + p_conc + p_duration + plot_layout(widths = c(0.02, 1, 1))
rm(df); invisible(gc())

# ══════════════════════════════════════════════════════════════════════════
# Panel C - BR (reg-up) vs Auxin (ARF), longitudinal root layout pair
# (= original "Panel A" from 14-BR-ARF-spatial-anticorrelation.R)
# ══════════════════════════════════════════════════════════════════════════
data("ggPm.At.longroot.longitudinal", package = "ggRootCellAtlas", envir = .GlobalEnv)

br_wide  <- read.csv(file.path(BR_DIR, "out/shahan_root_BR_scores.csv"),  check.names = FALSE)
aux_wide <- read.csv(file.path(BR_DIR, "out/shahan_auxin_scores.csv"),   check.names = FALSE)

get_score_df <- function(row, drop_cols) {
  cols <- setdiff(names(row), drop_cols)
  data.frame(Atlas = gsub("-", "_", cols),
             Gene  = as.numeric(row[1, cols]),
             stringsAsFactors = FALSE)
}
br_df  <- get_score_df(br_wide [br_wide$sensor  == "ATH-br-reg-4hr-up", ], c("sensor", "sensor_label", "inverse"))
aux_df <- get_score_df(aux_wide[aux_wide$sensor == "ATH-aux-trans-ARF", ], c("sensor", "aux_label"))

make_long_panel <- function(score_df, colours, legend_title, title_text) {
  merged <- ggPlantmap.merge(ggPm.At.longroot.longitudinal, score_df, "Atlas")
  lims   <- quantile(score_df$Gene, c(0.02, 0.98), na.rm = TRUE)
  ggPlantmap.heatmap(merged, Gene) +
    scale_fill_gradientn(colours = colours, limits = lims, oob = squish, na.value = "grey92", name = legend_title) +
    labs(title = title_text) +
    theme(plot.title        = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.position   = "right",
          legend.key.height = unit(1.4, "cm"),
          legend.title      = element_text(size = 8),
          legend.text       = element_text(size = 7),
          plot.margin       = margin(t = 16, r = 4, b = 4, l = 4))
}

p_br_long  <- make_long_panel(br_df,  colours = c("white", "#C7E9C0", "#74C476", "#238B45", "#00441B"),
                               legend_title = "BR\nactivity", title_text = "BR response\n(reg-late-up)")
p_aux_long <- make_long_panel(aux_df, colours = c("white", "#FFFFB2", "#FEB24C", "#F03B20", "#BD0026"),
                               legend_title = "Auxin\nARF", title_text = "Auxin Response\n(ARF)")

p_panelC <- wrap_elements(full = p_br_long) + wrap_elements(full = p_aux_long)

# ══════════════════════════════════════════════════════════════════════════
# Panel F - BR (reg-late-up) vs Auxin (ARF), epidermis layout pair
# (= original "Panel B" from 14-BR-ARF-spatial-anticorrelation.R)
# ══════════════════════════════════════════════════════════════════════════
REALISTIC_DIR <- "../RealisticLayouts"

br_epi_scores <- read.csv(file.path(BR_DIR, "out/shahan_epidermis_BR_scores.csv"), stringsAsFactors = FALSE)
br_epi_df <- br_epi_scores %>%
  filter(sensor == "ATH-br-reg-4hr-up") %>%
  select(Cell_type = atlas_stage, Gene = br_activity)

epi_cols <- c("Atrichoblast-d", "Atrichoblast-e1", "Atrichoblast-e2",
              "Atrichoblast-m1", "Atrichoblast-m2", "Atrichoblast-t",
              "Trichoblast-d", "Trichoblast-e1", "Trichoblast-e2",
              "Trichoblast-m1", "Trichoblast-m2", "Trichoblast-t")
aux_epi_row <- aux_wide[aux_wide$sensor == "ATH-aux-trans-ARF", ]
aux_epi_df  <- data.frame(
  Cell_type = epi_cols,
  Gene      = as.numeric(aux_epi_row[1, epi_cols]),
  stringsAsFactors = FALSE
)

epi_layout <- read.csv(file.path(REALISTIC_DIR, "out/new_ggPlantmap_epidermis.csv"), stringsAsFactors = FALSE)
epi_layout$x <- as.numeric(epi_layout$x)
epi_layout$y <- as.numeric(epi_layout$y)

make_epi_panel <- function(score_df, colours, legend_title, title_text) {
  merged <- epi_layout %>% left_join(score_df, by = "Cell_type")
  lims   <- quantile(score_df$Gene, c(0.02, 0.98), na.rm = TRUE)
  ggplot(merged, aes(x = x, y = y, group = ROI.id, fill = Gene)) +
    geom_polygon(colour = "grey70", linewidth = 0.15) +
    scale_fill_gradientn(colours = colours, limits = lims, oob = squish, na.value = "grey92", name = legend_title) +
    coord_fixed() +
    labs(title = title_text) +
    theme_void(base_size = 10) +
    theme(plot.title        = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.position   = "right",
          legend.key.height = unit(1.2, "cm"),
          legend.title      = element_text(size = 9.5),
          legend.text       = element_text(size = 8.5),
          plot.margin       = margin(t = 16, r = 4, b = 4, l = 4))
}

p_br_epi  <- make_epi_panel(br_epi_df,  colours = c("white", "#C7E9C0", "#74C476", "#238B45", "#00441B"),
                             legend_title = "BR\nactivity", title_text = "BR response\n(reg-late-up)")
p_aux_epi <- make_epi_panel(aux_epi_df, colours = c("white", "#FFFFB2", "#FEB24C", "#F03B20", "#BD0026"),
                             legend_title = "Auxin\nARF", title_text = "Auxin Response\n(ARF)")

p_panelF <- wrap_elements(full = p_br_epi) + wrap_elements(full = p_aux_epi)

# ══════════════════════════════════════════════════════════════════════════
# Panel D - Quadrant scatter (= original "Panel C")
# Panel E - Tissue composition (= original "Panel D")
# ══════════════════════════════════════════════════════════════════════════
BR_THRESH  <- 1.0
ARF_THRESH <- 0.1

joined <- read.csv(file.path(BR_DIR, "out/shahan_BR_Aux_joined.csv"), stringsAsFactors = FALSE)

dat_raw <- joined %>%
  filter(br_sensor == "ATH-br-reg-4hr-up", aux_sensor == "ATH-aux-trans-ARF") %>%
  select(cluster, br_score, aux_score, tissue) %>%
  distinct() %>%
  mutate(
    br_group  = ifelse(br_score  > BR_THRESH,  "BR high",  "BR low"),
    arf_group = ifelse(aux_score > ARF_THRESH, "ARF high", "ARF low"),
    quadrant  = paste(br_group, arf_group, sep = " / "),
    pattern   = case_when(
      quadrant %in% c("BR high / ARF low", "BR low / ARF high") ~ "Mutual exclusion",
      TRUE ~ "Co-localization"
    )
  )

dat <- dat_raw %>%
  mutate(tissue = ifelse(tissue == "Cell cycle", "Other", tissue),
         tissue = ifelse(tissue == "Root tip",   "Stem cell niche", tissue))

sp    <- cor.test(dat$br_score, dat$aux_score, method = "spearman")
r_obs <- unname(sp$estimate); p_obs <- sp$p.value
set.seed(42)
r_ci <- quantile(
  replicate(2000, { i <- sample(nrow(dat), replace = TRUE)
    cor(dat$br_score[i], dat$aux_score[i], method = "spearman") }),
  c(0.025, 0.975))

ctab <- table(BR_high  = dat$br_score  > BR_THRESH, ARF_high = dat$aux_score > ARF_THRESH)
fish <- fisher.test(ctab)

tissue_colors <- c(
  "Epidermis" = "#E67E22", "Cortex" = "#27AE60", "Endodermis" = "#2980B9",
  "Stele" = "#8E44AD", "LRC/Columella" = "#C0392B", "Stem cell niche" = "#16A085", "Other" = "#BDC3C7"
)
quad_colors <- c(
  "BR high / ARF low"  = "#7B2D8B",
  "BR low / ARF high"  = "#27AE60",
  "BR high / ARF high" = "#8B5E3C",
  "BR low / ARF low"   = "#E67E22"
)

stats_label <- sprintf("Fisher OR = %.3f, p = %.5f", fish$estimate, fish$p.value)

qn <- table(dat$quadrant)
xr <- range(dat$br_score); yr <- range(dat$aux_score)
dx <- diff(xr) * 0.03; dy <- diff(yr) * 0.03
corners <- data.frame(
  label = c(
    paste0("BR high\nARF high\nn=", qn["BR high / ARF high"]),
    paste0("BR low\nARF high\nn=",  qn["BR low / ARF high"]),
    paste0("BR high\nARF low\nn=",  qn["BR high / ARF low"]),
    paste0("BR low\nARF low\nn=",   qn["BR low / ARF low"])
  ),
  x = c(xr[2] - dx, xr[1] + dx, xr[2] - dx, xr[1] + dx),
  y = c(yr[2] - dy, yr[2] - dy, yr[1] + dy, yr[1] + dy),
  hjust = c(1, 0, 1, 0), vjust = c(1, 1, 0, 0)
)

p_panelD <- ggplot(dat, aes(br_score, aux_score)) +
  annotate("rect", xmin = BR_THRESH, xmax = Inf, ymin = -Inf, ymax = ARF_THRESH, fill = "#7B2D8B", alpha = 0.07) +
  annotate("rect", xmin = -Inf, xmax = BR_THRESH, ymin = ARF_THRESH, ymax = Inf, fill = "#27AE60", alpha = 0.07) +
  annotate("rect", xmin = BR_THRESH, xmax = Inf, ymin = ARF_THRESH, ymax = Inf, fill = "#8B5E3C", alpha = 0.07) +
  geom_vline(xintercept = BR_THRESH, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = ARF_THRESH, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = tissue), size = 3, alpha = 0.9) +
  geom_text(data = corners, aes(x, y, label = label, hjust = hjust, vjust = vjust),
            size = 2.7, colour = "grey30", lineheight = 0.9, inherit.aes = FALSE) +
  annotate("text", x = BR_THRESH, y = -Inf, label = paste0("BR = ", BR_THRESH),
           hjust = 0.5, vjust = -0.4, size = 2.8, colour = "grey50") +
  annotate("text", x = xr[2], y = ARF_THRESH, label = paste0("ARF = ", ARF_THRESH),
           hjust = 1.05, vjust = -0.4, size = 2.8, colour = "grey50") +
  scale_colour_manual(values = tissue_colors, name = "Tissue") +
  labs(x = "BR response (reg-late-up) cluster mean", y = "Auxin Response (ARF) cluster mean",
       title = "Quadrant analysis", subtitle = stats_label) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 9, colour = "grey20", family = "mono"),
        legend.position = "bottom", legend.key.size = unit(0.4, "cm"), legend.text = element_text(size = 8))

tissue_quad <- dat %>%
  count(tissue, quadrant) %>%
  group_by(tissue) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

tissue_order <- tissue_quad %>% filter(quadrant == "BR high / ARF low") %>% arrange(desc(pct)) %>% pull(tissue)
tissue_order <- c(tissue_order, setdiff(unique(dat$tissue), tissue_order))

tissue_fish <- dat %>%
  group_by(tissue) %>%
  filter(n() >= 4) %>%
  summarise(
    n_tot = n(), n_excl = sum(pattern == "Mutual exclusion"), pct_excl = n_excl / n_tot * 100,
    p_fish = tryCatch({
      tt <- table(BR_high = br_score > BR_THRESH, ARF_high = aux_score > ARF_THRESH)
      if (all(dim(tt) == c(2, 2))) fisher.test(tt)$p.value else NA_real_
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(sig = case_when(p_fish < 0.05 ~ "*", p_fish < 0.1 ~ ".", TRUE ~ ""),
         tissue = factor(tissue, levels = rev(tissue_order)))

p_panelE <- ggplot(
    tissue_quad %>% mutate(tissue = factor(tissue, levels = rev(tissue_order)),
                            quadrant = factor(quadrant, levels = names(quad_colors))),
    aes(pct, tissue, fill = quadrant)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = quad_colors, name = NULL) +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), labels = function(x) paste0(x, "%")) +
  labs(x = "% of tissue clusters", y = NULL, title = "Tissue composition") +
  guides(fill = guide_legend(nrow = 2)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        axis.text = element_text(size = AXIS_TEXT_SIZE),
        panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(0.35, "cm"), legend.text = element_text(size = 8))

# ══════════════════════════════════════════════════════════════════════════
# Combine: (A | B | C) / (D | E | F) and save
# ══════════════════════════════════════════════════════════════════════════
# A+B width fraction must match D+E width fraction so that C and F (each a
# BR | Auxin pair) line up in the same columns across both rows.
#
# Absolute widths (in), total canvas 12.64in:
#   A = 2.04 (-15% from 2.4)      D = 4.94 (nudged up to keep D+E = A+B)
#   B = 6.8  (unchanged)          E = 3.9  (-15% from 4.6, more gap before F)
#   C = 3.8  (unchanged)          F = 3.8  (= C, must match C for alignment)
w_A <- 2.04; w_B <- 6.8; w_C <- 3.8; w_E <- 3.9; w_F <- 3.8
w_D <- (w_A + w_B) - w_E
TOTAL_W <- w_A + w_B + w_C

p_row1 <- wrap_elements(full = p_schematic) + wrap_elements(full = p_panelB) + wrap_elements(full = p_panelC) +
  plot_layout(widths = c(w_A, w_B, w_C))
p_row2 <- p_panelD + p_panelE + wrap_elements(full = p_panelF) +
  plot_layout(widths = c(w_D, w_E, w_F))

p_final <- (p_row1 / p_row2) +
  plot_layout(heights = c(0.6, 0.4)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(output_dir, "Figure8.pdf"), plot = p_final, width = TOTAL_W, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure8.png"), plot = p_final, width = TOTAL_W, height = 11, dpi = 300, bg = "white")
ggsave(file.path(output_dir, "Figure8.svg"), plot = p_final, width = TOTAL_W, height = 11, dpi = 300, bg = "white")
cat("Saved Figure8.pdf/.png/.svg\n")

write.csv(res_conc %>% select(sensor, sensor_label, rho, lo, hi, n, p.value, p_adj, class),
          file.path(output_dir, "Figure8_concentration_stats.csv"), row.names = FALSE)
write.csv(res_dur %>% select(sensor, sensor_label, rho, lo, hi, n, p.value, p_adj, class),
          file.path(output_dir, "Figure8_duration_stats.csv"), row.names = FALSE)
cat("Saved stats CSVs.\n")

cat("\nEmpirical threshold (concentration axis):", round(thr_conc, 4), "\n")
cat("Empirical threshold (duration axis):     ", round(thr_dur, 4), "\n")
cat("BR-ARF Spearman r =", round(r_obs, 3), " Fisher OR =", round(fish$estimate, 3), "\n")
