library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(ggpubr)
library(ggrepel)
set.seed(100)


in_path  <- here("Manuscript-Figures", "in")
out_path <- here("Manuscript-Figures", "out")


### Figure 3C (Heatmap)

# Read iSensors scores and predefined row order
matrix_data_new <- as.matrix(read.table(file.path(in_path, "Figure3C_bulk_iSensors_matrix.txt"),header=TRUE,row.names=1)) # says first column are rownames     #TO BE UPDATED
row_order_names <- read.table(file.path(in_path, "Figure3C_HM_row_order_Z.txt"), header = TRUE, row.names = 1)[[1]]
col_order_vector <- read.table(file.path(in_path, "Figure3C_HM_col_order_Z.txt"), header = TRUE, row.names = 1)[[1]]


# Apply z-score normalization
mat_T <- t(matrix_data_new)
df_zscore <- as.data.frame(scale(mat_T))
matrix_data_new <- t(df_zscore)


# Define breaks and colors
mycols <- colorRamp2(breaks = c(min(matrix_data_new), 0, max(matrix_data_new)),
                     colors = c("#4575b4", "#ffffbf", "#d73027"))


# Prepare heatmap annotation (panel types)
Panel_name <- rownames(matrix_data_new)
Panel_type <- c()
i = 1
for (i in (1:nrow(matrix_data_new))) {
  if (grepl("cis", Panel_name[i], ignore.case = TRUE)) {
    Panel_type <- c(Panel_type, 'Cis')
  } else {
    if (grepl("trans", Panel_name[i], ignore.case = TRUE)) {
      Panel_type <- c(Panel_type, 'Trans')
    } else {
      if (grepl("reg", Panel_name[i], ignore.case = TRUE)) {
        Panel_type <- c(Panel_type, 'Reg')
      } else {
        Panel_type <- c(Panel_type, 'Negative')
      }
    }
  }
}

col_2 = list(
  Panel_type = c("Cis" = "#fdbf6f", "Trans" = "#ff7f00", "Reg" = "#b15928", "Negative" = "gray"))   # Define colors for each levels of qualitative variables
ha_2 <- rowAnnotation(
  Panel_type = Panel_type,
  col = col_2,
  annotation_legend_param = list(
      title = "Panel type",
      title_gp = gpar(fontsize = 14, fontface = "bold"), # Title font size
      labels_gp = gpar(fontsize = 12) # Labels font size
    )
)

# Plot heatmap
png(
  filename = file.path(out_path, "Fig3C_heatmap.png"),
  width = 31,
  height = 28, 
  units = "cm",
  res = 300
)
Heatmap(matrix_data_new,
        right_annotation = ha_2,
        column_order = col_order_vector,
        row_order = row_order_names,
        heatmap_legend_param = list(
          title = "Z-score", # Legend Title
          title_gp = gpar(fontsize = 14), # Title font size
          labels_gp = gpar(fontsize = 10), # Label font size
          direction = "vertical"         # "horizontal" or "vertical"
        ),
        col = mycols,
        column_title = "Pseudo cells",
        column_title_side = "bottom",
        column_names_side = "bottom",
        row_title = "iSensors",
        row_title_side = "left",
        row_names_side = "left",
        row_names_gp = gpar(fontsize = 7), # Text size for row names
        column_names_gp = gpar(fontsize = 7) # Text size for row names
)
dev.off()


### Figure 3C (forest plot)

# Read linear modeling results
res_limma <- read.table(file.path(in_path, "Figure3C_bulk-limma-results.txt"),header=TRUE,row.names=1)


# Prepare data to build a forest plot for effect size
res_forest <- res_limma %>%
  mutate(
    se = abs(effect / t),
    lo = effect - 1.96 * se,
    hi = effect + 1.96 * se,
    sig = case_when(
      p_adj >= 0.05 ~ "n.s.",
      effect < -0.05 ~ "down",
      effect > 0.05  ~ "up",
      TRUE           ~ "n.s."
    )
  )


# Change row order according to the heatmap
row_order_df <- tibble(iSensor = row_order_names)

setdiff(row_order_df$iSensor, res_forest$iSensor)   # Check point
setdiff(res_forest$iSensor, row_order_df$iSensor)   # Check point

res_forest <- left_join(row_order_df, res_forest, by = "iSensor")


# Generate forest plot
my_order <- rev(row_order_df$iSensor)
res_forest <- res_forest %>%
  mutate(iSensor = factor(iSensor, levels = my_order))

forest_plot <- ggplot(res_forest, aes(y = iSensor, x = effect)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  geom_errorbarh(
    aes(xmin = lo, xmax = hi),
    linewidth = 0.5,
    height = 0.18,
    color = "grey50"
  ) +
  geom_point(aes(color = sig), size = 2) +
  scale_color_manual(values = c("up"="#d73027","down"="#4575b4","n.s."="grey75"),
                     guide = "none") +                                             # Removes the legend for colors
  labs(
    x = "Effect size (Auxin – Control)",
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 9),
    #    axis.text.y  = element_blank(),       # Removes the text labels from the vertical Y-axis
    axis.ticks.y = element_blank(),
    plot.margin = margin(4, 6, 4, 2)
  )
forest_plot

ggsave(filename = file.path(out_path, "Fig3C_forest-plot.pdf"), plot = forest_plot, width = 8, height = 10, units = "in")


### Figure 3D

# Read iSensors score matrix
iSensor_test_3 <- as.matrix(read.table(file.path(in_path, "Figure3D_bulk_iSensors_matrix_concentration.txt"),header=TRUE,row.names=1))
fenoVecNotMean_3 <- as.numeric(iSensor_test_3["fenoVecNotMean_3", ])


# Specify iSensors to be visualized
vec <- c("ATH-aux-trans-PolarAuxinTransport", "ATH-aux-trans-IAA", "ATH-aux-reg-IR8-ARF2-up", "ATH-aux-trans-Receptors", "ATH-aux-cis-DR5-TGTCGG", "random3")

# Plot individual correlations (iSensors vs concentration)
for (element in vec) {
  panelSignal <- iSensor_test_3[element, ]
  spearman_test_result <- cor.test(fenoVecNotMean_3, panelSignal, method = "spearman")
  
  df <- data.frame(
    x = rank(fenoVecNotMean_3),  # X axis
    y = rank(panelSignal)      # Y axis
  )
  
  plot_title <- gsub("^ATH-aux-", "", element)
  
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = "black", alpha = 1, size = 2) +  # dots
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(method = "spearman", 
             label.x = 8,         # X coordinate for text
             label.y = 30,        # Y coordinate for text
             size = 4) +          # Annotation font size
    coord_cartesian(xlim = c(7.5, 26), ylim = c(0, 30)) +
    labs(
      title = plot_title,
      x = "Rank (Concentration)",
      y = "Rank (iSensors)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        size = 24,
      ),
      axis.text = element_text(size = 18),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks = element_line(color = "black", , linewidth = 0.8),
      axis.ticks.length = unit(0.2, "cm")
    )
  p
  
  file_name <- paste0("Fig3D_", plot_title, ".pdf")
  ggsave(filename = file.path(out_path, file_name), plot = p, width = 5, height = 3, units = "in")
}


### Figure 3E

# Read correlation data
df_stat_1 <- read.table(file.path(in_path, "Figure3E_Correlation_stat_duration.txt"), header=TRUE, row.names=1)
df_stat_3 <- read.table(file.path(in_path, "Figure3E_Correlation_stat_concentration.txt"), header=TRUE, row.names=1)


# Prepare a summary table
df_stat_1 <- rownames_to_column(df_stat_1, "iSensor")
df_stat_3 <- rownames_to_column(df_stat_3, "iSensor")

setdiff(df_stat_1$iSensor, res_forest$iSensor)
setdiff(res_forest$iSensor, df_stat_1$iSensor)
setdiff(df_stat_3$iSensor, res_forest$iSensor)
setdiff(res_forest$iSensor, df_stat_3$iSensor)

df_stat <- left_join(df_stat_1, df_stat_3, by = "iSensor")
df_stat <- left_join(df_stat, res_forest, by = "iSensor")

rownames(df_stat) <- df_stat$iSensor
rownames(df_stat) <- gsub("ATH-aux-", "", rownames(df_stat), fixed = TRUE)
df_stat$iSensor <- rownames(df_stat)


# Group the data according to a condition
cond <- character()
i <- 1
for (i in (1:nrow(df_stat))) {
  if (df_stat$bonferroni1[[i]] < 0.05 & df_stat$bonferroni3[[i]] < 0.05) {
    cond <- c(cond, 'Both')
  } else {
    if (df_stat$bonferroni1[[i]] < 0.05) {
      cond <- c(cond, 'Duration')
    } else {
      if (df_stat$bonferroni3[[i]] < 0.05) {
        cond <- c(cond, 'Concentration')
      } else {
        cond <- c(cond, 'Not significant')
      }
    }
  }
  i = i + 1
}

df_stat$condition1  <- cond

# Plot Rho (Duration) vs Rho (Concentration)
df <- data.frame(
  x = df_stat$Rho1,  # X axis
  y = df_stat$Rho3,      # Y axis
  CorrSign = factor(df_stat$condition1),
  Effect = factor(df_stat$sig)
)

rownames(df) <- df_stat$iSensor

df <- rownames_to_column(df, "iSensor")

sensors_to_label <- c("trans-PolarAuxinTransport", "trans-Transport", "trans-IAA", "trans-Synthesis", "trans-A-ARF", "reg-DR5-ARF1-up", "trans-ARF", "reg-DR5-ARF8-down", "reg-DR5-ARF2-down")

scarret_plot <- ggplot(df, aes(x = x, y = y, shape = CorrSign)) +
  geom_point(aes(color = Effect), size = 7, alpha = 1) +  # dots
  geom_text_repel(
    aes(label = ifelse(iSensor %in% sensors_to_label, as.character(iSensor), "")),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50"
  ) +
  scale_color_manual(values = c("up" = "#d73027", "down" = "#4575b4", "n.s." = "grey75")) +
  geom_hline(yintercept = 0, color = "black") +   # Add horizontal line
  geom_vline(xintercept = 0, color = "black") +
  labs(
    x = "Rho (Duration)",
    y = "Rho (Concentration)"
  ) +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 18)
  )

ggsave(filename = file.path(out_path, "Fig3E.pdf"), plot = scarret_plot, width = 10, height = 10, units = "in")




