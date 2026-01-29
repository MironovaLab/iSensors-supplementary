## Spearman correlation calculator

options("scipen" = 100, "digits" = 4)

library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(stringr)

# get object
# seurat_obj <- readRDS("C:/Users/cvand/OneDrive - Radboud Universiteit/Y3Q4 - Internship/Testing/iSensor_Shahan_output.rds")
seurat_obj <- readRDS("D:/MyPAPERS/iSensors/Data/iSensor_Shahan_output_v0.4_6meta_panels.rds")

# remove random panels
seurat_obj@assays

# assays
assay_list <- c("iSensor_mean_normed")
# assay_list <- c("iSensor_mean_normed")

# ranks known data
# auxin_rank <- c(
#   "QC" = 1,
#   "CSC" = 4,
#   "LRP" = 4,
#   "LRCEI" = 4,
#   "SI" = 4,
#   "CEI" = 4,
#   "Columella" = 9,
#   "CSCD" = 9,
#   "Young LRC" = 9,
#   "Protoxylem-m1" = 9,
#   "Metaxylem-m1" = 9,
#   "Protoxylem-m2" = 12.5,
#   "Metaxylem-m2" = 12.5,
#   "LRC" = 14,
#   "Dying LRC" = 15,
#   "Procambium-m1" = 16,
#   "Procambium-m2" = 17,
#   "Endodermis-m1" = 19,
#   "PPP-m1" = 19,
#   "XPP-m1" = 19,
#   "Endodermis-m2" = 22,
#   "PPP-m2" = 22,
#   "XPP-m2" = 22,
#   "Cortex-m1" = 24,
#   "Cortex-m2" = 25,
#   "Atrichoblast-m1" = 26.5,
#   "Trichoblast-m1" = 26.5,
#   "Atrichoblast-m2" = 28.5,
#   "Trichoblast-m2" = 28.5
# )

# new ranks
auxin_rank <- c(
  "QC" = 1,
  "CSC" = 2.5,
  "LRP" = 2.5,
  "LRCEI" = 8,
  "SI" = 8,
  "CEI" = 8,
  "Columella" = 5,
  "CSCD" = 5,
  "Young LRC" = 5,
  "Protoxylem-m1" = 10.5,
  "Metaxylem-m1" = 10.5,
  "Protoxylem-m2" = 13,
  "Metaxylem-m2" = 13,
  "LRC" = 13,
  #"Dying LRC" = 15,
  "Procambium-m1" = 15.5,
  "Procambium-m2" = 15.5,
  "Endodermis-m1" = 18,
  "PPP-m1" = 18,
  "XPP-m1" = 18,
  "Endodermis-m2" = 21,
  "PPP-m2" = 21,
  "XPP-m2" = 21,
  "Cortex-m1" = 23.5,
  "Cortex-m2" = 23.5,
  "Atrichoblast-m1" = 25.5,
  "Trichoblast-m1" = 25.5,
  "Atrichoblast-m2" = 27.5,
  "Trichoblast-m2" = 27.5
)

auxin_rank <- c(
  "QC" = 1,
  "CSC" = 2,
  "LRP" = 2,
  "LRCEI" = 4,
  "SI" = 4,
  "CEI" = 4,
  "Columella" = 3,
  "CSCD" = 3,
  "Young LRC" = 3,
  "Protoxylem-m1" = 5,
  "Metaxylem-m1" = 5,
  "Protoxylem-m2" = 6,
  "Metaxylem-m2" = 6,
  "LRC" = 6,
  "Procambium-m1" = 7,
  "Procambium-m2" = 7,
  "Endodermis-m1" = 8,
  "PPP-m1" = 8,
  "XPP-m1" = 8,
  "Endodermis-m2" = 9,
  "PPP-m2" = 9,
  "XPP-m2" = 9,
  "Cortex-m1" = 10,
  "Cortex-m2" = 10,
  "Atrichoblast-m1" = 11,
  "Trichoblast-m1" = 11,
  "Atrichoblast-m2" = 12,
  "Trichoblast-m2" = 12
)

#subset some panels for testing
# panels_test <- c(
  # "AT-aux-cis-DR5(TGTCGG)",
  # "AT-aux-cis-IR8(ARF3)",
  # "AT-aux-cistrans-DR5(ARF8)down",
  # "AT-aux-cistrans-IR8(ARF6)up",
  # "AT-aux-cistrans-DR5(ARF6)up",
  # "AT-aux-trans-A-ARF",
  # "AT-aux-trans-ARF",
  # "AT-aux-trans-Effluxinflux",
  # "AT-aux-trans-IAA",
  # "AT-aux-trans-Synthesis"
  # "AT-aux-trans-Transport",
  # "majortrend",
  # "random1",
  # "random2"
# )

# cleanpanel
panel_df <- tibble(panel = rownames(seurat_obj)) %>%
  mutate(
    reg_class = case_when(
      str_detect(panel, "R2D2") ~ "meta",
      str_detect(panel, "cistrans") ~ "cistrans",
      str_detect(panel, "cist") ~ "cist",
      str_detect(panel, "cis") ~ "cis",
      str_detect(panel, "trans") ~ "trans",
      TRUE ~ "other"
    ),
    reg_class = if_else(
      panel %in% c(
        "AT-aux-cist-IR8(ARF5(2))",
        "AT-aux-cist-IR8(ARF5(1))",
        "AT-aux-cist-DR5(ARF5(2))",
        "AT-aux-cist-DR5(ARF5(1))"
      ),
      "cis",
      reg_class
    ),
    clean_panel = str_replace(panel, "^AT-aux-[^-]+-", ""),
    clean_panel = case_when(
      str_detect(clean_panel, "^EffluxInflux") ~ str_replace(clean_panel, "^EffluxInflux", "AUX/LAX/PIN"),
      str_detect(clean_panel, "^R2D2")         ~ str_replace(clean_panel, "^R2D2", "ARFxIAA"),
      TRUE                                     ~ clean_panel
    ),
    clean_panel = paste0(reg_class, "-", clean_panel)
  )


all_spearman_df <- list()


#folder_out <- "Shahan/UPDATE/Gene_panel_spearman_plots"
folder_out <- "02-single-cell-data-analysis/out/Spearman-plots"
if (!dir.exists(folder_out)) dir.create(folder_out)

# loop
for (assay in assay_list) {
  avg_exp <- AverageExpression(seurat_obj, assays = assay, layer = "counts")
  mat <- avg_exp[[assay]]
  
  cell_types <- names(auxin_rank)
  filtered_mat <- mat[, cell_types]
  df_filtered <- as.data.frame(as.matrix(filtered_mat))
  
  panel_ranks <- NULL
  for (panel in rownames(df_filtered)) {
    one_panel <- df_filtered[panel, ]
    one_panel_ranked <- rank(-one_panel, ties.method = "average")
    one_panel_ranked <- as.data.frame(one_panel_ranked)
    colnames(one_panel_ranked) <- c(panel)
    if (is.null(panel_ranks)) {
      panel_ranks <- one_panel_ranked
    } else {
      panel_ranks <- cbind(panel_ranks, one_panel_ranked)
    }
  }
  
  cor_stats <- t(sapply(colnames(panel_ranks), function(panel) {
    test <- cor.test(panel_ranks[, panel], auxin_rank[rownames(panel_ranks)], method = "spearman", exact = FALSE)
    c(rho = test$estimate, p_value = test$p.value)
  }))
  spearman_df <- as.data.frame(cor_stats)
  spearman_df <- spearman_df[order(spearman_df$rho.rho), ]
  
  spearman_df$panel <- rownames(spearman_df)
  spearman_df$assay <- assay
  
  # Join with clean_panel
  spearman_df <- left_join(spearman_df, panel_df[, c("panel", "clean_panel")], by = "panel")
  
  all_spearman_df[[assay]] <- spearman_df
  
  for (panel in colnames(panel_ranks)) {
    axis_x <- auxin_rank
    axis_y <- panel_ranks[, panel]
    plot_cor <- data.frame(auxin_ranks = axis_x, panel_ranks = axis_y)
    
    rho_val <- round(spearman_df$rho.rho[spearman_df$panel == panel], 3)
    clean_panel_val <- spearman_df$clean_panel[spearman_df$panel == panel]
    
    p <- ggplot(plot_cor, aes(x = auxin_ranks, y = panel_ranks)) +
      geom_point(color = "#1f78b4", size = 3) +
      geom_smooth(method = "lm", color = "grey", se = FALSE) +
      scale_x_continuous(
        limits = c(1, 13),
        breaks = seq(2, 12, by = 2)
      ) +
      scale_y_continuous(
        limits = c(1, 30),
        breaks = seq(5, 30, by = 5)
      ) +
      labs(
        x = "Ground truth rank",
        y = "Predicted rank, iSensors"
      ) +
      annotate(
        "text", x = 12, y = max(axis_y),
        label = paste("Rho = ", rho_val),
        hjust = 1.2, vjust = 1.2,
        size = 4, color = "black"
      ) +
      NoLegend(TRUE) +
      theme_minimal()+
      theme(
        panel.grid = element_blank(),       # remove all grid lines
        plot.title = element_blank(),
        axis.line  = element_line(),        # add axis lines
        axis.ticks = element_line()
      )

    filename <- paste0(folder_out, "/", assay, "_", panel, "_spearman.svg")
    ggsave(filename, plot = p, width = 4, height = 3)
  }
  
  message("Finished: ", assay)
}




#################################################################################
## Make barplot for all panels -----
combined_df <- bind_rows(all_spearman_df)
rownames(combined_df) <- NULL

combined_df_sub <- subset(combined_df, combined_df$assay == "iSensor_mean_normed")

# Define significance
combined_df_sub$significance <- ifelse(
  combined_df_sub$rho.rho > 0.45 & combined_df_sub$p_value < 0.05,
  "significant", "not significant"
)



# Add reg class and clean panel
combined_df_sub <- combined_df_sub %>%
  mutate(
    reg_class = case_when(
      str_detect(panel, "R2D2") ~ "meta",
      str_detect(panel, "cistrans") ~ "reg",
      str_detect(panel, "cist") ~ "cis",
      str_detect(panel, "cis") ~ "cis",
      str_detect(panel, "trans") ~ "trans",
      TRUE ~ "other"
    ),
    reg_class = if_else(
      panel %in% c(
        "AT-aux-cist-IR8(ARF5(2))",
        "AT-aux-cist-IR8(ARF5(1))",
        "AT-aux-cist-DR5(ARF5(2))",
        "AT-aux-cist-DR5(ARF5(1))"
      ),
      "cis",
      reg_class
    ),
    clean_panel = str_replace(panel, "^AT-aux-[^-]+-", ""),
    clean_panel = case_when(
      str_detect(clean_panel, "^EffluxInflux") ~ str_replace(clean_panel, "^EffluxInflux", "PAT"),
      #str_detect(clean_panel, "^R2D2")         ~ str_replace(clean_panel, "^R2D2", "ARFxIAA"),
      TRUE                                     ~ clean_panel
    ),
    clean_panel = paste0(clean_panel)
  )

# Optional: set reg_class for panels that don't match cis/trans/cistrans
combined_df_sub$reg_class[is.na(combined_df_sub$reg_class)] <- "other"

combined_df_sub <- combined_df_sub %>%
  filter(!str_detect(reg_class, regex("meta", ignore_case = TRUE)))
# If your meta panels match a different pattern, edit "Response" accordingly
# e.g., regex("Response|^meta-", ignore_case = TRUE)


# Make reg_class a factor for consistent ordering
combined_df_sub$reg_class <- factor(
  combined_df_sub$reg_class,
  levels = c("cis", "trans", "reg", "other")
)

# Order panels in reverse order of rho
combined_df_sub <- combined_df_sub %>%
  arrange(rho.rho) %>%
  mutate(clean_panel = factor(clean_panel, levels = unique(clean_panel)))

output_dir <- "02-single-cell-data-analysis/out/"
filename <- paste0(output_dir, "/iSensors_spearman_sc-endo.csv")
write.csv(file = filename, combined_df_sub)



# Create the bar plot
barPlot <- ggplot(combined_df_sub, aes(
  x = rho.rho,
  y = clean_panel,
  fill = reg_class,
  alpha = significance
)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c(
    cis      = "#fdbf6f",
    trans    = "#ff7f00",
    reg = "#b15928",
    other   = "gray"
  )) +
  scale_alpha_manual(values = c(
    "significant" = 1,
    "not significant" = 0.3
  )) +
  scale_x_continuous(breaks = seq(-1, 1, by = 0.1)) +
  labs(
    x = "Correlation (Rho)",
    y = "iSensors Panel",
    fill = "iSensors type",
  ) +
  guides(
    fill = guide_legend(order = 1),
    alpha = guide_legend(order = 2)
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(angle = 0, hjust = 1),
    axis.title.y = element_text(margin = margin(r = 10))
  )

print(barPlot)
# Save your ggplot object


filename <- paste0(output_dir, "/iSensors_spearman_nometapanels_rho0.45.svg")

filename
ggsave(filename, plot = barPlot, width = 8, height = 10, dpi = 300)
