## ---- Packages ----
install.packages(c("car"))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
library(car)
library(tidyverse)


## ---- Load ARF expression file ----
rice <- read.csv(
  "03-single-cell-rice/in/Rice_iSensors_mean_normed_290126.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

## Drop the unnamed index column

str(rice)

polish_isensor <- function(x) {
  x <- as.character(x)
  x <- gsub("^Os-aux-cistrans-", "", x)
  x <- gsub("^OS-aux-cis-", "", x)
  x <- gsub("^OS-aux-trans-", "", x)
  x
}

df <- rice %>%
  mutate(iSensor = polish_isensor(iSensor))

library(tidyverse)

df_long <- df %>%
  pivot_longer(
    cols = -iSensor,
    names_to = "Tissue",
    values_to = "Expression"
  )

df_long <- df_long %>%
  separate(
    Tissue,
    into = c("Tissue_Cleaned", "Organ"),
    sep = "_(?=[^_]+$)",   # split on LAST underscore
    remove = TRUE
  )


df_wide <- df_long %>%
  filter(iSensor %in% c("ARF", "synthesis", "ConjugationDeconjugation", "PAT")) %>%
  pivot_wider(
    id_cols  = c(Organ, Tissue_Cleaned),
    names_from  = iSensor,
    values_from = Expression
  ) %>%
  drop_na()


str(df_wide)
colnames(df_wide)[colnames(df_wide) == "ConjugationDeconjugation"] <- "Conjugation"


stopifnot(all(c("ARF","synthesis","Conjugation","PAT") %in% colnames(df_wide)))
summary(df_wide)

m1 <- lm(ARF ~ synthesis + Conjugation + PAT, data = df_wide)
summary(m1)
m2 <- lm(ARF ~ 0 + synthesis + Conjugation + PAT, data = df_wide)
summary(m2)


#Visualizing trend plot ----

df_wide$ARF_hat_m2 <- predict(m2)
rmse <- sqrt(mean(residuals(m2)^2))

# Build a ribbon around the identity line across the observed range
band <- tibble(
  x = seq(min(df_wide$ARF, na.rm = TRUE),
          max(df_wide$ARF, na.rm = TRUE),
          length.out = 200)
) %>%
  mutate(
    y  = x,
    lo = x - rmse,
    hi = x + rmse
  )

p<- ggplot(df_wide, aes(x = ARF, y = ARF_hat_m2)) +
  geom_ribbon(
    data = band,
    aes(x = x, ymin = lo, ymax = hi),
    inherit.aes = FALSE,
    fill = "grey80",
    alpha = 0.4
  ) +
  geom_point(size = 2.2, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 1) +
  coord_equal() +
  coord_cartesian(ylim = c(0, 5), xlim = c(0, 5)) +
  labs(
    x = "Observed ARF iSensor",
    y = "Predicted Auxin",
    title = "Auxin(ARF) = β1*Synthesis +\n +β2*Conjugation + β3*PAT"
  ) +
  theme_classic(base_size = 12)  +
  annotate(
    "text",
    x = 1, y = 4.5,
    label = "Adjusted R^2 = 0.89\np-value < 2.2e-16",
    size = 4
  )   +
  theme(plot.title = element_text(face = "bold"))
p
ggsave("03-single-cell-rice/out/ARF_mixedmodel_scatter_RSE.pdf", plot = p, width = 4, height = 4.0, dpi = 300)


#Vizualizing impact ----

coef_tbl <- summary(m2)$coefficients

coef_df <- tibble(
  Term     = rownames(coef_tbl),
  Estimate = coef_tbl[, "Estimate"],
  SE       = coef_tbl[, "Std. Error"],
  p_value  = coef_tbl[, "Pr(>|t|)"]
) %>%
  mutate(
    signif = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    # place label above positive bars and below negative bars
    label_y = if_else(
      Estimate >= 0,
      Estimate + SE + 0.05 * max(abs(Estimate + SE), na.rm = TRUE),
      Estimate - SE - 0.05 * max(abs(Estimate - SE), na.rm = TRUE)
    )
  )

p1 <- ggplot(coef_df, aes(x = Term, y = Estimate)) +
  geom_col(fill = "grey70") +
  geom_errorbar(
    aes(ymin = Estimate - SE, ymax = Estimate + SE),
    width = 0.2
  ) +
  geom_text(
    aes(y = label_y, label = signif),
    size = 5
  ) +
  ylab("Scaling coefficient") +
  xlab("") +
  theme_classic(base_size = 12)

p1

ggsave("03-single-cell-rice/out/ARF_mixedmodel_coeff.pdf", plot = p1, width = 3, height = 3.0, dpi = 300)




df_long_contrib <- df_wide |>
  mutate(
    Synthesis_c   = coef(m2)["synthesis"]   * synthesis,
    Conjugation_c = coef(m2)["Conjugation"] * Conjugation,
    PAT_c         = coef(m2)["PAT"]         * PAT
  ) |>
  select(Organ, Tissue_Cleaned, Synthesis_c, Conjugation_c, PAT_c) |>
  tidyr::pivot_longer(
    cols = ends_with("_c"),
    names_to = "Component",
    values_to = "Contribution"
  )

p2 <- ggplot(df_long_contrib,
             aes(x = Tissue_Cleaned, y = Contribution, fill = Component)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic()+
  theme(axis.title.y = element_blank())

p2  
ggsave("03-single-cell-rice/out/ARF_mixedmodel_contributions_tissue.pdf", plot = p2, width = 6, height = 8.0, dpi = 300)



p3<- ggplot(df_long_contrib,
            aes(x = Organ, y = Contribution, fill = Component)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_classic()+
  theme(axis.title.y = element_blank())

p3
ggsave("03-single-cell-rice/out/ARF_mixedmodel_contributions_stage.pdf", plot = p3, width = 6, height = 8.0, dpi = 300)


#Per tissue and stage performance ----
df_wide$resid_m2   <- resid(m2)
organ_summary <- df_wide %>%
  group_by(Organ) %>%
  summarise(
    mean_resid   = mean(resid_m2),
    median_resid = median(resid_m2),
    sd_resid     = sd(resid_m2),
    n            = n()
  )



## A-intervals based on the residual standard error (RSE) of m2
rse <- summary(m2)$sigma

ggplot(organ_summary, aes(x = Organ, y = mean_resid)) +
  
  # A2 interval (extended acceptable)
  geom_rect(
    ymin = -1 * rse, ymax =  1 * rse,
    xmin = -Inf, xmax = Inf,
    fill = "grey90", alpha = 0.6
  ) +
  
  # A1 interval (core acceptable)
  geom_rect(
    ymin = -0.3 * rse, ymax =  0.3 * rse,
    xmin = -Inf, xmax = Inf,
    fill = "grey80", alpha = 0.8
  ) +
  
  # zero reference
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # points and uncertainty
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_resid - sd_resid / sqrt(n),
      ymax = mean_resid + sd_resid / sqrt(n)
    ),
    width = 0.2
  ) +
  
  coord_flip() +
  
  # Title only
  labs(
    title = "Organ-Level Residuals of the Auxin model",
    #    title = " ",
    
    y = "Mean ARF (observed − predicted)",
    x = NULL
  ) +
  
  # Annotation in top-right corner
  # annotate(
  #   "text",
  #   x = Inf, y = Inf,
  #   label = "Adjusted R^2 = 0.78\np-value < 2.2e-16",
  #   hjust = 1.05,
  #   vjust = 1.1,
  #   size = 4
  # ) +
  
theme_classic() +
  theme(
    plot.title = element_text(face = "bold")
  )


ggsave("03-single-cell-rice/out/ARF_mixedmodel_stage.pdf", last_plot(), width = 6, height = 8.0, dpi = 300)


# ggplot(df_wide,
#        aes(x = ARF_hat_m2, y = ARF, color = Organ)) +
#   geom_point(alpha = 0.6) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   facet_wrap(~ Organ) +
#   xlab("Predicted ARF") +
#   ylab("Observed ARF") +
#   theme_classic()
# 
# organ_means <- df_wide %>%
#   group_by(Organ) %>%
#   summarise(
#     Observed = mean(ARF),
#     Predicted = mean(ARF_hat_m2)
#   ) %>%
#   tidyr::pivot_longer(
#     cols = c(Observed, Predicted),
#     names_to = "Type",
#     values_to = "ARF"
#   )
# 
# ggplot(organ_means, aes(x = Organ, y = ARF, fill = Type)) +
#   geom_col(position = "dodge") +
#   coord_flip() +
#   ylab("Mean ARF expression") +
#   theme_classic()


#Same for tissues ----

df_wide$ARF_hat_m2 <- predict(m2)


library(dplyr)
library(stringr)

df_wide2 <- df_wide %>%
  mutate(
    Tissue_Broad = case_when(
      str_detect(Tissue_Cleaned, regex("epidermis|guard-cell|atrichoblast|trichoblast", ignore_case = TRUE)) ~ "Epidermis",
      str_detect(Tissue_Cleaned, regex("vascular-cylinder|procambium|xylem|phloem|stele", ignore_case = TRUE)) ~ "Vascular",
      str_detect(Tissue_Cleaned, regex("parenchyma", ignore_case = TRUE)) ~ "Parenchyma",
      str_detect(Tissue_Cleaned, regex("mesophyll", ignore_case = TRUE)) ~ "Mesophyll",
      str_detect(Tissue_Cleaned, regex("cortex", ignore_case = TRUE)) ~ "Cortex",
      str_detect(Tissue_Cleaned, regex("endodermis|exodermis", ignore_case = TRUE)) ~ "Endodermal",
      str_detect(Tissue_Cleaned, regex("columella|root-cap|lrc", ignore_case = TRUE)) ~ "Root cap / LRC",
      str_detect(Tissue_Cleaned, regex("meristem|initials|proliferating", ignore_case = TRUE)) ~ "Meristematic",
      str_detect(Tissue_Cleaned, regex("fiber|sclerenchyma", ignore_case = TRUE)) ~ "Fibers",
      str_detect(Tissue_Cleaned, regex("tapetum|pollen|ovule|ovary|stigmata|filament|lodicule|scutellum|plumule", ignore_case = TRUE)) ~ "Reproductive",
      str_detect(Tissue_Cleaned, regex("unknown", ignore_case = TRUE)) ~ "Unknown",
      TRUE ~ Tissue_Cleaned
    )
  )

tissue_summary <- df_wide2 %>%
  group_by(Tissue_Broad) %>%
  summarise(
    mean_resid   = mean(resid_m2, na.rm = TRUE),
    median_resid = median(resid_m2, na.rm = TRUE),
    sd_resid     = sd(resid_m2, na.rm = TRUE),
    n            = n(),
    .groups = "drop"
  ) %>%
  filter(n >= 3) %>%
  arrange(mean_resid)

tissue_summary


#same but with the intervals ----
rse <- summary(m2)$sigma

ggplot(tissue_summary, aes(x = Tissue_Broad, y = mean_resid)) +
  
  # A2 interval (extended acceptable)
  geom_rect(
    ymin = -1 * rse, ymax =  1 * rse,
    xmin = -Inf, xmax = Inf,
    fill = "grey90", alpha = 0.6
  ) +
  
  # A1 interval (core acceptable)
  geom_rect(
    ymin = -0.3 * rse, ymax =  0.3 * rse,
    xmin = -Inf, xmax = Inf,
    fill = "grey80", alpha = 0.8
  ) +
  
  # zero reference
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # points and uncertainty
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_resid - sd_resid / sqrt(n),
      ymax = mean_resid + sd_resid / sqrt(n)
    ),
    width = 0.2
  ) +
  
  coord_flip() +
  ylab("Mean ARF (observed − predicted)") +
  xlab("") +
  theme_classic()+
  
  # Title only
  labs(
    title = "Tissue-Level Residuals of the Auxin Model",
    y = "Mean ARF (observed − predicted)",
    x = NULL
  ) +
  # 
  # # Annotation in top-right corner
  # annotate(
  #   "text",
  #   x = Inf, y = Inf,
  #   label = "Adjusted R^2 = 0.78\np-value < 2.2e-16",
  #   hjust = 1.05,
  #   vjust = 1.1,
  #   size = 4
  # ) +
  
theme_classic() +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave("03-single-cell-rice/out/ARF_mixedmodel_tissue.pdf", last_plot(), width = 6, height = 6.0, dpi = 300)

ggsave("03-single-cell-rice/out/ARF_mixedmodel_tissue.svg", last_plot(), width = 6, height = 6.0, dpi = 300)


##no axes lables:
ggplot(tissue_summary, aes(x = Tissue_Broad, y = mean_resid)) +
  
  # A2 interval (extended acceptable)
  geom_rect(
    ymin = -1 * rse, ymax =  1 * rse,
    xmin = -Inf, xmax = Inf,
    fill = "grey90", alpha = 0.6
  ) +
  
  # A1 interval (core acceptable)
  geom_rect(
    ymin = -0.3 * rse, ymax =  0.3 * rse,
    xmin = -Inf, xmax = Inf,
    fill = "grey80", alpha = 0.8
  ) +
  
  # zero reference
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  # points and uncertainty
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = mean_resid - sd_resid / sqrt(n),
      ymax = mean_resid + sd_resid / sqrt(n)
    ),
    width = 0.2
  ) +
  
  coord_flip() +
  ylab("") +
  xlab("") +
  theme_classic()+
  theme(
    axis.text.y  = element_blank(),
    axis.text.x = element_blank(),
  )
  
ggsave("03-single-cell-rice/out/ARF_mixedmodel_tissue_noaxeslabels.svg", last_plot(), width = 6, height = 6.0, dpi = 300)
write.csv(df_wide2,
          file = "03-single-cell-rice/out/rice_residual_summary.csv")

# ggplot(df_wide, aes(x = ARF_hat_m2, y = ARF)) +
#   geom_point(alpha = 0.6) +
#   geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
#   facet_wrap(~ Tissue_Cleaned, scales = "free", ncol = 4) +
#   xlab("Predicted ARF (m2)") +
#   ylab("Observed ARF") +
#   theme_classic()


# tissue_means <- df_wide %>%
#   group_by(Tissue_Cleaned) %>%
#   summarise(
#     Observed  = mean(ARF),
#     Predicted = mean(ARF_hat_m2),
#     n = n(),
#     .groups = "drop"
#   ) %>%
#   filter(n >= 3) %>%
#   pivot_longer(
#     cols = c(Observed, Predicted),
#     names_to = "Type",
#     values_to = "ARF"
#   )
# 
# ggplot(tissue_means, aes(x = Tissue_Cleaned, y = ARF, fill = Type)) +
#   geom_col(position = "dodge") +
#   coord_flip() +
#   ylab("Mean ARF expression") +
#   xlab("") +
#   theme_classic()



