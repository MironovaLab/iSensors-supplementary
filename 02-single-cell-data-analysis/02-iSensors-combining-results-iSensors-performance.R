library(tidyverse)
library(readr)
library(stringr)
in_dir <- "02-single-cell-data-analysis/in"
getwd()


# ---------- 1) Robust reader ----------
library(readr)

read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  
  if (ext %in% c("txt", "tsv")) {
    # Try tab first
    df <- tryCatch(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
                   error = function(e) NULL)
    if (!is.null(df) && ncol(df) > 1) return(df)
    
    # Fallback: whitespace + quoted strings (your format)
    df <- tryCatch(
      read.table(
        file = path,
        header = TRUE,
        sep = "",
        quote = "\"",
        comment.char = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(df) && ncol(df) > 1) return(as.data.frame(df))
    
    stop("Could not parse txt/tsv file: ", path)
  }
  
  stop("Unsupported extension: ", ext, " for file: ", path)
}


# ---------- 2) Ensure iSensor column is correct ----------
ensure_isensor_col <- function(df) {
  nms <- names(df)
  
  # Case-insensitive exact match
  idx <- which(tolower(nms) == "isensor")
  if (length(idx) == 1) {
    names(df)[idx] <- "iSensor"
    return(df)
  }
  
  # If first col is a simple row index and second col looks like iSensor IDs, use col2
  if (ncol(df) >= 2) {
    c1 <- df[[1]]
    c2 <- df[[2]]
    
    is_rowindex <- is.numeric(c1) || all(str_detect(as.character(c1), "^\\d+$"))
    looks_isensor <- all(str_detect(as.character(c2)[!is.na(c2)][1:min(20, sum(!is.na(c2)))],
                                    "^(AT-aux-(cis|trans|cistrans)-|AT\\dG\\d{5})"))
    
    if (is_rowindex && looks_isensor) {
      df <- df %>% select(-1) %>% rename(iSensor = 1)
      return(df)
    }
  }
  
  # Otherwise: assume first column holds iSensor
  df %>% rename(iSensor = 1)
}

# ---------- 3) Polish iSensor names ----------
polish_isensor <- function(x) {
  x <- as.character(x)
  x <- gsub("^AT-aux-cistrans-", "", x)
  x <- gsub("^AT-aux-cis-", "", x)
  x <- gsub("^AT-aux-trans-", "", x)
  x
}

# ---------- 4) Clean one table ----------
clean_one <- function(df, file_tag) {
  df <- ensure_isensor_col(df)
  
  df <- df %>%
    mutate(
      iSensor_raw = iSensor,
      iSensor = polish_isensor(iSensor)
    )
  
  # If polishing collapses distinct rows, aggregate numeric columns by mean
  # (and keep first for non-numeric).
  df <- df %>%
    group_by(iSensor) %>%
    summarise(
      across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
      across(where(~!is.numeric(.x)), ~ dplyr::first(.x)),
      .groups = "drop"
    )
  
  # Avoid name collisions on merge by tagging columns (except key)
  df %>%
    rename_with(
      ~ ifelse(.x %in% c("iSensor", "iSensor_raw"),
               .x,
               paste0(.x, "__", file_tag)),
      -iSensor
    )
}

# ---------- 5) Load, clean, merge ----------


files <- list.files(in_dir, pattern = "\\.(csv|txt|tsv)$", full.names = TRUE)
stopifnot(length(files) >= 1)

tables <- lapply(files, function(f) {
  tag <- tools::file_path_sans_ext(basename(f))
  df  <- read_any_table(f)
  
  # quick diagnostic (optional): print file + ncol
  message(basename(f), " -> ", ncol(df), " columns")
  
  clean_one(df, file_tag = tag)
})

merged <- Reduce(function(x, y) full_join(x, y, by = "iSensor"), tables) %>%
  arrange(iSensor)


merged_f <- merged %>%
  filter(!is.na(`id__bulk-limma-results`))


write_csv(merged, "02-single-cell-data-analysis/out/statistics_merged_by_iSensor.csv")

#DotPlot

library(tidyr)
library(ggplot2)
library(scales)

cols_use <- c(
    "effect__Statistics-exo-scdata",
  "Rho1__PerformanseTest1_stat",
  "rho.rho__iSensors_spearman_sc-endo"

)

cols_use <- c(
  "effect__Statistics-exo-scdata",
  #    "effect__bulk-limma-results", 
  "Rho1__PerformanseTest1_stat",
  "rho.rho__iSensors_spearman_sc-endo"
  
)

# keep only columns that actually exist (prevents cryptic errors)
cols_use <- cols_use[cols_use %in% colnames(merged_f)]
stopifnot(length(cols_use) > 0)

df_long <- merged_f %>%
  select(iSensor, all_of(cols_use)) %>%
  pivot_longer(
    cols = all_of(cols_use),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(metric, levels = cols_use)
  )

# robust limits so one outlier doesn't dominate
vals <- df_long$value
vals <- vals[is.finite(vals)]
lim <- quantile(abs(vals), probs = 0.98, na.rm = TRUE)
lims <- c(-lim, lim)

p_dot <- ggplot(df_long, aes(x = metric, y = iSensor)) +
  geom_point(
    aes(color = value, size = abs(value)),
    alpha = 0.9
  ) +
  scale_color_gradient2(
    low = muted("#4575B4"),
    mid = "white",
    high = muted("#D73027"),
    midpoint = 0,
    limits = lims,
    oob = squish,
    name = "Effect / ρ"
  ) +
  scale_size(
    range = c(0.6, 4.2),
    name = "|value|"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 7),
    axis.ticks  = element_blank(),
    legend.position = "right",
    plot.margin = margin(4, 6, 4, 6)
  )

p_dot


# ArrowPlot ----

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

rho_cutoff <- 0.45



merged_f<- read.csv("02-single-cell-data-analysis/out/statistics_merged_by_iSensor_polished.csv")
merged_f <- merged_f %>%
  filter(
    rho.rho__iSensors_spearman_sc.endo > rho_cutoff
  )

metric_cols <- c(
  "rho.rho__iSensors_spearman_sc.endo",
  "Rho1__PerformanseTest1_stat",
  "estimate__Statistics.exo.scdata"
)

names(merged_f)

# Map each metric -> p-value column used for significance


pval_map <- c(
"rho.rho__iSensors_spearman_sc.endo" = "p_value__iSensors_spearman_sc.endo",
 "Rho1__PerformanseTest1_stat"        = "p.value__PerformanseTest1_stat",
 "estimate__Statistics.exo.scdata"      = "p_adj__Statistics.exo.scdata"
)


rho_cutoff <- 0.45
alpha <- 0.05



# keep only columns that exist (prevents errors)
metric_cols <- metric_cols[metric_cols %in% names(merged_f)]
pval_map <- pval_map[names(pval_map) %in% metric_cols]
pval_map <- pval_map[pval_map %in% names(merged_f)]
stopifnot(length(metric_cols) > 0, length(pval_map) == length(metric_cols))

# ---- long table ----
df_long <- bind_rows(lapply(metric_cols, function(m) {
  pcol <- unname(pval_map[m])
  merged_f %>%
    transmute(
      iSensor,
      metric = m,
      value  = suppressWarnings(as.numeric(.data[[m]])),
      pval   = suppressWarnings(as.numeric(.data[[pcol]]))
    )
})) %>%
  mutate(
    metric = factor(metric, levels = metric_cols),
    sig = !is.na(pval) & pval < alpha,
    arrow = case_when(
      is.na(value) ~ NA_character_,
      value > 0    ~ "\u2191",  # ↑
      value < 0    ~ "\u2193",  # ↓
      TRUE         ~ "\u2192"   # →
    ),
    mag = abs(value),
    class = case_when(
      !sig ~ "n.s.",
      value > 0 ~ "up",
      value < 0 ~ "down",
      TRUE ~ "zero"
    )
  )

# size cap (avoid a few huge effects dominating)
mag_lim <- quantile(df_long$mag[is.finite(df_long$mag)], 0.98, na.rm = TRUE)

# ---- plot ----
p_arrow_sig <- ggplot(df_long, aes(x = metric, y = iSensor)) +
  geom_text(
    aes(
      label = arrow,
      size  = pmin(mag, mag_lim),
      color = class
    ),
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "up"   = "#D73027",  # strong red
      "down" = "#4575B4",  # strong blue

      "n.s." = "grey75"
    )
  ) +
  scale_size(range = c(2.5, 7.5), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.ticks  = element_blank(),
    legend.position = "none"   # <<< remove legend
  )

p_arrow_sig

file_path <- file.path("02-single-cell-data-analysis/out/statistics_merged_by_iSensor.svg")

ggsave(file_path, plot = p_arrow_sig,
       width = 2.5, height = 3, dpi = 300)
