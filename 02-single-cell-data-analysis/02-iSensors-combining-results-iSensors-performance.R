library(tidyverse)
library(readr)
library(stringr)
in_dir <- "02-single-cell-data-analysis/in"


# ---------- 1) Robust reader ----------
read_any_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  
  if (ext == "csv") {
    return(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  }
  
  if (ext %in% c("txt", "tsv")) {
    # Try common delims first, fall back to whitespace
    tries <- list(
      function() readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
      function() readr::read_delim(path, delim = ",", show_col_types = FALSE, progress = FALSE),
      function() readr::read_delim(path, delim = ";", show_col_types = FALSE, progress = FALSE),
      function() readr::read_table2(path, show_col_types = FALSE, progress = FALSE) # whitespace
    )
    
    for (f in tries) {
      df <- tryCatch(f(), error = function(e) NULL)
      if (is.null(df)) next
      
      # Accept if we got >1 column OR if first column contains separators (rare)
      if (ncol(df) > 1) return(df)
    }
    
    stop("Could not parse txt/tsv file with common delimiters: ", path)
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

dir.create("out", showWarnings = FALSE)
write_csv(merged, "out/merged_by_iSensor.csv")
