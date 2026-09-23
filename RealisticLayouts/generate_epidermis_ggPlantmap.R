# Generate new_ggPlantmap_epidermis.csv from EpidermisT_AT_Carmen.svg
# Run with working directory = d:/!GitHub/DigitalSensor-Toolbox/
# Output: iSensors-supplementary/Manuscript-Figures/in/new_ggPlantmap_epidermis.csv
#
# SVG parsing adapted from CW_scripts_all/CW_AT_T_atlas.R (Carmen Vandenberghe)

library(xml2)
library(stringr)
library(ggPlantmap)

svg_path <- "RealisticLayouts/EpidermisT_AT_Carmen.svg"
out_path <- "RealisticLayouts/out/new_ggPlantmap_epidermis.csv"

# ── Bezier helper ─────────────────────────────────────────────────────────────
calculate_bezier_point <- function(t, P0, P1, P2, P3) {
  (1 - t)^3 * P0 + 3 * (1 - t)^2 * t * P1 + 3 * (1 - t) * t^2 * P2 + t^3 * P3
}

# ── Parse SVG ─────────────────────────────────────────────────────────────────
new.ggPlantmap <- data.frame(
  matrix(vector(), 0, 7,
         dimnames = list(c(), c("ROI.name","Cell_type","Colour","ROI.id","point","x","y"))),
  stringsAsFactors = FALSE
)

svg_file    <- readLines(svg_path, warn = FALSE)
svg_content <- paste(svg_file, collapse = "\n")
svg         <- read_xml(svg_content)
ns          <- c(svg = "http://www.w3.org/2000/svg")
paths       <- xml_find_all(svg, "//svg:path", ns)
paths       <- paths[grepl("[zZ]", xml_attr(paths, "d"))]

for (i in seq_along(paths)) xml_set_attr(paths[[i]], "id", paste0("poly", i))

style_node <- xml_find_first(svg, "//svg:style", ns)
css_text   <- gsub("[\r\n\t]", "", xml_text(style_node))
rules      <- str_split(css_text, "\\.st")[[1]]
rules      <- rules[nzchar(rules)]

css_list <- lapply(rules, function(x) {
  parts <- str_match(x, "(\\d+)\\{([^}]*)\\}")
  if (!is.na(parts[1])) return(setNames(parts[3], paste0("st", parts[2])))
  NULL
})
css_list <- unlist(css_list)

for (path_node in paths) {
  class_val <- xml_attr(path_node, "class")
  style_val <- css_list[class_val]
  if (!is.na(style_val)) xml_set_attr(path_node, "style", style_val)
}

for (path in paths) {
  fill_value <- xml_attr(path, "style")
  fill_value <- sub(".*fill:([^;]*);.*", "\\1", fill_value)
  if (is.na(fill_value) || fill_value == "#231f20") next

  path_id     <- xml_attr(path, "id")
  path_d      <- xml_attr(path, "d")
  coordinates <- str_extract_all(path_d, "[A-Za-z]|-?\\d*\\.?\\d+|-?\\d+\\.?\\d*")[[1]]

  x <- 0; y <- 0; l <- 1; k <- 1; point <- 0
  h <- 0; v <- 0; c_flag <- 0

  for (i in seq_along(coordinates)) {
    co <- coordinates[i]
    if      (co == "h") { h <- 1; l <- 0; v <- 0; c_flag <- 0; k <- 1
    } else if (co == "v") { h <- 0; l <- 0; v <- 1; c_flag <- 0; k <- 1
    } else if (co %in% c("m","l")) { h <- 0; l <- 1; v <- 0; c_flag <- 0; k <- 1
    } else if (co %in% c("M","L")) { h <- 0; l <- 2; v <- 0; c_flag <- 0; k <- 1
    } else if (co == "H") { h <- 2; l <- 0; v <- 0; c_flag <- 0; k <- 1
    } else if (co == "V") { h <- 0; l <- 0; v <- 2; c_flag <- 0; k <- 1
    } else if (co == "C") { h <- 0; l <- 0; v <- 0; c_flag <- 2; k <- 1
    } else if (co == "c") { h <- 0; l <- 0; v <- 0; c_flag <- 1; k <- 1
    } else if (!co %in% c("l","h","v")) {
      val <- suppressWarnings(as.numeric(co))
      if (is.na(val)) next

      if      (l == 1 && k == 1) { x <- x + val; y <- y + as.numeric(coordinates[i+1]); point <- point+1; new.ggPlantmap[nrow(new.ggPlantmap)+1,] <- c(path_id,NA,fill_value,path_id,point,x,y); k <- 2
      } else if (l == 2 && k == 1) { x <- val;      y <- as.numeric(coordinates[i+1]);         point <- point+1; new.ggPlantmap[nrow(new.ggPlantmap)+1,] <- c(path_id,NA,fill_value,path_id,point,x,y); k <- 2
      } else if (k > 1)  { k <- k - 1
      } else if (h == 1) { x <- x + val; point <- point+1; new.ggPlantmap[nrow(new.ggPlantmap)+1,] <- c(path_id,NA,fill_value,path_id,point,x,y)
      } else if (h == 2) { x <- val;     point <- point+1; new.ggPlantmap[nrow(new.ggPlantmap)+1,] <- c(path_id,NA,fill_value,path_id,point,x,y)
      } else if (v == 1) { y <- y + val; point <- point+1; new.ggPlantmap[nrow(new.ggPlantmap)+1,] <- c(path_id,NA,fill_value,path_id,point,x,y)
      } else if (v == 2) { y <- val;     point <- point+1; new.ggPlantmap[nrow(new.ggPlantmap)+1,] <- c(path_id,NA,fill_value,path_id,point,x,y)
      } else if (c_flag == 2 && k == 1) {
        P0x <- x; P0y <- y
        P1x <- as.numeric(coordinates[i]);   P1y <- as.numeric(coordinates[i+1])
        P2x <- as.numeric(coordinates[i+2]); P2y <- as.numeric(coordinates[i+3])
        P3x <- as.numeric(coordinates[i+4]); P3y <- as.numeric(coordinates[i+5])
        for (t in seq(0, 1, length.out = 10)) {
          sx <- calculate_bezier_point(t, P0x, P1x, P2x, P3x)
          sy <- calculate_bezier_point(t, P0y, P1y, P2y, P3y)
          point <- point + 1
          new.ggPlantmap[nrow(new.ggPlantmap)+1,] <- c(path_id,NA,fill_value,path_id,point,sx,sy)
        }
        x <- P3x; y <- P3y; k <- 6
      } else if (c_flag == 1 && k == 1) {
        P0x <- x; P0y <- y
        P1x <- x + as.numeric(coordinates[i]);   P1y <- y + as.numeric(coordinates[i+1])
        P2x <- x + as.numeric(coordinates[i+2]); P2y <- y + as.numeric(coordinates[i+3])
        P3x <- x + as.numeric(coordinates[i+4]); P3y <- y + as.numeric(coordinates[i+5])
        for (t in seq(0, 1, length.out = 10)) {
          sx <- calculate_bezier_point(t, P0x, P1x, P2x, P3x)
          sy <- calculate_bezier_point(t, P0y, P1y, P2y, P3y)
          point <- point + 1
          new.ggPlantmap[nrow(new.ggPlantmap)+1,] <- c(path_id,NA,fill_value,path_id,point,sx,sy)
        }
        x <- P3x; y <- P3y; k <- 6
      }
    }
  }
}

# ── Type conversion + y-flip ──────────────────────────────────────────────────
new.ggPlantmap$ROI.name  <- as.character(new.ggPlantmap$ROI.name)
new.ggPlantmap$Cell_type <- as.character(new.ggPlantmap$Cell_type)
new.ggPlantmap$Colour    <- as.character(new.ggPlantmap$Colour)
new.ggPlantmap$ROI.id    <- as.character(new.ggPlantmap$ROI.id)
new.ggPlantmap$point     <- as.numeric(new.ggPlantmap$point)
new.ggPlantmap$x         <- as.numeric(new.ggPlantmap$x)
new.ggPlantmap$y         <- as.numeric(new.ggPlantmap$y)

new.ggPlantmap$y <- new.ggPlantmap$y + 2 * (
  mean(c(max(new.ggPlantmap$y, na.rm = TRUE), min(new.ggPlantmap$y, na.rm = TRUE))) -
    new.ggPlantmap$y
)

# ── Inspect colours before assigning cell types ───────────────────────────────
message("Unique colours in SVG (check against assignments below):")
print(unique(new.ggPlantmap$Colour))

# ── Assign cell types by fill colour ─────────────────────────────────────────
# Colours from EpidermisT_AT_Carmen.svg as used in CW_AT_T_atlas.R
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#FFFFFF"] <- "Atrichoblast-m1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#DBD3D2"] <- "Atrichoblast-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#EBF5EA"] <- "Atrichoblast-t"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#D1E9D1"] <- "Atrichoblast-e1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#ADBCAF"] <- "Atrichoblast-e2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#868F89"] <- "Atrichoblast-d"

new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#CCCADB"] <- "Trichoblast-m1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#A2A1BF"] <- "Trichoblast-m2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#908FB3"] <- "Trichoblast-t"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#85859C"] <- "Trichoblast-e1"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#72727F"] <- "Trichoblast-e2"
new.ggPlantmap$Cell_type[new.ggPlantmap$Colour == "#56565C"] <- "Trichoblast-d"

# Check for unassigned polygons
unassigned <- unique(new.ggPlantmap$Colour[is.na(new.ggPlantmap$Cell_type)])
if (length(unassigned)) {
  message("WARNING: unassigned colours (add mappings above): ",
          paste(unassigned, collapse = ", "))
} else {
  message("All polygons assigned to cell types.")
}

# ── Quick preview ─────────────────────────────────────────────────────────────
ggPlantmap.plot(data = new.ggPlantmap, layer = Cell_type)

# ── Save ──────────────────────────────────────────────────────────────────────
write.csv(new.ggPlantmap, out_path, row.names = FALSE)
message("Saved: ", out_path)
getwd()
