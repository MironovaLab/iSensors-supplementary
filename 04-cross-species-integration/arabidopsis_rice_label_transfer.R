#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

options(Seurat.object.assay.version = "v5")
set.seed(100)

# Run from iSensors-supplementary/. Each phase is deliberately restartable so
# the two original multi-GB objects never need to coexist in memory.
defaults <- list(
  phase = "pilot",
  arabidopsis = "00-iSensors-objects/data/02_ra_integrated_Shahan_20250114.rds",
  rice = "00-iSensors-objects/data/Wang_rice.Rds",
  outdir = "00-iSensors-objects/data/arabidopsis_rice_transfer",
  orthologs = "00-iSensors-objects/data/arabidopsis_rice_transfer/arabidopsis_rice_one2one.tsv",
  cells_per_label = 500L,
  pilot_cells = 20000L,
  chunk_size = 25000L,
  nfeatures = 3000L,
  dims = 30L,
  min_score = 0.50,
  min_margin = 0.10,
  seed = 100L
)

parse_args <- function(defaults) {
  x <- commandArgs(trailingOnly = TRUE)
  out <- defaults
  i <- 1L
  while (i <= length(x)) {
    key <- sub("^--", "", x[[i]])
    if (!startsWith(x[[i]], "--")) stop("Unexpected argument: ", x[[i]])
    if (i == length(x)) stop("Missing value for --", key)
    key <- gsub("-", "_", key)
    if (!key %in% names(out)) stop("Unknown option: --", key)
    template <- out[[key]]
    value <- x[[i + 1L]]
    if (is.integer(template)) value <- as.integer(value)
    if (is.numeric(template) && !is.integer(template)) value <- as.numeric(value)
    out[[key]] <- value
    i <- i + 2L
  }
  out
}

cfg <- parse_args(defaults)
set.seed(cfg$seed)
dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  orthologs = cfg$orthologs,
  reference = file.path(cfg$outdir, "shahan_balanced_ortholog_reference.rds"),
  query = file.path(cfg$outdir, "wang_rice_root_ortholog_query.rds"),
  pilot_cells = file.path(cfg$outdir, "pilot_cells.txt"),
  pilot_predictions = file.path(cfg$outdir, "pilot_predictions.tsv"),
  pilot_confusion = file.path(cfg$outdir, "pilot_existing_vs_transferred.tsv"),
  pilot_label_summary = file.path(cfg$outdir, "pilot_label_summary.tsv"),
  full_predictions = file.path(cfg$outdir, "wang_rice_root_predictions.tsv"),
  final_query = file.path(cfg$outdir, "wang_rice_root_label_transferred.rds")
)

log_step <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}

strip_version <- function(x) sub("\\.[0-9]+$", "", as.character(x))

rna_counts <- function(obj) {
  assay <- obj[["RNA"]]
  layers <- Layers(assay)
  count_layers <- grep("^counts($|\\.)", layers, value = TRUE)
  if (!length(count_layers)) {
    stop("RNA assay has no counts layer. Available layers: ",
         paste(layers, collapse = ", "))
  }
  if (length(count_layers) > 1L) {
    log_step("Joining ", length(count_layers), " RNA count layers")
    assay <- JoinLayers(assay, layers = count_layers, new = "counts")
  }
  LayerData(assay, layer = "counts")
}

read_orthologs <- function(path) {
  if (!file.exists(path)) stop("Ortholog table not found: ", path)
  tab <- fread(path)
  required <- c("arabidopsis_gene", "rice_gene")
  if (!all(required %in% names(tab))) {
    stop("Ortholog table requires columns: ", paste(required, collapse = ", "))
  }
  tab <- unique(tab[, ..required])
  tab[, arabidopsis_gene := strip_version(arabidopsis_gene)]
  tab[, rice_gene := strip_version(rice_gene)]
  tab <- tab[nzchar(arabidopsis_gene) & nzchar(rice_gene)]
  # Defensively enforce one-to-one mappings even for user-supplied files.
  tab <- tab[!duplicated(arabidopsis_gene) &
             !duplicated(arabidopsis_gene, fromLast = TRUE)]
  tab <- tab[!duplicated(rice_gene) &
             !duplicated(rice_gene, fromLast = TRUE)]
  if (nrow(tab) < 1000L) {
    stop("Only ", nrow(tab), " one-to-one orthologs remain; check gene IDs.")
  }
  tab
}

make_orthologs <- function(path) {
  if (!requireNamespace("biomaRt", quietly = TRUE)) {
    stop("Package 'biomaRt' is required for --phase orthologs. ",
         "See 04-cross-species-integration/README.md.")
  }
  biomart_cache <- file.path(cfg$outdir, "biomart_cache")
  dir.create(biomart_cache, recursive = TRUE, showWarnings = FALSE)
  Sys.setenv(BIOMART_CACHE = normalizePath(biomart_cache, winslash = "/"))
  log_step("Querying Ensembl Plants for Arabidopsis-rice homology")
  plants <- biomaRt::useEnsemblGenomes(biomart = "plants_mart")
  ath <- biomaRt::useDataset("athaliana_eg_gene", mart = plants)
  # Query homology attributes from one dataset. This is more reliable than
  # getLDS(), which frequently returns HTTP 500 from Ensembl Plants.
  raw <- biomaRt::getBM(
    attributes = c(
      "ensembl_gene_id",
      "osativa_eg_homolog_ensembl_gene",
      "osativa_eg_homolog_orthology_type",
      "osativa_eg_homolog_orthology_confidence"
    ),
    mart = ath,
    uniqueRows = TRUE
  )
  setDT(raw)
  setnames(raw, c(
    "arabidopsis_gene", "rice_gene", "orthology_type",
    "orthology_confidence"
  ))
  raw <- raw[
    orthology_type == "ortholog_one2one" &
      orthology_confidence == 1
  ]
  raw[, arabidopsis_gene := strip_version(arabidopsis_gene)]
  raw[, rice_gene := strip_version(rice_gene)]
  one <- raw[!duplicated(arabidopsis_gene) &
             !duplicated(arabidopsis_gene, fromLast = TRUE)]
  one <- one[!duplicated(rice_gene) &
             !duplicated(rice_gene, fromLast = TRUE)]
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(one, path, sep = "\t")
  log_step("Saved ", nrow(one), " one-to-one orthologs to ", path)
}

atlas_to_broad <- function(x) {
  z <- tolower(as.character(x))
  out <- rep("Other", length(z))
  out[grepl("trichoblast|atrichoblast|epiderm", z)] <- "Epidermis"
  out[grepl("cortex", z)] <- "Cortex"
  out[grepl("endoderm", z)] <- "Endodermis"
  out[grepl("pericycle|ppp|xpp", z)] <- "Pericycle"
  out[grepl("xylem|phloem|procamb|vascular|pcc|pse", z)] <- "Vascular"
  out[grepl("columella|root.?cap|lrc|cscd|dying lrc|young lrc", z)] <- "Root_cap"
  out[grepl("(^|[^a-z])qc([^a-z]|$)|(^|[^a-z])csc([^a-z]|$)|initial|cei|lrp|(^|[^a-z])si([^a-z]|$)|meristem", z)] <- "Meristem_initials"
  out
}

rice_to_broad <- function(x) {
  z <- tolower(as.character(x))
  out <- rep("Other", length(z))
  out[grepl("epiderm|root.?hair", z)] <- "Epidermis"
  out[grepl("cortex|sclerenchyma|exoderm", z)] <- "Cortex"
  out[grepl("endoderm", z)] <- "Endodermis"
  out[grepl("pericycle", z)] <- "Pericycle"
  out[grepl("vascular|xylem|phloem|procamb", z)] <- "Vascular"
  out[grepl("root.?cap|columella|\\blrc\\b", z)] <- "Root_cap"
  out[grepl("meristem|initial|\\bqc\\b", z)] <- "Meristem_initials"
  out
}

choose_existing_annotation <- function(meta) {
  candidates <- c("cluster_names", "tissue_cluster", "tissue_cluster_names")
  hit <- candidates[candidates %in% colnames(meta)]
  if (!length(hit)) return(rep(NA_character_, nrow(meta)))
  as.character(meta[[hit[[1L]]]])
}

make_reference <- function() {
  orth <- read_orthologs(paths$orthologs)
  log_step("Loading full Shahan root object: ", cfg$arabidopsis)
  obj <- readRDS(cfg$arabidopsis)
  if (!inherits(obj, "Seurat")) stop("Arabidopsis input is not a Seurat object")
  if (!"Atlas_new" %in% colnames(obj[[]])) stop("Missing metadata: Atlas_new")

  labels <- as.character(obj$Atlas_new)
  eligible <- which(!is.na(labels) & nzchar(labels))
  split_cells <- split(colnames(obj)[eligible], labels[eligible])
  sampled <- unlist(lapply(split_cells, function(cells) {
    if (length(cells) <= cfg$cells_per_label) cells
    else sample(cells, cfg$cells_per_label)
  }), use.names = FALSE)
  balance <- data.table(
    Atlas_new = names(split_cells),
    available = lengths(split_cells),
    sampled = pmin(lengths(split_cells), cfg$cells_per_label)
  )
  fwrite(balance, file.path(cfg$outdir, "reference_sampling.tsv"), sep = "\t")
  log_step("Sampled ", length(sampled), " cells across ",
           length(split_cells), " Atlas_new labels")

  counts <- rna_counts(obj)
  ath_ids <- strip_version(rownames(counts))
  idx <- match(orth$arabidopsis_gene, ath_ids)
  keep <- !is.na(idx)
  if (sum(keep) < 1000L) stop("Too few Arabidopsis ortholog genes matched")
  mat <- counts[idx[keep], sampled, drop = FALSE]
  rownames(mat) <- orth$arabidopsis_gene[keep]
  meta <- obj[[]][sampled, , drop = FALSE]
  rm(counts, obj); invisible(gc())

  ref <- CreateSeuratObject(mat, assay = "RNA", meta.data = meta,
                            project = "Shahan_full_root_balanced")
  ref$Atlas_broad <- atlas_to_broad(ref$Atlas_new)
  ref <- NormalizeData(ref, verbose = FALSE)
  ref <- FindVariableFeatures(ref, selection.method = "vst",
                              nfeatures = min(cfg$nfeatures, nrow(ref)),
                              verbose = FALSE)
  ref <- ScaleData(ref, features = VariableFeatures(ref), verbose = FALSE)
  ref <- RunPCA(ref, features = VariableFeatures(ref),
                npcs = cfg$dims, verbose = FALSE)
  saveRDS(ref, paths$reference, compress = FALSE)
  log_step("Saved balanced reference: ", paths$reference)
}

make_query <- function() {
  orth <- read_orthologs(paths$orthologs)
  log_step("Loading Wang rice object: ", cfg$rice)
  obj <- readRDS(cfg$rice)
  if (!inherits(obj, "Seurat")) stop("Rice input is not a Seurat object")
  if (!"tissue" %in% colnames(obj[[]])) stop("Missing rice metadata: tissue")
  root_cells <- rownames(obj[[]])[tolower(as.character(obj$tissue)) == "root"]
  if (!length(root_cells)) stop("No cells matched tissue == 'Root'")
  log_step("Keeping ", length(root_cells), " rice root cells")

  counts <- rna_counts(obj)
  rice_ids <- strip_version(rownames(counts))
  idx <- match(orth$rice_gene, rice_ids)
  keep <- !is.na(idx)
  if (sum(keep) < 1000L) stop("Too few rice ortholog genes matched")
  mat <- counts[idx[keep], root_cells, drop = FALSE]
  # Arabidopsis IDs become common feature names in both objects.
  rownames(mat) <- orth$arabidopsis_gene[keep]
  meta <- obj[[]][root_cells, , drop = FALSE]
  rm(counts, obj); invisible(gc())

  query <- CreateSeuratObject(mat, assay = "RNA", meta.data = meta,
                             project = "Wang_rice_root")
  query$rice_annotation_original <- choose_existing_annotation(query[[]])
  query$rice_broad_original <- rice_to_broad(query$rice_annotation_original)
  saveRDS(query, paths$query, compress = FALSE)
  log_step("Saved ortholog-only rice root query: ", paths$query)
}

score_predictions <- function(score_assay, cells, label_key) {
  scores <- GetAssayData(score_assay, layer = "data")
  scores <- scores[, cells, drop = FALSE]
  score_ids <- sub("^prediction\\.score\\.", "", rownames(scores))
  keep <- score_ids %in% names(label_key)
  scores <- scores[keep, , drop = FALSE]
  score_ids <- score_ids[keep]
  if (nrow(scores) < 2L) stop("Prediction assay contains fewer than two labels")
  ord <- apply(scores, 2L, order, decreasing = TRUE)
  if (is.null(dim(ord))) ord <- matrix(ord, ncol = 1L)
  top1_i <- ord[1L, ]
  top2_i <- ord[2L, ]
  decoded <- unname(label_key[score_ids])
  if (anyNA(decoded)) {
    stop("Could not decode transferred labels: ",
         paste(head(rownames(scores)[is.na(decoded)]), collapse = ", "))
  }
  data.table(
    cell = cells,
    predicted.Atlas_new_raw = decoded[top1_i],
    prediction.score.max = scores[cbind(top1_i, seq_along(cells))],
    prediction.score.second = scores[cbind(top2_i, seq_along(cells))]
  )[, prediction.score.margin :=
       prediction.score.max - prediction.score.second]
}

transfer_one_lineage <- function(ref, query, cells, rice_broad) {
  q <- subset(query, cells = cells)
  q <- NormalizeData(q, verbose = FALSE)
  ref_broad <- atlas_to_broad(ref$Atlas_new)
  allowed <- rice_broad
  if (rice_broad == "Vascular") allowed <- c("Vascular", "Pericycle")
  state_pattern <- paste(
    c("^G1/S$", "^G2/M$", "^S$", "^Endocycle$", "^Unknown$",
      "^Stressed", "^C1$"),
    collapse = "|"
  )
  ref_cells <- colnames(ref)[
    ref_broad %in% allowed &
      !grepl(state_pattern, as.character(ref$Atlas_new), ignore.case = TRUE)
  ]
  if (length(ref_cells) < 100L) {
    stop("Only ", length(ref_cells), " reference cells for rice broad class ",
         rice_broad)
  }
  r <- subset(ref, cells = ref_cells)
  features <- intersect(VariableFeatures(r), rownames(q))
  if (length(features) < 500L) stop("Only ", length(features),
                                    " transfer features available")
  anchors <- FindTransferAnchors(
    reference = r,
    query = q,
    normalization.method = "LogNormalize",
    reference.reduction = "pca",
    reduction = "pcaproject",
    features = features,
    dims = seq_len(cfg$dims),
    k.anchor = 5,
    verbose = FALSE
  )
  original_labels <- as.character(r$Atlas_new)
  label_levels <- sort(unique(original_labels))
  safe_levels <- sprintf("AtlasLabel%03d", seq_along(label_levels))
  label_key <- setNames(label_levels, safe_levels)
  encoded_labels <- safe_levels[match(original_labels, label_levels)]
  anchor_pairs <- slot(anchors, "anchors")
  query_anchor_cells <- length(unique(anchor_pairs[, 2L]))
  if (query_anchor_cells < 3L) {
    stop("Only ", query_anchor_cells, " query anchor cells for ", rice_broad)
  }
  k_weight <- min(25L, query_anchor_cells - 1L)
  score_assay <- TransferData(
    anchorset = anchors,
    refdata = encoded_labels,
    dims = seq_len(cfg$dims),
    k.weight = k_weight,
    prediction.assay = TRUE,
    verbose = FALSE
  )
  pred <- score_predictions(score_assay, colnames(q), label_key)
  pred[, predicted.Atlas_broad := atlas_to_broad(predicted.Atlas_new_raw)]
  existing <- query[[]][pred$cell, , drop = FALSE]
  pred[, rice_annotation_original :=
         as.character(existing$rice_annotation_original)]
  pred[, rice_broad_original := as.character(existing$rice_broad_original)]
  pred[, broad_compatible :=
         rice_broad_original == predicted.Atlas_broad |
         (rice_broad_original == "Vascular" &
            predicted.Atlas_broad == "Pericycle")]
  pred[, predicted.Atlas_new := fifelse(
    prediction.score.max >= cfg$min_score &
      prediction.score.margin >= cfg$min_margin &
      broad_compatible,
    predicted.Atlas_new_raw,
    "Uncertain"
  )]
  pred[]
}

transfer_chunk <- function(ref, query, cells) {
  broad <- as.character(query$rice_broad_original)
  names(broad) <- colnames(query)
  groups <- split(cells, broad[cells])
  unsupported <- setdiff(names(groups), c(
    "Epidermis", "Cortex", "Endodermis", "Vascular", "Root_cap",
    "Meristem_initials"
  ))
  if (length(unsupported)) {
    stop("Unsupported rice broad classes: ", paste(unsupported, collapse = ", "))
  }
  rbindlist(lapply(names(groups), function(g) {
    log_step("  lineage ", g, ": ", length(groups[[g]]), " cells")
    transfer_one_lineage(ref, query, groups[[g]], g)
  }), use.names = TRUE)
}

run_transfer <- function(pilot = TRUE) {
  log_step("Loading compact reference and query")
  ref <- readRDS(paths$reference)
  query <- readRDS(paths$query)
  common <- intersect(rownames(ref), rownames(query))
  if (length(common) < 1000L) stop("Only ", length(common),
                                   " shared ortholog features")
  log_step(length(common), " shared one-to-one orthologs")

  all_cells <- colnames(query)
  if (pilot) {
    n <- min(cfg$pilot_cells, length(all_cells))
    # Stratify the pilot across existing rice labels where possible.
    strata <- split(all_cells, query$rice_annotation_original)
    quota <- max(1L, ceiling(n / length(strata)))
    cells <- unlist(lapply(strata, function(z) sample(z, min(length(z), quota))),
                    use.names = FALSE)
    if (length(cells) > n) cells <- sample(cells, n)
    if (length(cells) < n) {
      cells <- c(cells, sample(setdiff(all_cells, cells), n - length(cells)))
    }
    writeLines(cells, paths$pilot_cells)
    chunks <- split(cells, ceiling(seq_along(cells) / cfg$chunk_size))
    out_path <- paths$pilot_predictions
  } else {
    chunks <- split(all_cells, ceiling(seq_along(all_cells) / cfg$chunk_size))
    out_path <- paths$full_predictions
  }

  if (file.exists(out_path)) file.remove(out_path)
  for (i in seq_along(chunks)) {
    log_step("Transfer chunk ", i, "/", length(chunks), " (",
             length(chunks[[i]]), " cells)")
    pred <- transfer_chunk(ref, query, chunks[[i]])
    fwrite(pred, out_path, sep = "\t", append = i > 1L,
           col.names = i == 1L)
    rm(pred); invisible(gc())
  }
  log_step("Saved predictions: ", out_path)

  pred <- fread(out_path)
  confusion <- pred[, .N, by = .(
    rice_annotation_original, rice_broad_original,
    predicted.Atlas_broad, predicted.Atlas_new
  )][order(rice_annotation_original, -N)]
  label_summary <- pred[, .(
    cells = .N,
    median_score = as.numeric(median(prediction.score.max)),
    median_margin = as.numeric(median(prediction.score.margin)),
    accepted_fraction = mean(predicted.Atlas_new != "Uncertain"),
    broad_compatible_fraction = mean(broad_compatible)
  ), by = predicted.Atlas_new_raw][order(-cells)]

  if (pilot) {
    fwrite(confusion, paths$pilot_confusion, sep = "\t")
    fwrite(label_summary, paths$pilot_label_summary, sep = "\t")
    log_step("Pilot complete. Inspect ", paths$pilot_confusion)
  } else {
    meta <- as.data.frame(pred)
    rownames(meta) <- meta$cell
    meta$cell <- NULL
    query <- AddMetaData(query, metadata = meta)
    saveRDS(query, paths$final_query, compress = FALSE)
    log_step("Saved final annotated rice root object: ", paths$final_query)
  }
}

valid_phases <- c("orthologs", "reference", "query", "pilot", "full", "all")
if (!cfg$phase %in% valid_phases) {
  stop("--phase must be one of: ", paste(valid_phases, collapse = ", "))
}

if (cfg$phase %in% c("orthologs", "all")) make_orthologs(paths$orthologs)
if (cfg$phase %in% c("reference", "all")) make_reference()
if (cfg$phase %in% c("query", "all")) make_query()
if (cfg$phase %in% c("pilot", "all")) run_transfer(pilot = TRUE)
if (cfg$phase == "full") run_transfer(pilot = FALSE)
