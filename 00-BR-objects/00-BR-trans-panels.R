# 00-BR-trans-panels.R
# Creates BR iSensors trans-panels from gene annotation (hormones_gene_annotation.xlsx, BR sheet).
#
# ── Original panels (v1) ──────────────────────────────────────────────────────
#   ATH-br-trans-Homeostasis       - biosynthesis + metabolism + inactivation + transport (21)
#   ATH-br-trans-PositiveSignaling - positive components + receptors + general signalling (40)
#   ATH-br-trans-NegativeSignaling - negative components (12)
#   ATH-br-trans-TF                - all BR-regulated TFs (23)
#
# ── Alternative panels (v2, empirically informed) ─────────────────────────────
# Time-course analysis (02-BR-bulk-analysis.R / explore_trans_genes.R) revealed:
#   - Biosynthesis genes are STRONGLY REPRESSED by BL (negative feedback, mean sum_diff = -5.4)
#   - Metabolism + Inactivation genes are INDUCED by BL (homeostatic BR inactivation, +1.5)
#   - TF panel is bimodal: ~8 induced (PRE/BEE/AIF) vs ~7 repressed (AtIBH1/BIM/CESTA)
#   => Mixing opposing signals in one panel dilutes sensor power
#
# v2 panels split by direction of BL response:
#   ATH-br-trans-Biosynthesis    - 12 biosynthesis genes (REPRESSED by BL, inverse sensor)
#   ATH-br-trans-Homeostasis     - 9 genes: inactivation + metabolism + transport (INDUCED by BL)
#   ATH-br-trans-PositiveSignaling  (unchanged, 40 genes)
#   ATH-br-trans-NegativeSignaling  (unchanged, 12 genes)
#   ATH-br-trans-TF-induced      - 12 TF genes induced by BL (PRE6, PRE1, BEE1, AIF1, MYB30...)
#   ATH-br-trans-TF-repressed    - 11 TF genes repressed by BL (AtIBH1, ATBS1, BIM1, CESTA...)

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(iSensors)
})

setwd("D:/!GitHub/DigitalSensor-Toolbox/BRiSensors")

# ── 1. Load BR gene annotation ───────────────────────────────────────────────
br_genes <- read_excel("in/hormones_gene_annotation.xlsx", sheet = "BR") %>%
  as.data.frame()

cat("BR annotation:", nrow(br_genes), "genes\n")
cat("Types:\n"); print(sort(table(br_genes$Type), decreasing = TRUE))

# ── 2. Gene-set definitions ──────────────────────────────────────────────────

# v1 panels (unchanged grouping)
biosynthesis <- br_genes %>%
  filter(Type == "BR-biosynthesis") %>% pull(TAIR_ID)

# Homeostasis v2 = inactivation + metabolism + transport (biosynthesis removed)
homeostasis_v2 <- br_genes %>%
  filter(grepl("metabolism|inactivation|transport", Type, ignore.case = TRUE)) %>%
  pull(TAIR_ID)

negative <- br_genes %>%
  filter(grepl("negative component", Type, ignore.case = TRUE)) %>%
  pull(TAIR_ID)

tf_all <- br_genes %>%
  filter(grepl("TF|transcription factor", Type, ignore.case = TRUE)) %>%
  pull(TAIR_ID)

positive <- br_genes %>%
  filter(!grepl("biosynthesis|metabolism|inactivation|transport", Type, ignore.case = TRUE)) %>%
  filter(!grepl("negative component", Type, ignore.case = TRUE)) %>%
  filter(!grepl("TF|transcription factor", Type, ignore.case = TRUE)) %>%
  pull(TAIR_ID)

# TF panel split by BL response direction (from time-course analysis):
# Induced TFs: sum_diff > +0.5 across 6 timepoints
tf_induced_ids <- c(
  "AT1G26945",  # PRE6/KDR         sum_diff = +5.42
  "AT3G05800",  # AIF1/bHLH150     sum_diff = +4.23
  "AT5G39860",  # PRE1/BNQ1        sum_diff = +3.76
  "AT1G18400",  # BEE1/bHLH044     sum_diff = +3.76
  "AT3G28910",  # MYB30            sum_diff = +2.45
  "AT4G36540",  # BEE2/bHLH058     sum_diff = +1.95
  "AT3G06590",  # AIF2/bHLH148     sum_diff = +1.38
  "AT1G73830",  # BEE3/bHLH050     sum_diff = +0.70
  "AT3G28857",  # PRE5             sum_diff = +0.21
  "AT3G47710",  # PRE4/BNQ3        sum_diff = +0.39
  "AT3G48430",  # REF6             sum_diff = +0.16
  "AT1G32130"   # AtIWS1           sum_diff = +0.01 (low but in induced direction)
)
tf_induced <- intersect(tf_induced_ids, br_genes$TAIR_ID)

# Repressed TFs: sum_diff < -0.5 across 6 timepoints
tf_repressed_ids <- c(
  "AT2G43060",  # AtIBH1/bHLH158   sum_diff = -6.46
  "AT1G74500",  # ATBS1/PRE3       sum_diff = -1.84
  "AT5G08130",  # BIM1/bHLH046     sum_diff = -1.72
  "AT5G15160",  # PRE2/BNQ2        sum_diff = -0.87
  "AT1G25330",  # CESTA            sum_diff = -0.79
  "AT1G09250",  # AIF4/bHLH149     sum_diff = -0.74
  "AT3G17100",  # AIF3/bHLH147     sum_diff = -0.57
  "AT1G69010",  # BIM2/bHLH102     sum_diff = -0.14 (borderline)
  "AT5G04240",  # ELF6             sum_diff = -0.33
  "AT5G38860",  # BIM3/bHLH141     sum_diff = +0.06 (borderline, kept with repressed cluster)
  "AT1G67260"   # TCP1             sum_diff = 0 (weak)
)
tf_repressed <- intersect(tf_repressed_ids, br_genes$TAIR_ID)

# ── 3. Report gene counts ────────────────────────────────────────────────────
cat("\n── Panel gene counts ────────────────────────────────────────────────────\n")
cat("v1 Homeostasis (all):      ", length(c(biosynthesis, homeostasis_v2)), "\n")
cat("\nv2 panels:\n")
cat("  Biosynthesis:            ", length(biosynthesis), "\n")
cat("  Homeostasis (v2):        ", length(homeostasis_v2),
    "  (inactivation + metabolism + transport)\n")
cat("  PositiveSignaling:       ", length(positive), "\n")
cat("  NegativeSignaling:       ", length(negative), "\n")
cat("  TF-induced:              ", length(tf_induced), "\n")
cat("  TF-repressed:            ", length(tf_repressed), "\n")
cat("  TF-all (v1):             ", length(tf_all),
    "  (induced + repressed + borderline)\n")

tf_unassigned <- setdiff(tf_all, c(tf_induced, tf_repressed))
if (length(tf_unassigned) > 0) {
  cat("  TF genes unassigned to induced/repressed:",
      paste(tf_unassigned, collapse = ", "), "\n")
}

cat("\nGenes in Biosynthesis:\n")
print(br_genes$GeneName[br_genes$TAIR_ID %in% biosynthesis])

cat("\nGenes in Homeostasis v2:\n")
print(br_genes$GeneName[br_genes$TAIR_ID %in% homeostasis_v2])

cat("\nGenes in TF-induced:\n")
print(br_genes$GeneName[br_genes$TAIR_ID %in% tf_induced])

cat("\nGenes in TF-repressed:\n")
print(br_genes$GeneName[br_genes$TAIR_ID %in% tf_repressed])

# ── 4. Create all panels ─────────────────────────────────────────────────────
panels_def <- list(
  # ─ v1 panels (updated Homeostasis, rest unchanged) ─
  list(
    name  = "ATH-br-trans-Homeostasis",
    genes = homeostasis_v2,
    desc  = paste(
      "BR homeostatic response: inactivation, catabolism and transport genes.",
      "These genes are INDUCED by BL (BAS1, UGT73C5-6, BEN1, ST4A, CYP72C1, PVA12).",
      "Score increases with BR activity."
    )
  ),
  list(
    name  = "ATH-br-trans-PositiveSignaling",
    genes = positive,
    desc  = "Positive components of BR signalling: receptors, BSK, BES1/BZR1 pathway, kinases"
  ),
  list(
    name  = "ATH-br-trans-NegativeSignaling",
    genes = negative,
    desc  = "Negative components of BR signalling: BIN2, BKI1, 14-3-3, AtSK kinases"
  ),
  list(
    name  = "ATH-br-trans-TF",
    genes = tf_all,
    desc  = "All BR-regulated transcription factors: BIM, BEE, PRE, AIF, CESTA and others"
  ),
  # ─ v2 additional panels ─
  list(
    name  = "ATH-br-trans-Biosynthesis",
    genes = biosynthesis,
    desc  = paste(
      "BR biosynthesis genes (REPRESSED by BL via negative feedback).",
      "BR6ox2/DWF4/CYP90D1/CPD/ROT3/DWF1 etc.",
      "Score DECREASES with BR activity (inverse sensor)."
    )
  ),
  list(
    name  = "ATH-br-trans-TF-induced",
    genes = tf_induced,
    desc  = paste(
      "BR-induced transcription factors (positive BL response in bulk RNA-Seq).",
      "PRE6, PRE1, AIF1, BEE1, MYB30, BEE2, AIF2, BEE3, PRE5, PRE4.",
      "Score increases with BR activity."
    )
  ),
  list(
    name  = "ATH-br-trans-TF-repressed",
    genes = tf_repressed,
    desc  = paste(
      "BR-repressed transcription factors (negative BL response in bulk RNA-Seq).",
      "AtIBH1, ATBS1/PRE3, BIM1, PRE2, CESTA, AIF4, AIF3, BIM2, ELF6.",
      "Score DECREASES with BR activity (inverse sensor)."
    )
  )
)

for (p in panels_def) {
  cat("\nCreating:", p$name, "(", length(p$genes), "genes)\n")
  iSensorsTransPanelCreate(
    panel_name        = p$name,
    species           = "Arabidopsis Thaliana",
    gene_list         = p$genes,
    panel_description = p$desc
  )
}

cat("\nAll panels created and saved to ./iSensors/\n")

# ── 5. Verify ────────────────────────────────────────────────────────────────
cat("\nVerifying installation...\n")
BRpanels <- LoadSensors(setName = "BR", defaultPanels = FALSE, customPanels = TRUE)
cat("Loaded", length(BRpanels$panels), "panels:\n")
print(sort(names(BRpanels$panels)))
