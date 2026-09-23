# Brassinosteroid (BR) iSensor objects

Scripts that build the BR gene panels and the scored single-cell object used by
Figure 8 and Extended Data Figure 8.

The scored object (`out/GSE212230_iSensors_obj.rds`, ~322 MB) is **not stored in
this repository** — it exceeds GitHub's file size limit. Regenerate it with the
scripts below, or request it from the corresponding author
(victoria.mironova@ru.nl).

---

## Pipeline

Run from the **repository root**, not from this directory.

### 1. Build the BR gene panels

```r
source("00-BR-objects/00-BR-trans-panels.R")   # 7 trans-panels
source("00-BR-objects/01-BR-reg-panels.R")     # 4 reg-panels
```

Writes 11 panel definitions to `00-BR-objects/iSensors/*.rda`. These are small
and **are** committed, so this step can be skipped unless you want to rebuild
the panels from the source annotation.

| Panel type | Panels |
|---|---|
| `trans` (7) | Biosynthesis, Homeostasis, PositiveSignaling, NegativeSignaling, TF, TF-induced, TF-repressed |
| `reg` (4) | 2hr-up, 2hr-down, 4hr-up, 4hr-down |

`TF-induced` and `TF-repressed` are BL-response subsets of the `TF` panel,
derived from the GSE147589 bulk time course (Clark et al. 2021).

### 2. Score the single-cell data

```r
source("00-BR-objects/04-GSE212230-sc-analysis.R")
```

| | |
|---|---|
| **Input** | `00-iSensors-objects/data/GSE212230_BR_time_course_inner.rds.gz` |
| **Output** | `00-BR-objects/out/GSE212230_iSensors_obj.rds` (~322 MB, gitignored) |
| **Also writes** | `out/GSE212230_line_plots.pdf`, `out/GSE212230_zone_heatmap.pdf`, `out/GSE212230_celltype_heatmap.pdf` |

Computes all 11 BR panels plus the Auxin trans-panels on the 79,982-cell
GSE212230 dataset (Nolan et al. 2023), covering BRZ, mock, and a 0.5–8 h
brassinolide time course.

### 3. Generate the figures

```r
source("Manuscript-Figures/Figure8.R")
source("Manuscript-Figures/SupplementaryFigure8_BR_duration_celltype_layouts.R")
```

Both read `00-BR-objects/out/GSE212230_iSensors_obj.rds` from step 2.
`Figure8.R` additionally reads the root layout template from
`RealisticLayouts/out/new_ggPlantmap_epidermis.csv`.

---

## Source data and a provenance caveat

The starting point, `GSE212230_BR_time_course_inner.rds.gz`, is a **processed
Seurat object** derived from GEO accession
[GSE212230](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212230)
(Nolan et al. 2023).

**No script in this repository reconstructs that object from the raw GEO
submission** — the pipeline here begins from the processed object. Anyone
wishing to rebuild it from raw counts must repeat the original authors'
preprocessing as described in Nolan et al. 2023. The object is ~33 GB and is
therefore gitignored; it is available from the corresponding author on request.

This is a known limitation of the current deposit rather than an oversight in
the scripts, and is noted here so that users are not left searching for a step
that does not exist.
