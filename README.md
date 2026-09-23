# iSensors — supplementary code and data

Scripts, input tables and figure outputs for:

> **Computational reporters turn plant single-cell atlases into quantitative hormone response maps**
> Zemlyanskaya E., Rybakov M., van der Wijst C., Dolgikh V., Savina M., Gidding F.,
> van der Wijk E., Pasternak T., Shaw R., Wiebe D., Xu J., Mironova V.
> *(manuscript submitted — DOI to be added)*

This repository reproduces every figure and supplementary table in the paper.
The iSensors method itself lives in a separate package.

| | |
|---|---|
| **iSensors R package** | https://github.com/MironovaLab/iSensors |
| **Tutorial** | https://mironovalab.github.io/isensors-tutorial.html |
| **ggRootCellAtlas** (visualisation) | https://github.com/MironovaLab/ggRootCellAtlas |
| **This repository** | https://github.com/MironovaLab/iSensors-supplementary |

---

## What iSensors does

Plant hormones are metabolites, not gene products, so they are invisible in
single-cell transcriptomes. iSensors infers hormone response from expression data
by aggregating curated gene panels into a per-cell score — a computational
analogue of a genetic reporter. Panels are validated against the criteria used to
judge genetic reporters: response to exogenous hormone, dose and duration
sensitivity, recovery of endogenous gradients, and performance above randomised
controls.

---

## Repository layout

```
00-iSensors-objects/    Scripts building the scored Seurat objects (Arabidopsis + rice)
00-BR-objects/          Brassinosteroid panels and GSE212230 scoring pipeline
08-bulk-data-obtain-and-analyze/
                        Bulk microarray preparation and modelling (Figure 3)
Manuscript-Figures/     Figure scripts (*.R), inputs (in/), outputs (out/)
RealisticLayouts/       Root and epidermis layout templates
Supplementary-Tables/   Supplementary Tables S1–S7
Tutorial-iSensors/      Standalone worked example with bundled data
FIGURES.md              Figure → script → output index
```

**[FIGURES.md](FIGURES.md) is the place to start** if you want to regenerate a
specific figure.

---

## Getting started

### Install

```r
install.packages("remotes")
remotes::install_github("MironovaLab/iSensors")
remotes::install_github("MironovaLab/ggRootCellAtlas")
```

### Try it without any large downloads

`Tutorial-iSensors/` ships a small bundled dataset (3.3 MB) and runs end to end:

```r
source("Tutorial-iSensors/iSensors-tutorial-install.R")
source("Tutorial-iSensors/iSensors-tutorial-calculation.R")
```

### Reproduce a figure

Always run **from the repository root**:

```r
setwd("/path/to/iSensors-supplementary")
source("Manuscript-Figures/Figure4.R")
```

Scripts use paths relative to the repository root. Running them from inside
`Manuscript-Figures/` will fail to find their inputs.

---

## Data

Figure outputs and all input tables are committed. The **scored Seurat objects are
not** — they range from hundreds of megabytes to tens of gigabytes and exceed
GitHub's limits. Two ways to obtain them:

1. **Rebuild** them from the public accessions listed below, using the scripts in
   `00-iSensors-objects/` and `00-BR-objects/`. This is the fully self-contained
   route — every source dataset is publicly available.
2. **Request the pre-computed objects** from the corresponding author
   (victoria.mironova@ru.nl). These are provided as a convenience to save
   recomputation; they are derived from the public datasets below and are the
   subject of a separate forthcoming publication.

### Source datasets

| Dataset | Accession | Used for |
|---|---|---|
| Arabidopsis root, mock vs 1 µM IAA (Martin-Arevalillo et al. 2025) | GSE241573 | Figure 2 |
| Arabidopsis root atlas (Shahan et al. 2022) | see publication | Figures 4, 5, 8 |
| Arabidopsis life-cycle atlas (Guo et al. 2025) | see publication | Figure 6 |
| Rice life-cycle atlas (Wang et al. 2025) | see publication | Figure 7 |
| BR single-cell time course (Nolan et al. 2023) | GSE212230 | Figure 8 |
| BR bulk time course (Clark et al. 2021) | GSE147589 | BR reg-panels |
| 20 auxin microarray datasets | listed in Supplementary Table S3 | Figure 3 |

---

## Software environment

Figures in the paper were produced with:

| | |
|---|---|
| R | 4.5.x |
| Seurat / SeuratObject | 5.4.0 / 5.3.0 |
| limma / edgeR / sva | 3.66.0 / 4.8.2 / 3.58.0 |
| affy | 1.88.0 |
| universalmotif | 1.28.0 |
| AUCell / UCell / pROC | 1.32.0 / 2.14.0 / 1.19.0.1 |
| ggplot2 / pheatmap / patchwork | 4.0.3 / 1.0.13 / 1.3.2 |
| iSensors / ggRootCellAtlas / ggPlantmap | 1.2.4 / 0.0.0.9000 / 1.1.0 |

Root reannotation (Figure 4) additionally used `slingshot` 2.0.0 and `CytoTRACE`
0.3.3 under R 4.1.

---

## Licence and citation

Released under the MIT licence. If you use this code or the
iSensors framework, please cite the paper (to be provided) and the
[iSensors package](https://github.com/MironovaLab/iSensors).

Questions and issues: https://github.com/MironovaLab/iSensors-supplementary/issues
