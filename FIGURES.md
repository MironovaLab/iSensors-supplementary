# Figure index

Maps every manuscript figure to the script that produces it, the data objects it
requires, and its output files.

Run all scripts **from the repository root**, not from inside `Manuscript-Figures/`:

```r
setwd("/path/to/iSensors-supplementary")
source("Manuscript-Figures/Figure4.R")
```

Objects under `00-iSensors-objects/data/` are too large for GitHub and are not
included. Build them with the scripts in `00-iSensors-objects/`, or download them
from the Zenodo deposit (DOI: *to be added*). See the README for details.

---

## Main figures

### Figure 1 — The iSensors framework
Schematic, assembled by hand in Illustrator. No generating script.

| | |
|---|---|
| **Output** | `Manuscript-Figures/out/2026-06-19_Fig1.pdf`, `2026-06-19_Fig1_p.svg` |

The `_p` suffix marks the polished (hand-edited) version.

### Figure 2 — iSensors infer auxin responses from single-cell data

| | |
|---|---|
| **Script** | `Manuscript-Figures/Figure2.R` |
| **Input** | `00-iSensors-objects/data/iSensors-Martin-Arevalillo-auxin-root2025_mean.rds` |
| **Also reads** | `Manuscript-Figures/in/Statistics-exo-scdata.csv` |
| **Output** | `out/Fig2A_AuxinTreated_DimPlot_integrated.pdf`, `out/Fig2A_AuxinTreated_ARF_FeaturePlot_integrated.pdf`, `out/Fig2B_AuxinTreated_ARF_FeaturePlot_split.pdf`, `out/Fig2CDE_iSensors_replica_heatmap_plus_forest.pdf` |

### Figure 3 — Evaluation on bulk transcriptomic data
Assembled by hand from the analyses in `08-bulk-data-obtain-and-analyze/`.

| | |
|---|---|
| **Analysis** | `08-bulk-data-obtain-and-analyze/01–03-*.ipynb` (data preparation, integration, iSensor scoring), then `04-bulk-linear-modelling.R`, `05-MA-ComplexHeatmap.R`, `06-UMAPs-diversity.R` |
| **Output** | `Manuscript-Figures/out/2026-06-24_Fig3.pdf`, `2026-06-24_Fig3_p.svg` |

The `.rds` intermediates produced by the notebooks are gitignored; regenerate them
by running the notebooks in order.

### Figure 4 — Validation against endogenous auxin gradients

| | |
|---|---|
| **Script** | `Manuscript-Figures/Figure4.R` |
| **Input** | `00-iSensors-objects/data/shahan-iSensors-obj-groundtruth.rds`, `shahan-roottip-iSensors-obj.rds` |
| **Output** | `out/Figure4B_dimplot*.{pdf,svg}`, `out/Figure4B_featureplot_ARF.*`, `out/Figure4C_spearman_2x2.*`, `out/Figure4C_*_spearman.svg`, `out/Figure4D_spearman_barplot.*`, `out/Figure4D_spearman_sc-endo.csv`, `out/Figure4E_*_root.*`, `out/Figure4F_performance_summary.*`, `out/Figure4H_ARF_root.*` |

### Figure 5 — Fine-scale gradients in the root epidermis

| | |
|---|---|
| **Script** | `Manuscript-Figures/Figure5.R` |
| **Input** | `00-iSensors-objects/data/shahan-roottip-iSensors-obj.rds`, `shahan_epidermis_iSensors_obj.rds` |
| **Output** | `out/Figure5A_ARF_barplot_by_celltype.*`, `out/Figure5A_ARF_root_layout.*`, `out/Figure5B_epidermis_DimPlot.*`, `out/Figure5B_*_FeaturePlot.*` |

Panels D–E (PIN2 immunolocalisation, AUX1::AUX1–YFP) are microscopy, not script output.

### Figure 6 — Systems-level map across Arabidopsis development

| | |
|---|---|
| **Script** | `Manuscript-Figures/Figure6.R` |
| **Input** | `00-iSensors-objects/data/Guo/guo_avg_exp_iSensors_mean.rds` |
| **Output** | `out/Figure6A_stacked_barplot_4sensors.*`, `out/Figure6C_ARF_scatter_observed_predicted.*`, `out/Figure6D_ARF_model_coefficients.*`, `out/Figure6E*`, `out/Figure6F*` |
| **Variant** | `Figure6_noroot.R` — same analysis excluding root tissues, writing `out/Figure6*_noroot_*` |

### Figure 7 — Cross-species comparison with rice

Built from several component scripts, then compiled. Run in this order:

| Step | Script | Produces |
|---|---|---|
| 1 | `Figure7.R` | `out/Figure7A_rice_ARF_scatter.*`, `out/Figure7B_rice_organ_residuals.*`, `out/Figure7C_rice_tissue_residuals.*` |
| 2 | `Figure7_root_species_stacked_barplot.R` | `out/Figure7D_root_Arabidopsis_vs_rice_stacked_barplot_4sensors.*`, `in/Figure7D_root_Arabidopsis_vs_rice_iSensor_means.csv` |
| 3 | `Figure7_root_ARF_mirrored.R` | `out/Figure7D_root_Arabidopsis_vs_rice_ARF_mirrored.*` |
| 4 | `Figure7_root_species_tissue_average.R` | `out/Figure7D_root_species_rank_difference_heatmap.*`, `out/Figure7D_root_epidermis_subtype_means.*` |
| 5 | `Figure7_CD_root_comparison_one_row.R` | `out/Figure7_CD_root_comparison_one_row.*` |
| 6 | `Figure7_compile_with_root_panels.R` | `out/Figure7_combined_with_root_comparison.{pdf,png}` — **the published figure** |

`Figure7_root_model_diagnostics.R` produces supporting model-diagnostic panels
(`out/Figure7_root_model_contributions.*`, `out/Figure7_root_model_coefficient_sensitivity.*`).

**Input:** `00-iSensors-objects/data/wang_rice_avg_exp_iSensors_mean.rds` and the
Shahan root object.

### Figure 8 — BR iSensors and BR–auxin cross-hormone analysis

| | |
|---|---|
| **Script** | `Manuscript-Figures/Figure8.R` |
| **Input** | `00-BR-objects/out/GSE212230_iSensors_obj.rds` (~322 MB, not in repo — see `00-BR-objects/README.md`) |
| **Also reads** | `RealisticLayouts/out/new_ggPlantmap_epidermis.csv` |
| **Output** | `out/Figure8.{pdf,png,svg}` — **the published figure**; plus `out/Figure8_concentration_stats.csv`, `out/Figure8_duration_stats.csv` |
| **Panel B alone** | `Figure8_BR_concentration_duration.R` → `out/Figure8_BR_concentration_duration.*` |

---

## Supplementary / Extended Data figures

| Figure | Script | Output |
|---|---|---|
| Auxin-treated DimPlot, split | `SupplementaryFigures.R` | `out/Supplementary_AuxinTreated_DimPlot_splitted.pdf` |
| log2 AUX vs CTR, masked | `SupplementaryFigures.R` | `out/FigureS3_log2_AUX_vs_CTR_masked.pdf` |
| Guo atlas heatmap | `FigureS6_Guo_atlas_heatmap.R` | `out/FigureS6_Guo_atlas_heatmap.{pdf,svg}` |
| Gradient correlation | `FigureS_gradient_correlation.R` | `out/FigureS_C_gradient_spearman.*`, `out/FigureS_C_gradient_extended.*` |
| Scoring benchmark vs AUCell/UCell | `FigureS_scoring_comparison.R` **then** `FigureS_scoring_final.R` | `out/scoring_comparison/FigureS_scoring_final.{pdf,svg}` |
| BR panel validation | `SupplementaryFigure8_BR_duration_celltype_layouts.R` | `out/SupplementaryFigure8_BR_duration_celltype_layouts.{pdf,png,svg}` |

> **Order matters for the scoring benchmark.** `FigureS_scoring_comparison.R` computes
> the AUROC and Spearman tables (`auroc_all_panels.csv`, `auroc_selected_comparison.csv`,
> `gradient_spearman.csv`) and the per-cell score objects that `FigureS_scoring_final.R`
> then reads to render the published panel. Running `_final` alone will fail.
> The intermediate `.rds` files it writes are gitignored; the `.csv` outputs are committed.

---

## Supporting analyses

| Directory | Contents |
|---|---|
| `00-iSensors-objects/` | Scripts that build the scored Seurat objects for Arabidopsis (Shahan, Guo, Martin-Arevalillo) and rice (Wang), plus the Arabidopsis→rice orthology and label-transfer tables |
| `00-BR-objects/` | BR gene panels and the GSE212230 scoring pipeline — see its README |
| `08-bulk-data-obtain-and-analyze/` | Bulk microarray preparation, integration and modelling behind Figure 3 |
| `RealisticLayouts/` | Root and epidermis layout templates for `ggRootCellAtlas` |
| `Tutorial-iSensors/` | Standalone worked example with a small bundled dataset |
| `Supplementary-Tables/` | Supplementary Tables S1–S7 |
