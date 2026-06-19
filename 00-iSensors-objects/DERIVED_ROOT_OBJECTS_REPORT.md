# Derived Shahan and Wang root objects

## Purpose

This report documents the smaller objects created after Arabidopsis-to-rice
label transfer. They were designed to support downstream visualization and
iSensor analysis without repeatedly loading the original multi-gigabyte
datasets.

The transfer methodology itself is documented in:

```text
../04-cross-species-integration/REANNOTATION_REPORT.md
```

## Balanced full-gene Shahan object

Generation script:

```text
create-balanced-shahan-full-genes.R
```

Output:

```text
data/shahan_balanced_full_root_full_genes.rds
```

Properties:

| Property | Value |
|---|---:|
| Cells | 23,494 |
| RNA genes | 39,578 |
| `Atlas_new` labels | 84 |
| Maximum cells per label | 500 |
| Random seed | 100 |

Unlike the orthologue-only transfer reference, this object preserves the full
Arabidopsis RNA gene set. Integrated assays, graphs, reductions and scaled
matrices were removed to reduce its size.

The exact sampling is recorded in:

```text
data/shahan_balanced_full_root_sampling.tsv
```

## All-gene confidently reannotated Wang object

Generation script:

```text
create-confident-wang-iSensors.R
```

Input:

```text
data/arabidopsis_rice_transfer/
wang_rice_root_label_transferred_all_genes_umap.rds
```

Output:

```text
data/wang_rice_root_confident_transferred_all_genes.rds
```

Only cells for which `predicted.Atlas_new != "Uncertain"` were retained.

Properties:

| Property | Value |
|---|---:|
| Cells | 17,201 |
| RNA genes | 37,960 |
| RNA layers | counts and normalized data |
| Reductions | PCA and `umap.allgenes` |

The object retains the original Wang metadata, transfer scores, raw and
accepted transferred labels, broad annotations and de novo cluster identity.

This is the preferred compact object when analyses require genes outside the
auxin panels.

## Auxin-panel Wang object

Output:

```text
data/wang_rice_root_confident_Auxin_OSA_iSensors.rds
```

Properties:

| Property | Value |
|---|---:|
| Cells | 17,201 |
| RNA genes retained | 216 |
| iSensor assays | RNA and `iSensors_mean` |
| Reductions | PCA and `umap.allgenes` |
| Auxin trans-panels | 7 |

Panel sizes:

| Panel | Genes |
|---|---:|
| ARF | 30 |
| IAA | 30 |
| PAT | 18 |
| Receptors | 58 |
| Synthesis | 19 |
| Transport | 41 |
| Local ConjugationDeconjugation | 20 |

All panel genes listed in the construction summary were present in the object:

```text
data/wang_rice_root_confident_Auxin_OSA_panel_genes.tsv
```

## Local rice ConjugationDeconjugation panel

The released OSA iSensors registry did not contain a
ConjugationDeconjugation panel at the time of this work. A local rice panel was
therefore generated from the installed Arabidopsis panel.

Source panel:

```text
ATH-aux-trans-ConjugationDeconjugation
```

Mapping:

```text
data/OSA-aux-trans-ConjugationDeconjugation_all_high_confidence.tsv
```

Selection included high-confidence:

- one-to-one orthologues;
- one-to-many orthologues;
- many-to-many orthologues.

The final local panel contained 20 rice genes representing 14 Arabidopsis
source genes. This broader mapping was selected because a one-to-one-only
panel contained too few genes.

The object stores panel provenance in:

```text
rice_panel@misc$local_OSA_ConjugationDeconjugation
```

Paralog expansion should be considered when interpreting the score because
multiple rice genes may represent one Arabidopsis source gene.

## FeaturePlots

FeaturePlots for all seven OSA auxin trans-panels were created using
`umap.allgenes`. The colour maximum for each panel was set to its 98th
percentile to reduce the influence of extreme cells.

Outputs:

```text
data/wang_rice_root_confident_Auxin_OSA_FeaturePlots/
```

This folder contains:

- a combined seven-panel PDF and PNG;
- individual PDF and PNG files for each panel.

## Recommended object selection

- Use `shahan_balanced_full_root_full_genes.rds` for balanced Arabidopsis
  full-root analyses requiring the complete gene set.
- Use `wang_rice_root_confident_transferred_all_genes.rds` for rice analyses
  requiring genes beyond the auxin panels.
- Use `wang_rice_root_confident_Auxin_OSA_iSensors.rds` for rapid iSensor
  plotting and modelling.
- Return to the 31,108-cell all-gene transfer object when uncertain cells or
  transfer-confidence sensitivity analyses are required.

## Limitations

The 17,201-cell Wang objects contain only confidently transferred cells.
Analyses based on them therefore describe the confidently annotatable subset,
not the entire Wang root atlas. Vascular cells are underrepresented because
their fine-label transfer acceptance was lower than for most other tissues.

