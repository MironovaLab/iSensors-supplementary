# Arabidopsis-to-rice root label transfer

This pipeline transfers the full-root Shahan `Atlas_new` annotation to root
cells from the Wang rice atlas without integrating the two complete Seurat
objects in memory.

The Arabidopsis reference is balanced before normalization and PCA. By default,
at most 500 cells are retained per `Atlas_new` label; labels with fewer cells
are retained completely. All developmental zones present in the source object
remain eligible.

## Inputs

- `../00-iSensors-objects/data/02_ra_integrated_Shahan_20250114.rds`
- `../00-iSensors-objects/data/Wang_rice.Rds`
- A two-column, one-to-one Arabidopsis/rice ortholog table. The script can
  create this with Ensembl Plants and `biomaRt`, or use a cached TSV.

## Run

From `iSensors-supplementary/`:

```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" `
  "04-cross-species-integration/arabidopsis_rice_label_transfer.R" `
  --phase orthologs

& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" `
  "04-cross-species-integration/arabidopsis_rice_label_transfer.R" `
  --phase reference

& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" `
  "04-cross-species-integration/arabidopsis_rice_label_transfer.R" `
  --phase query

& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" `
  "04-cross-species-integration/arabidopsis_rice_label_transfer.R" `
  --phase pilot
```

Inspect the pilot outputs in
`00-iSensors-objects/data/arabidopsis_rice_transfer/`. If they are sensible,
run:

```powershell
& "C:\Program Files\R\R-4.5.1\bin\Rscript.exe" `
  "04-cross-species-integration/arabidopsis_rice_label_transfer.R" `
  --phase full
```

The `orthologs` phase requires `biomaRt`. Install it once if necessary:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("biomaRt")
```

Alternatively, place a tab-separated file at
`00-iSensors-objects/data/arabidopsis_rice_transfer/arabidopsis_rice_one2one.tsv`
with columns `arabidopsis_gene` and `rice_gene`. Version suffixes such as
`.1` are removed automatically.

## Important options

Options can follow the phase argument:

```text
--cells-per-label 500
--pilot-cells 20000
--chunk-size 25000
--nfeatures 3000
--dims 30
--min-score 0.50
--min-margin 0.10
--seed 100
```

Lower `--chunk-size` if transfer still approaches the machine's memory limit.
Lower `--cells-per-label` to 300 for a smaller reference.

## Outputs

- Cached one-to-one ortholog table
- Balanced, ortholog-only Arabidopsis reference
- Ortholog-only rice root query
- Pilot predictions and validation tables
- Full rice-root prediction metadata
- Final slim rice-root Seurat object with transferred labels

`predicted.Atlas_new_raw` is the best label regardless of confidence.
`predicted.Atlas_new` is set to `Uncertain` unless its score and score margin
pass the configured thresholds and its broad lineage is compatible with the
existing rice annotation.

