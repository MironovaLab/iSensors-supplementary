# Arabidopsis-to-rice root reannotation report

## Objective

The purpose of this analysis was to obtain comparable root cell subtypes in
Arabidopsis and rice. The Wang rice atlas had broad root annotations, whereas
the Shahan Arabidopsis dataset contained the more detailed `Atlas_new`
annotation. The Shahan labels were therefore transferred to Wang rice root
cells using one-to-one orthologous genes.

The workflow was designed to avoid loading and integrating both complete,
multi-gigabyte Seurat objects simultaneously.

## Source datasets

### Arabidopsis reference

```text
00-iSensors-objects/data/02_ra_integrated_Shahan_20250114.rds
```

The complete Shahan root dataset was used, not only the root tip. The
`Atlas_new` metadata field supplied the reference labels.

### Rice query

```text
00-iSensors-objects/data/Wang_rice.Rds
```

Only cells with `tissue == "Root"` were selected. This produced 31,108 Wang
root cells.

## Memory-conscious strategy

The workflow was divided into restartable phases:

1. Retrieve and cache orthologues.
2. Build a balanced Arabidopsis reference.
3. Build an orthologue-only rice root query.
4. Run a stratified pilot transfer.
5. Run the complete transfer in chunks.
6. Reconstruct an all-gene rice object and calculate rice-only PCA/UMAP.

This prevented the two original datasets from coexisting in memory.

## Arabidopsis reference balancing

Cells were sampled separately for every non-empty `Atlas_new` label:

- Maximum cells per label: 500
- Random seed: 100
- Available annotated cells: 53,767
- Retained cells: 23,494
- Retained `Atlas_new` labels: 84

Labels represented by fewer than 500 cells were retained completely. This
preserved the complete set of developmental zones and rare labels while
preventing abundant labels from dominating anchor identification.

The label-level sampling table is:

```text
00-iSensors-objects/data/arabidopsis_rice_transfer/reference_sampling.tsv
```

## Orthologue mapping

High-confidence one-to-one Arabidopsis-rice orthologues were downloaded from
Ensembl Plants through `biomaRt`.

Selection criteria:

- Orthology type: `ortholog_one2one`
- Orthology confidence: 1
- Arabidopsis and rice identifiers unique after removing version suffixes

The resulting table contained 3,283 one-to-one orthologue pairs:

```text
00-iSensors-objects/data/arabidopsis_rice_transfer/
arabidopsis_rice_one2one.tsv
```

The balanced reference contained 23,494 cells and 3,283 genes. The rice query
contained 31,108 cells and 3,282 matched genes.

Rice orthologues were renamed with their corresponding Arabidopsis identifiers
inside the query object, creating a shared feature namespace for transfer.

## Label-transfer procedure

Transfer was performed independently within broad Wang-derived root lineages:

- Epidermis
- Cortex
- Endodermis
- Vascular, with Arabidopsis pericycle labels also permitted
- Root cap
- Meristem/initials

Broad lineage mappings were constructed from text patterns in the existing
Wang annotation and the Shahan `Atlas_new` labels.

For each lineage:

1. Relevant reference cells were selected.
2. State-like labels such as cell-cycle states, stressed states and unknown
   states were excluded from the reference.
3. Query cells were normalized.
4. Transfer anchors were calculated using PCA projection.
5. `TransferData` generated label scores.
6. The highest and second-highest scores were retained.

Main parameters:

| Parameter | Value |
|---|---:|
| Variable features | up to 3,000 |
| PCA dimensions | 30 |
| `k.anchor` | 5 |
| Maximum `k.weight` | 25 |
| Full-transfer chunk size | 25,000 cells |
| Pilot size | 20,000 cells |
| Random seed | 100 |

## Confidence filtering

Four principal fields were stored:

- `predicted.Atlas_new_raw`: highest-scoring label
- `prediction.score.max`: highest label score
- `prediction.score.margin`: difference between the first and second scores
- `predicted.Atlas_new`: accepted label or `Uncertain`

A fine label was accepted when:

```text
maximum score >= 0.50
margin >= 0.10
broad lineage compatible
```

Because transfer was run separately inside the Wang-derived broad lineages,
all 31,108 cells were broad-compatible by construction. Broad compatibility
was therefore a safety check, not an independent validation. Acceptance was
primarily determined by the score and margin thresholds.

## Transfer results

| Result | Value |
|---|---:|
| Wang root cells | 31,108 |
| Confident fine-label assignments | 17,201 |
| Uncertain assignments | 13,907 |
| Accepted fraction | 55.3% |
| Raw transferred labels | 52 |
| Accepted transferred labels | 29 |

Acceptance varied strongly among the original Wang annotations:

| Original Wang type | Cells | Accepted | Accepted percent |
|---|---:|---:|---:|
| Exodermis | 2,902 | 2,521 | 86.9% |
| Root cap | 903 | 765 | 84.7% |
| Meristem | 3,770 | 2,777 | 73.7% |
| Endodermis | 1,378 | 902 | 65.5% |
| Epidermis near root hair | 5,008 | 3,251 | 64.9% |
| Sclerenchyma | 2,348 | 1,517 | 64.6% |
| Cortex | 2,733 | 1,380 | 50.5% |
| Epidermis | 8,489 | 3,347 | 39.4% |
| Vascular cylinder | 3,577 | 741 | 20.7% |

The lower vascular acceptance indicates that vascular fine-subtype transfer
was considerably less certain than transfer in root cap, cortex-associated or
meristematic lineages.

The complete prediction table is:

```text
00-iSensors-objects/data/arabidopsis_rice_transfer/
wang_rice_root_predictions.tsv
```

## Rice all-gene PCA, UMAP and clustering

The orthologue-only UMAP was useful for checking transfer mechanics, but it did
not represent the complete rice transcriptome. The final visualization was
therefore rebuilt from all 37,960 rice RNA genes.

Processing:

- Log normalization
- 3,000 variable genes
- 50 principal components calculated
- PCA dimensions 1-30 used downstream
- UMAP with 30 neighbours, minimum distance 0.15 and cosine metric
- Graph clustering at resolution 0.6
- Random seed 100

This produced 17 de novo rice RNA clusters.

The resulting object is:

```text
00-iSensors-objects/data/arabidopsis_rice_transfer/
wang_rice_root_label_transferred_all_genes_umap.rds
```

It contains:

- 31,108 cells
- 37,960 RNA genes
- RNA counts and normalized data
- PCA
- `umap.allgenes`
- transferred labels and confidence metadata
- 17 de novo clusters

## Annotation comparison statistics

| Comparison | Cells | Adjusted Rand index | Normalized mutual information |
|---|---:|---:|---:|
| All-gene clusters vs original Wang annotation | 31,108 | 0.531 | 0.655 |
| All-gene clusters vs accepted transferred labels | 17,201 | 0.301 | 0.530 |
| Original Wang annotation vs accepted transferred labels | 17,201 | 0.556 | 0.728 |

These measures describe association between label systems. They are not an
independent test of biological truth because:

- transferred labels were constrained by Wang-derived broad lineages;
- only confident cells were used in some comparisons;
- transferred labels originated from a different species.

## Main visual outputs

All-gene plots and statistics are stored in:

```text
00-iSensors-objects/data/arabidopsis_rice_transfer/all_gene_umap/
```

Important plots:

- `rice_all_genes_clusters_vs_original.pdf`
- `rice_all_genes_broad_vs_fine_transfer.pdf`
- `rice_all_genes_transferred_fine_faceted.pdf`
- `rice_all_genes_de_novo_clusters.pdf`
- `rice_all_genes_original_annotation.pdf`
- `rice_all_genes_transferred_broad.pdf`

The orthologue-only diagnostic UMAPs are stored in:

```text
00-iSensors-objects/data/arabidopsis_rice_transfer/dimplots/
```

## Interpretation and limitations

The procedure supplied a useful Arabidopsis-compatible fine annotation for
17,201 Wang rice root cells. It allowed direct comparison of corresponding
epidermal, cortex, endodermal, pericycle, xylem, root-cap, lateral-root and
stem-cell-niche subtypes.

The following limitations should remain explicit:

1. Fine labels are inferred by cross-species transcriptomic similarity rather
   than experimentally validated lineage tracing.
2. Broad lineage agreement is partly imposed by the transfer design.
3. Approximately 44.7% of rice root cells did not satisfy the confidence
   thresholds.
4. Vascular fine labels had particularly low coverage.
5. One-to-one orthologues improve interpretability but omit lineage-specific
   paralogues and species-specific genes.
6. The all-gene UMAP is a visualization of rice transcriptional structure; it
   was not used to generate the transferred labels.

## Reproducibility

Primary scripts:

- `arabidopsis_rice_label_transfer.R`
- `plot_transferred_rice_dimplot.R`
- `plot_transferred_rice_all_genes.R`

The scripts use fixed seed 100 and preserve intermediate objects so that
orthologue retrieval, reference preparation, query preparation and transfer
can be rerun independently.

