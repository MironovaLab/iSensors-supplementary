# Session report: rice auxin modelling and Conjugation panel correction

## Objective

This session investigated why the original rice model gave misleading auxin
predictions in root tissues, particularly cortex. The analysis focused on
whether the OSA PAT panel was incomplete and whether the model had used an
appropriate ConjugationDeconjugation panel.

The model considered was:

```text
ARF = beta1 * Synthesis + beta2 * Conjugation + beta3 * PAT
```

where ARF iSensor expression was treated as the response and the three
auxin-related trans-panels as predictors.

## Data used

Two related Wang rice datasets were analysed:

1. The complete rice atlas represented by
   `00-iSensors-objects/data/Wang_rice_iSensors.rds`.
   It contains 115,395 cells and 56 organ-cell-type groups.
2. The confidently reannotated root dataset represented by
   `00-iSensors-objects/data/wang_rice_root_confident_Auxin_OSA_iSensors.rds`.
   It contains 17,201 cells with Arabidopsis `Atlas_new` labels transferred to
   rice.

The existing normalized RNA data and previously calculated ARF, Synthesis and
PAT scores were reused. Only the corrected ConjugationDeconjugation score had
to be introduced.

## Initial diagnostic

The original Figure 7 rice model used:

```text
OSA-aux-trans-IAA
```

as the `Conjugation` predictor. This is not equivalent to the expanded
ConjugationDeconjugation panel developed during the cross-species analysis.

In the original model, root cortex had:

| Quantity | Value |
|---|---:|
| Observed ARF | 0.154 |
| Predicted ARF | 0.227 |
| Residual, observed - predicted | -0.073 |

Thus, the original model **overpredicted**, rather than underpredicted, root
cortex ARF. The apparent direction depended on correctly interpreting the
residual definition.

## PAT panel audit

The installed OSA PAT panel contains 18 rice genes. All 18 were present in the
Wang RNA assay, so the panel was technically complete relative to the iSensors
registry.

The panel predominantly represents PIN efflux and AUX/LAX influx transporters.
However, its effective single-cell coverage was limited:

- Four genes were not detected in any confidently reannotated root cell.
- Most panel genes were detected in fewer than 1% of cells.
- Approximately 65% of the cortex PAT expression signal came from OsPIN2 and
  PIN1A.
- OsLAX1 contributed approximately another 19%.

Consequently, the nominal 18-gene PAT panel behaved largely as a three-gene
signal in cortex. No smaller PIN/AUX-LAX panels were created or retained.

Audit tables:

- `00-iSensors-objects/data/OSA_PAT_panel_annotation.tsv`
- `00-iSensors-objects/data/OSA_PAT_Wang_expression_coverage.tsv`
- `00-iSensors-objects/data/OSA_PAT_Wang_cortex_gene_contributions.tsv`

## Corrected ConjugationDeconjugation panel

The corrected local panel contains 20 rice genes representing 14 Arabidopsis
source genes. It includes high-confidence Ensembl Plants one-to-one,
one-to-many and many-to-many orthologues of the Arabidopsis
ConjugationDeconjugation panel.

Mapping source:

```text
00-iSensors-objects/data/
OSA-aux-trans-ConjugationDeconjugation_all_high_confidence.tsv
```

## Complete rice-atlas model

The model was rerun on the same 56 organ-cell-type means used for Figure 7.
The averaging procedure was matched exactly to Seurat's
`AverageExpression(layer = "data")` convention by exponentiating normalized
values before averaging.

### Model comparison

| Model | Adjusted R2 | RMSE | Residual SE | AIC |
|---|---:|---:|---:|---:|
| Old OSA IAA panel | 0.843 | 0.0514 | 0.0528 | -165.6 |
| Expanded ConjugationDeconjugation | 0.883 | 0.0444 | 0.0457 | -181.9 |

The expanded panel improved model fit across all reported measures.

### Corrected model coefficients

| Predictor | Estimate | Significance |
|---|---:|---:|
| Synthesis | 2.198 | p = 3.63e-05 |
| Conjugation | 0.781 | p = 1.14e-07 |
| PAT | 0.554 | p = 0.000925 |

Unlike the old model, PAT became statistically significant after replacing
the IAA panel with the expanded ConjugationDeconjugation panel.

At the organ level, the mean root residual was approximately zero:

```text
mean root residual = +0.0032
```

For the root-cortex atlas group:

| Quantity | Old model | Corrected model |
|---|---:|---:|
| Observed ARF | 0.154 | 0.154 |
| Predicted ARF | 0.227 | 0.131 |
| Residual | -0.073 | +0.023 |

This demonstrates that the old IAA panel was responsible for much of the
misleading root-cortex overprediction.

## Reannotated root model

The corrected model was then fitted to mean iSensor expression for transferred
`Atlas_new` root subtypes.

Subtypes represented by at least 10 cells were retained:

- 23 robust subtypes
- 17,180 of 17,201 cells

Six very small labels were excluded from the primary model:

- PPP_t
- Protoxylem_e2
- LRCEI
- PCC_e1
- PSE_t
- Trichoblast_m2

### No-intercept root model

The fitted coefficients were:

| Predictor | Estimate | p-value |
|---|---:|---:|
| Synthesis | -0.224 | 0.774 |
| Conjugation | 3.035 | 0.015 |
| PAT | 0.050 | 0.861 |

The model reported:

```text
adjusted R2 = 0.884
RMSE = 0.0624
```

All three cortex zones had modest positive residuals:

| Subtype | Residual |
|---|---:|
| Cortex_m1 | +0.0179 |
| Cortex_m2 | +0.0078 |
| Cortex_t | +0.0339 |

The mean cortex residual was +0.0198. The largest model discrepancies occurred
in xylem, especially Protoxylem_t and Metaxylem_d.

## Why the apparent R2 was misleading

The root model was deliberately fitted without an intercept, following the
original Figure 6/7 formulation. In a no-intercept model, R2 is calculated
relative to zero rather than relative to the mean response.

All iSensor scores are positive and share a substantial positive baseline.
This inflates the uncentred R2.

For the reannotated root data:

| Statistic | Value |
|---|---:|
| No-intercept adjusted R2 | 0.884 |
| Observed-predicted correlation | 0.611 |
| Squared observed-predicted correlation | 0.374 |
| Adjusted R2 with an intercept | 0.344 |

The conventional predictive interpretation is therefore approximately
34-37% explained variation, not 88%.

### Predictor collinearity

The predictors were strongly correlated:

| Predictor pair | Correlation |
|---|---:|
| Synthesis-Conjugation | 0.850 |
| Conjugation-PAT | 0.717 |
| Synthesis-PAT | 0.560 |

This multicollinearity makes the individual coefficients unstable because the
three panels carry overlapping information.

When an intercept was included, none of the individual predictors was
significant:

| Predictor | Estimate | p-value |
|---|---:|---:|
| Synthesis | 0.476 | 0.474 |
| Conjugation | 1.486 | 0.164 |
| PAT | -0.198 | 0.418 |

The intercept itself was significant and estimated at 0.109.

## Interpretation

The expanded ConjugationDeconjugation panel is clearly preferable to using the
OSA IAA panel as a surrogate for conjugation. It improves the complete-atlas
model and removes the strong root-cortex overprediction.

However, the root-only model does not support the conclusion that
Conjugation alone mechanistically determines ARF activity. Its apparent
dominance in the no-intercept model is affected by:

- a shared positive baseline;
- strong correlations among the predictors;
- a relatively small number of aggregated subtype observations;
- the use of transcript abundance as a proxy for transport and metabolism;
- lack of information about PIN protein localization and polarity.

The defensible conclusion is:

> The corrected auxin panels collectively reproduce part of the variation in
> ARF iSensor expression across rice root subtypes. The expanded
> ConjugationDeconjugation panel improves prediction substantially, but the
> independent contributions of Synthesis, Conjugation and PAT cannot be
> resolved reliably because their expression scores are strongly correlated.

## Generated analysis files

Main scripts:

- `rerun_rice_linear_model.R`
- `rerun_reannotated_root_model.R`

Complete-atlas plots:

- `out/Rice_ARF_model_expanded_conjugation_combined.pdf`
- `out/Rice_ARF_scatter_expanded_conjugation.pdf`
- `out/Rice_ARF_coefficients_expanded_conjugation.pdf`
- `out/Rice_ARF_residual_forests_expanded_conjugation.pdf`

Reannotated-root plots:

- `out/Reannotated_root_ARF_model_combined.pdf`
- `out/Reannotated_root_ARF_scatter.pdf`
- `out/Reannotated_root_ARF_coefficients.pdf`
- `out/Reannotated_root_ARF_residual_forests.pdf`

Detailed model tables and summaries are available in the `in/` directory.

## Recommended next steps

For manuscript-level reporting:

1. Retain the expanded ConjugationDeconjugation panel.
2. Report ordinary observed-predicted correlation or centred R2 rather than
   the uncentred no-intercept R2.
3. Treat individual coefficients as descriptive because of collinearity.
4. Present residuals by transferred subtype to show where the model succeeds
   and fails.
5. Avoid interpreting PAT transcript expression as directional auxin flux.

