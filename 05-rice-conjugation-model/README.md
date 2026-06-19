# Rice auxin model with expanded ConjugationDeconjugation panel

This analysis refits the Figure 7 rice linear model after replacing the old
`OSA-aux-trans-IAA` predictor with the local expanded
`OSA-aux-trans-ConjugationDeconjugation` panel.

Run from `iSensors-supplementary/`:

```r
source("05-rice-conjugation-model/rerun_rice_linear_model.R")
```

Inputs are read from `00-iSensors-objects/data/`. Tables are written to `in/`
and plots to `out/`. Existing Figure 7 files are not modified.
