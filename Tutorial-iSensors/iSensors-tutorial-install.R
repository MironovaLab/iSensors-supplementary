# Core packages
install.packages(c("Seurat", "Matrix", "tidyverse", "magrittr", "purrr", "stringr", devtools))

# Bioconductor dependency used by iSensors for generating cis-sensors
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")

devtools::install_github("MironovaLab/iSensors")