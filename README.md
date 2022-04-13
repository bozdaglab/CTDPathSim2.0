# CTDPathSim2.0
R package CTDPathSim2.0 computes the similarity scores between patient tumor samples and cancer cell lines at genetic, genomic, and epigenetic levels integrating multi-omics datasets.
             It has five computational steps:
             Step 1. Computing sample-specific deconvoluted DNA methylation profile.
             Step 2. Computing sample-specific deconvoluted expression profile.
             Step 3: Computing sample and cell line-specific differentially expressed (DE) genes and biological pathways.
             Step 4: Computing sample and cell line-specific differentially methylated (DM) and differentially aberrated (DA) genes.
             Step 5: Computing sample-cell line pathway activity-based similarity score.
             
# Installing the package
install.packages("remotes")

remotes::install_github("boseb/CTDPathSim2.0")

# Check package vignette to run the tool
library(CTDPathSim2.0)

browseVignettes("CTDPathSim2.0")
