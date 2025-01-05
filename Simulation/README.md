# Simulation

- 1.1 & 1.2 & 1.3_generation_refX.R Generate synthetic data based on three reference datasets.

- reference Store reference real datasets for synthetic data generation.

- results_summary.R Calculate statistics based on the results.

- 30 synthetic datasets is available at https://doi.org/10.5281/zenodo.7227771

## Dependency
The scripts in this repository are excuted in R (v4.3.1) and use the following R packages:
ggplot2==3.3.4, Seurat==4.4.0, cpm=2.3, PRROC==1.3.1, ggpubr==0.6.0, cowplot==1.1.1, ggforce==0.4.1, SRTsim==0.99.6, reticulate==1.34.0, latex2exp==0.9.6

The scripts in 'method' folder are excuted in R(v4.0.5) and python(v3.7.10),and use the following packages:
### R packages
Loading data recording time memory: Seurat==4.0.1, Giotto==1.0.4, tidyverse==1.3.1, patchwork==1.1.1, dplyr==1.0.7, ggplot2==3.3.5, cowplot==1.1.1, parallel==4.0.5, peakRAM==1.0.2, SeuratData==0.2.1

### python packages
Loading data recording time memory:numpy==1.21.5, pandas==1.5.3, memory-profiler==0.61.0
