# Simulation

- 1.1 & 1.2 & 1.3 Generate synthetic data based on three reference datasets,the reference data is in 'reference' folder .
- plot_code Figure 2 synthetic data.
- tm_plot Figure 4 D E
- method 11 methods code
- results_summary Calculate statistics based on the results.

## Dependency
The scripts in this repository are excuted in R (v4.3.1) and use the following R packages:
ggplot2==3.3.4, Seurat==4.4.0, cpm=2.3, PRROC==1.3.1, ggpubr==0.6.0, cowplot==1.1.1, ggforce==0.4.1, SRTsim==0.99.6, reticulate==1.34.0, latex2exp==0.9.6

The scripts in 'method' folder are excuted in R(v4.0.5) and python(v3.7.10),and use the following packages:
### R packages
Loading data recording time memory: Seurat==4.0.1, Giotto==1.0.4, tidyverse==1.3.1, patchwork==1.1.1, dplyr==1.0.7, ggplot2==3.3.5, cowplot==1.1.1, parallel==4.0.5, peakRAM==1.0.2, SeuratData==0.2.1
Run spark/sparkX: SPARK==1.1.1
Run RV: FactoMineR==2.4
Run meringue: MERINGUE==1.0
Run Dcor: Rfast==2.0.6
Run HSIC: dHSIC==2.1

### python packages
Loading data recording time memory:numpy==1.21.5, pandas==1.5.3, memory-profiler==0.61.0
Run sepal: sepal==1.0.0
Run scGCO: scGCO==1.1.0
Run spatialDE:SpatialDE==1.1.3
Run SOMDE: somde==0.1.8