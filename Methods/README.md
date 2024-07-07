# Methods
We include 15 SVG identification methods,

12 SVG identification methods: BinSpect, BOOST-GP, BOOST-MI, GPcounts, MERINGUE, scGCO, sepal, SOMDE, SPARK, SPARK-X, SpatialDE, and trendsceek,

3 additional general-purpose multivariate-correlation methods: RV-coefficient, distance correlation (dCor) coefficient and Hilbert Schmidt Independent Criterion (HSIC),

More details see section Methods in our article. Below are the code sources and references for these SVG methods.

## BinSpect
Source code: https://rubd.github.io/Giotto_site/

Reference to: Dries, R., Zhu, Q., Dong, R., Eng, C. H. L., Li, H., Liu, K., ... & Yuan, G. C. (2021). Giotto: a toolbox for integrative analysis and visualization of spatial expression data. *Genome biology*, 22, 1-31.

## BOOST-GP
Source code: https://github.com/Minzhe/BOOST-GP

Reference to: Li, Q., Zhang, M., Xie, Y., & Xiao, G. (2021). Bayesian modeling of spatial molecular profiling data via Gaussian process. *Bioinformatics, 37(22)*, 4129-4136.

## BOOST-MI
Source code: https://github.com/Xijiang1997/BOOST-MI

Reference to: Jiang, X., Xiao, G., & Li, Q. (2022). A Bayesian modified Ising model for identifying spatially variable genes from spatial transcriptomics data. *Statistics in medicine, 41(23)*, 4647-4665.

## GPcounts
Source code: https://github.com/ManchesterBioinference/GPcounts

Reference to: BinTayyash, N., Georgaka, S., John, S. T., Ahmed, S., Boukouvalas, A., Hensman, J., & Rattray, M. (2021). Non-parametric modelling of temporal and spatial counts data from RNA-seq experiments. *Bioinformatics, 37(21)*, 3788-3795.

## MERINGUE
Source code: https://github.com/JEFworks-Lab/MERINGUE

Reference to: Miller, B. F., Bambah-Mukku, D., Dulac, C., Zhuang, X., & Fan, J. (2021). Characterizing spatial gene expression heterogeneity in spatially resolved single-cell transcriptomic data with nonuniform cellular densities. *Genome research, 31(10)*, 1843-1855.

## scGCO
Source code: https://github.com/WangPeng-Lab/scGCO

Reference to: Zhang, K., Feng, W., & Wang, P. (2022). Identification of spatially variable genes with graph cuts. *Nature Communications, 13(1)*, 5488.

## sepal
Source code: https://github.com/almaan/sepal

Reference to: Andersson, A., & Lundeberg, J. (2021). sepal: Identifying transcript profiles with spatial patterns by diffusion-based modeling. *Bioinformatics, 37(17)*, 2644-2650.

## SOMDE
Source code: https://github.com/WhirlFirst/somde

Reference to: Hao, M., Hua, K., & Zhang, X. (2021). SOMDE: a scalable method for identifying spatially variable genes with self-organizing map. *Bioinformatics, 37(23)*, 4392-4398.

## SPARK
Source code: https://github.com/xzhoulab/SPARK

Reference to: Sun, S., Zhu, J., & Zhou, X. (2020). Statistical analysis of spatial expression patterns for spatially resolved transcriptomic studies. *Nature methods, 17(2)*, 193-200.

## SPARK-X
Source code: https://github.com/xzhoulab/SPARK

Reference to: Zhu, J., Sun, S., & Zhou, X. (2021). SPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies. *Genome biology, 22(1)*, 184.

## SpatialDE
Source code: https://github.com/Teichlab/SpatialDE

Reference to: Svensson, V., Teichmann, S. A., & Stegle, O. (2018). SpatialDE: identification of spatially variable genes. *Nature methods, 15(5)*, 343-346.

## trendsceek
Source code: https://github.com/edsgard/trendsceek

Reference to: Edsgard, D., Johnsson, P., & Sandberg, R. (2018). Identification of spatial expression trends in single-cell gene expression data. *Nature methods, 15(5)*, 339-342.

## RV
Source code: https://github.com/cran/FactoMineR

Reference to: Escoufier, Y. (1973). Le traitement des variables vectorielles. *Biometrics*, 751-760.

## dCor
Source code: https://github.com/cran/Rfast

Reference to: Sz¨¦kely, G. J., Rizzo, M. L., & Bakirov, N. K. (2007). Measuring and testing dependence by correlation of distances.

## HSIC
Source code: https://github.com/cran/dHSIC

Reference to: Gretton, A., Herbrich, R., Smola, A., Bousquet, O., Scholkopf, B., & Hyvarinen, A. (2005). Kernel methods for measuring independence. *Journal of Machine Learning Research, 6(12)*.

## Dependency
For the algorithms that specifically developed for SVG detection, the softwares we used are 

BinSpect (Giotto v1.0.4), SPARK (SPARK v1.1), MERINGUE (MERINGUE v1.0), SpatialDE (SpatialDE v1.1.3), SOMDE (somde v0.1.8), SPARK-X (SPARK v1.1.1), sepal (sepal v1.0.0), scGCO (scGCO v1.1.0), trendsceek (trendsceek v1.0.0), BOOST-GP (Github documentation download required), BOOST-MI (Github documentation download required), GPcounts (GPcounts v0.1). For the three general-purpose correlation methods, we use R packages FactoMineR (v2.4), Rfast (v2.0.6) and dHSIC (v2.1) to calculate the correlations RV, dCor and HSIC between spatial locations¡¯ coordinates and gene expression.

For time and memory recording, R packages parallel==4.0.5, peakRAM==1.0.2; python package memory-profiler==0.61.0
