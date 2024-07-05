Source code to reproduce the results described in the article: "Benchmarking algorithms for spatially variable gene identification in spatial transcriptomics".

We include 15 SVG identification methods:

12 SVG identification methods: BinSpect, BOOST-GP, BOOST-MI, GPcounts, MERINGUE, scGCO, sepal, SOMDE, SPARK, SPARK-X, SpatialDE, and trendsceek,

3 additional general-purpose multivariate-correlation methods: RV-coefficient, distance correlation (dCor) coefficient and Hilbert Schmidt Independent Criterion (HSIC).

The codes include spatially variable gene identification algorithms ("./Methods"), synthetic data generation ("./Simulation"), running time and memory usage records ("./Scalability"), and analysis of real data involved ("./Real analysis").

Data description: The reference datasets for synthetic data generation are under "./Simulation/Reference", while the real-world datasets are available at https://doi.org/10.5281/zenodo.7227772 .
