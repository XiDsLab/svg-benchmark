# Real analysis

- adjacent_reproducibility.R Obtain reproducibility of SVG identification algorithms by calculating the Jaccard similarity between the top-K identified SVGs in adjacent slices of spatial transcriptomics.

- generation_of_pseudo_datasets.R To estimate the type I errors of SVG identification algorithms, we generate pseudo datasets (take 10x Visium as an example) by randomly permuting each gene¡¯s expression values across the spatial locations.

- similarity_between_methods.R Obtain between-method-similarity by calculating the Jaccard similarity of the top-K genes identified by different SVG identification algorithms on real datasets.

- preprocess.R Quality control and Normalize process.

## clustering 

There are five clustering algorithms, among which BayesSpace and spaGCN are spatial clustering algorithms, while Lovain, LVM, and SLM are scRNA-seq clustering algorithms. Given the expert annotations, we tune the parameters of these five clustering methods, including the ¡°resolution¡± parameter of the scRNA-seq clustering methods, the ¡°q¡± parameter of BayesSpace, the ¡°target_num¡± parameter of SpacGCN, such that the resulting number of clusters equals to the number of clusters in the expert annotations.

## construction_of_silver_standards
- histological_signals.R Silver standards constructed by correlating with histology images. Taking the results from LABTrans.py as input, we calculate the Pearson's correlation and then obtain silver standard SVGs.

- LABTrans.py Convert the RGB image data to lightness vectors through the LAB transform.

- Moran_index.R Silver standards constructed using spatial aggregation. Dependency: R packages spdep==1.2_4, moranfast==1.0

- NBregression.R Silver standards constructed using the NB regression. Dependency: R packages MASS==7.3_54, stats==4.0.5

- Wilcoxon_test.R Silver standards constructed using the Wilcoxon test. Dependency: R packages stats==4.0.5, ICSKAT==0.2.0

## stability
- generation_of_perturbed_datasets.R To evaluate SVG methods¡¯ robustness against the ¡°spot-swapping¡± perturbation, we generate perturbed datasets (take 10x Visium as an example).