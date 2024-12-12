# MetabolicExpressR: Metabolic subtyping of patient tumors from gene expression data  <img src="badge.png" align="right" height="139" />

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A detailed tutorial can be found [here](https://github.com/dBenedek/genexp_metabolic_subtyping/blob/main/vignettes/tutorial.Rmd).

## Installation

```r
devtools::install_github("dBenedek/MetabolicExpressR")
```

## Features

1. Perform Gene Set Variation Analysis ([GSVA](https://doi.org/10.1186/1471-2105-14-7)) on cancer patient gene expression data using [KEGG metabolic pathways](https://www.genome.jp/kegg/pathway.html#metabolism). Potentially non-cancer data can be used as well.
2. Perform k-means clustering on the GSVA matrix, identify metabolic subtypes. Optimal number of k is defined based on data or user-specified k is used.
3. Summarize pathway "expression" per cluster: i.e. mean pathway activity per cluster.
4. Perform GSVA on cancer patient gene expression data using MSigDB Hallmarks and compare between the identified clusters.    
5. Run [PROGENy](https://saezlab.github.io/progeny/) on the gene expression data of tumors and compare signature activity scores between the identified clusters. 
6. Perform Kaplan-Meier (KM) analysis with clusters if survival data is present.

## Usage 

### Gene Set Variation Analysis 

Gene Set Variation Analysis (GSVA) allows us to measure the enrichment of gene sets on the single sample level. It is recommended to have at least around 10 samples in the data set. In all the examples we use the KEGG metabolic pathways collection, although any other curated gene sets might be used.    
The input of GSVA is a normalized gene expression matrix, where row names are gene names (with the same gene nomenclature as in the gene set used), and column names correspond to sample IDs. The parameter `kcdf` defines the kernel to use during the non-parametric estimation of the cumulative distribution function (see GSVA documentation for details).    

```r
gsva_results <- run_gsva_metabolic(gene_exp_data = gene_exp_cptac_hnscc,
                                   kcdf = "Gaussian",
                                   kegg_gs = kegg_gs,
                                   n_cores = 1L)
dim(gsva_results)
```

The output is an N x M matrix with row names as pathways and column names as sample IDs (same as in the gene expression matrix). The data values are the enrichment scores of the pathways per sample.

### K-means clustering on the GSVA data

The next step is the k-means clustering on the GSVA matrix derived from the previous step using the function `kmeans_gsva_metabolic()`.    
This requires the GSVA matrix from the previous step, the gene set collection (`kegg_gs`, also used in the previous step). The user can either use the data set-specific optimal k for the clustering (default) when `user_def_k = FALSE`, or when `user_def_k = TRUE` then the user has to specify the parameter `k` used for k-means clustering. This way the resulting number of clusters is going to be already pre-defined.   

```r
kmeans_results <- kmeans_gsva_metabolic(gsva_data = gsva_results,
                                        kegg_gs = kegg_gs,
                                        user_def_k = FALSE,
                                        k = NULL)

head(kmeans_results$kmean_res)
kmeans_results$heatmap
``` 

The function returns a list with 1) a data frame with cluster membership per sample (`kmean_res`), and 2) a heatmap (`heatmap`) with GSVA enrichment scores visualized as a heatmap for the clusters. 

### Comparison of pathway "expression" between the identified clusters

Here we compare the pathway enrichment scores (pathway "expression") between the identified clusters. For this, we use the function `plot_mean_pathway_activity()` which compares and visualizes the mean pathway enrichment scores between the clusters. This gives us an overview on the pathway-level differences of the sample clusters.

```r
plot_mean_pathway_activity(gsva_data = gsva_results,
                           kegg_gs = kegg_gs,
                           kmeans_res = kmeans_results$kmean_res)
```

### MSigDB Hallmarks analysis and comparison between the clusters

```r
gsva_hallmarks <- run_gsva_hallmarks(gene_exp_cptac_hnscc,
                                     kcdf = "Gaussian",
                                     n_cores = 1L)
                                     
plot_mean_hallmark_activity(gsva_data = gsva_hallmarks,
                            kmeans_res = kmeans_results$kmean_res,
                            hallmarks_data = hallmarks_data)                                     
```

### PROGENy analysis and signature activity comparison between the identified clusters

From the website of PROGENy: "PROGENy is resource that leverages a large compendium of publicly available signaling perturbation experiments to yield a common core of pathway responsive genes for human and mouse. These, coupled with any statistical method, can be used to infer pathway activities from bulk or single-cell transcriptomics." Thus, it aims to quantify the activity of 14 signatures which are very often disturbed in several cancer entities. These signatures are androgen, EGFR, estrogen, hypoxia, JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb, TNFa, Trail, VEGF, and WNT signaling.

```r
progeny_results <- run_progeny(gene_exp_data = gene_exp, 
				kmeans_res = kmeans_clusters$kmean_res)
progeny_results$progeny_boxplot
```

The output is a boxplot (`progeny_boxplot`) with signaling activity scores per pathway compared between the identified clusters. In addition, differential testing P-value (Wilcoxon or Kruskal-Wallis test) is indicated in the facet header of each signature. The function outputs a data frame (`progeny_res`) as well, with PROGENy pathway activity scores per sample and signature. 

### Kaplan-Meier analysis with the identified clusters

Finally, we perform Kaplan-Meier (KM) analysis between the identified clusters. The function `kmeans_metab_clust_surv()` does exactly this, for which we need survival/clinical endpoint data with e.g. overall survival time and overall survival status of the samples (patients).    
We have to specify the k-means clustering results with the `kmeans_res` parameter (just as in the previous step), we have to provide the clinical data table (`clinical_data` parameter). Moreover, we specify the column name of the sample IDs in the clinical data table (parameter `sample_id_col`), the endpoint status column name (`surv_status_col` parameter), and the endpoint time column name (`surv_time_col` parameter).

```r
kmeans_metab_clust_surv(kmeans_res = kmeans_results$kmean_res, 
                        clinical_data = clinical_data_table,
                        sample_id_col = "case_id", 
                        surv_status_col = "overall_free_status", 
                        surv_time_col = "overall_survival")
```

The output is a KM survival curve of the sample clusters for the selected endpoint.

## Test data sets

[LinkedOmics](https://www.linkedomics.org/login.php) is a great resource of several cancer types omics data. Data sets in the `data/` folder were downloaded from there.

- [CPTAC-HNSCC](https://www.linkedomics.org/data_download/CPTAC-HNSCC/): `data/CPTAC_HNSCC/`
- [CPTAC-COAD](https://www.linkedomics.org/data_download/CPTAC-COAD/): `data/CPTAC_COAD/`
- [TCGA-COAD](https://www.linkedomics.org/data_download/TCGA-COADREAD/): `data/TCGA_COAD/`
- [TCGA-LAML](https://www.linkedomics.org/data_download/TCGA-LAML/): `data/TCGA_LAML/`
- [TCGA-UVM](https://www.linkedomics.org/data_download/TCGA-UVM/): `data/TCGA_UVM/`
- [TCGA-ACC](https://www.linkedomics.org/data_download/TCGA-ACC/): `data/TCGA_ACC/`
