# Metabolic subtyping of patient tumors from gene expression data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A detailed tutorial can be found [here](https://github.com/dBenedek/genexp_metabolic_subtyping/blob/main/vignettes/tutorial.Rmd).

## Features

1. Perform [GSVA](https://doi.org/10.1186/1471-2105-14-7) on cancer patient gene expression data using [KEGG metabolic pathways](https://www.genome.jp/kegg/pathway.html#metabolism).
2. Perform k-means clustering on the GSVA matrix, identify metabolic subtypes. Optimal number of k is defined based on data or user specified k is used.
3. Summarize pathway activity per cluster: i.e. mean pathway activity per cluster.
4. Perform Kaplan-Meier (KM) analysis with clusters if survival data is present.

## Test data sets

[LinkedOmics](https://www.linkedomics.org/login.php) is a great resource of several cancer types omics data. Data sets in the `data/` folder were downloaded from there.

- [CPTAC-HNSCC](https://www.linkedomics.org/data_download/CPTAC-HNSCC/): `data/CPTAC_HNSCC/`
- [CPTAC-COAD](https://www.linkedomics.org/data_download/CPTAC-COAD/): `data/CPTAC_COAD/`
- [TCGA-COAD](https://www.linkedomics.org/data_download/TCGA-COADREAD/): `data/TCGA_COAD/`
- [TCGA-LAML](https://www.linkedomics.org/data_download/TCGA-LAML/): `data/TCGA_LAML/`
- [TCGA-UVM](https://www.linkedomics.org/data_download/TCGA-UVM/): `data/TCGA_UVM/`
- [TCGA-ACC](https://www.linkedomics.org/data_download/TCGA-ACC/): `data/TCGA_ACC/`
