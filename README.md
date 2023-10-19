# Metabolic subtyping of patient tumors from gene expression data

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[LinkedOmics](https://www.linkedomics.org/login.php) is a great resource of several cancer types omics data. Data sets in the `data/` folder were downloaded from there.

## Features

1. Run GSVA data on cancer patient gene expression data using KEGG metabolic pathways.
2. Perform k-means clustering on the GSVA matrix, identify metabolic subtypes. Optimal number of k is defined based on data.
3. Planned: Summarize pathway activity per cluster: e.g. mean pathway activity per cluster, or do differential testing.
4. Perform KM analysis with clusters if survival data is present.

## Test data sets

- [CPTAC-HNSCC](https://www.linkedomics.org/data_download/CPTAC-HNSCC/): `data/CPTAC_HNSCC/`
- [CPTAC-COAD](https://www.linkedomics.org/data_download/CPTAC-COAD/): `data/CPTAC_COAD/`
- [TCGA-COAD](https://www.linkedomics.org/data_download/TCGA-COADREAD/): `data/TCGA_COAD/`
- [TCGA-LAML](https://www.linkedomics.org/data_download/TCGA-LAML/): `data/TCGA_LAML/`
- [TCGA-UVM](https://www.linkedomics.org/data_download/TCGA-UVM/): `data/TCGA_UVM/`
- [TCGA-ACC](https://www.linkedomics.org/data_download/TCGA-ACC/): `data/TCGA_ACC/`
