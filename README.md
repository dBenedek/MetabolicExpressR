# Metabolic subtyping of cancer patient data from gene expression data

## Plans

1. Run GSVA data on cancer patient gene expression data using KEGG metabolic pathways.
2. Perform k-means clustering on the GSVA matrix, identify metabolic subtypes. Optimal number of k has to be defined based on data.
3. Summarize pathway activity per cluster: e.g. mean pathway activity per cluster, or do differential testing.
4. Perform KM analysis with clusters if survival data is present.
