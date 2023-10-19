###############################################################################
# Run metabolic subtype identification from transcriptomic data               #
###############################################################################

# Benedek Danko
# 2023-06-09

setwd("/home/benedek_danko/R/genexp_metabolic_subtyping")

###############################################################################
# Load functions and libraries
###############################################################################

library("tidyverse")

source("R/run_gsva_metabolic.R")
source("R/kmeans_gsva_metabolic.R")
source("R/kmeans_metab_clust_survival.R")

###############################################################################
# Load data
###############################################################################

# CPTAC-HNSCC
gene_exp_cptac_hnscc <- read.table("data/CPTAC_HNSCC/HS_CPTAC_HNSCC_RNAseq_RSEM_UQ_log2_Normal.cct",
                             sep="\t", stringsAsFactors = F, header = T) %>%
  column_to_rownames("Idx")
clinical_cptac_hnscc  <- read.table("data/CPTAC_HNSCC/HS_CPTAC_HNSCC_CLI.tsi",
                             sep="\t", header = T, stringsAsFactors = F)
clinical_cptac_hnscc  <- clinical_cptac_hnscc [-1,] %>% 
  mutate(case_id=str_replace_all(case_id, "-", "\\."),
         overall_survival=as.numeric(overall_survival),
         overall_free_status=as.numeric(overall_free_status),
         progression_free_status=as.numeric(progression_free_status),
         progression_free_survival=as.numeric(progression_free_survival))


# TCGA-ACC
gene_exp_tcga_acc <- read.table("data/TCGA_ACC/ACC.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                                sep="\t", stringsAsFactors = F, header = T)
gene_exp_tcga_acc <- gene_exp_tcga_acc[-1,] %>% 
  mutate(Hybridization.REF=str_replace(Hybridization.REF, "\\|\\d+$", "")) %>% 
  filter(Hybridization.REF!="?" & Hybridization.REF!= "SLC35E2") %>% 
  column_to_rownames("Hybridization.REF") %>% 
  mutate_if(is.character, as.numeric)
colnames(gene_exp_tcga_acc) <- colnames(gene_exp_tcga_acc) %>% 
  str_extract("TCGA\\.[\\w\\d]{2}\\.[\\w\\d]{4}")
clinical_cptac_acc  <- read.table("data/TCGA_ACC/HS_CPTAC_ACC_CLI.tsi",
                                    sep="\t", header = T, stringsAsFactors = F) %>% 
  column_to_rownames("attrib_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("case_id")
clinical_cptac_acc  <- clinical_cptac_acc %>% 
  mutate(overall_survival=as.numeric(overall_survival),
         overall_free_status=as.numeric(status))

# TCGA-UVM
gene_exp_tcga_uvm <- read.table("data/TCGA_UVM/UVM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                                sep="\t", stringsAsFactors = F, header = T)
gene_exp_tcga_uvm <- gene_exp_tcga_uvm[-1,] %>% 
  mutate(Hybridization.REF=str_replace(Hybridization.REF, "\\|\\d+$", "")) %>% 
  filter(Hybridization.REF!="?" & Hybridization.REF!= "SLC35E2") %>% 
  column_to_rownames("Hybridization.REF") %>% 
  mutate_if(is.character, as.numeric)
colnames(gene_exp_tcga_uvm) <- colnames(gene_exp_tcga_uvm) %>% 
  str_extract("TCGA\\.[\\w\\d]{2}\\.[\\w\\d]{4}")
clinical_cptac_uvm  <- read.table("data/TCGA_UVM/HS_CPTAC_UVM_CLI.tsi",
                                  sep="\t", header = T, stringsAsFactors = F) %>% 
  column_to_rownames("attrib_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("case_id")
clinical_cptac_uvm  <- clinical_cptac_uvm %>% 
  mutate(overall_survival=as.numeric(overall_survival),
         overall_free_status=as.numeric(status))

# TCGA-LAML
gene_exp_tcga_laml <- read.table("data/TCGA_LAML/LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",
                                sep="\t", stringsAsFactors = F, header = T)
gene_exp_tcga_laml <- gene_exp_tcga_laml[-1,] %>% 
  mutate(Hybridization.REF=str_replace(Hybridization.REF, "\\|\\d+$", "")) %>% 
  filter(Hybridization.REF!="?" & Hybridization.REF!= "SLC35E2") %>% 
  column_to_rownames("Hybridization.REF") %>% 
  mutate_if(is.character, as.numeric)
colnames(gene_exp_tcga_laml) <- colnames(gene_exp_tcga_laml) %>% 
  str_extract("TCGA\\.[\\w\\d]{2}\\.[\\w\\d]{4}")
clinical_cptac_laml  <- read.table("data/TCGA_LAML/HS_CPTAC_LAML_CLI.tsi",
                                  sep="\t", header = T, stringsAsFactors = F) %>% 
  column_to_rownames("attrib_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("case_id")
clinical_cptac_laml  <- clinical_cptac_laml %>% 
  mutate(overall_survival=as.numeric(overall_survival),
         overall_free_status=as.numeric(status))

# CPTAC-COAD data:
gene_exp_cptac_coad <- read.table("data/CPTAC_COAD/Human__CPTAC_COAD__UNC__RNAseq__HiSeq_RNA__03_01_2017__BCM__Gene__BCM_RSEM_UpperQuartile_log2.cct",
                                   sep="\t", stringsAsFactors = F, header = T) %>%
  column_to_rownames("attrib_name") 
# CPTAC-COAD: There is no survival days in survival data

# TCGA-COAD data:
gene_exp_tcga_coad <- read.table(gzfile("data/TCGA_COAD/Human__TCGA_COADREAD__UNC__RNAseq__HiSeq_RNA__01_28_2016__BI__Gene__Firehose_RSEM_log2.cct.gz"), 
                                 sep="\t", stringsAsFactors = F, header = T) %>% 
  column_to_rownames("attrib_name")
clinical_tcga_coad  <- read.table("data/TCGA_COAD/HS_CPTAC_TCGA_COAD_CLI.tsi",
                                   sep="\t", header = T, stringsAsFactors = F) %>% 
  column_to_rownames("attrib_name") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("case_id")
clinical_tcga_coad  <- clinical_tcga_coad %>% 
  mutate(overall_survival=as.numeric(overall_survival),
         overall_free_status=as.numeric(status))

###############################################################################
# Process data
###############################################################################

# CPTAC-HNSC data:
gsva_metab_cptac_hnscc  <- run_gsva_metabolic(gene_exp_cptac_hnscc ,
                                       kcdf="Gaussian")
kmeans_clusters_cptac_hnscc <- kmeans_gsva_metabolic(gsva_metab_cptac_hnscc )
kmeans_clusters_cptac_hnscc$heatmap
kmeans_metab_clust_surv(kmeans_res=kmeans_clusters_cptac_hnscc $kmean_res,
                        clinical_data=clinical_cptac_hnscc,
                        sample_id_col="case_id",
                        surv_status_col="overall_free_status",
                        surv_time_col="overall_survival")

# TCGA-ACC data:
gsva_metab_acc <- run_gsva_metabolic(gene_exp_tcga_acc,
                                       kcdf="Gaussian")
kmeans_clusters_acc <- kmeans_gsva_metabolic(gsva_metab_acc)
kmeans_clusters_acc$heatmap
kmeans_metab_clust_surv(kmeans_res=kmeans_clusters_acc$kmean_res,
                        clinical_data=clinical_cptac_acc,
                        sample_id_col="case_id",
                        surv_status_col="overall_free_status",
                        surv_time_col="overall_survival")

# TCGA-UVM data:
gsva_metab_uvm <- run_gsva_metabolic(gene_exp_tcga_uvm,
                                     kcdf="Gaussian")
kmeans_clusters_uvm <- kmeans_gsva_metabolic(gsva_metab_uvm)
kmeans_clusters_uvm$heatmap
kmeans_metab_clust_surv(kmeans_res=kmeans_clusters_uvm$kmean_res,
                        clinical_data=clinical_cptac_uvm,
                        sample_id_col="case_id",
                        surv_status_col="overall_free_status",
                        surv_time_col="overall_survival")

# TCGA-LAML data:
gsva_metab_laml <- run_gsva_metabolic(gene_exp_tcga_laml,
                                     kcdf="Gaussian")
kmeans_clusters_laml <- kmeans_gsva_metabolic(gsva_metab_laml)
kmeans_clusters_laml$heatmap
kmeans_metab_clust_surv(kmeans_res=kmeans_clusters_laml$kmean_res,
                        clinical_data=clinical_cptac_laml,
                        sample_id_col="case_id",
                        surv_status_col="overall_free_status",
                        surv_time_col="overall_survival")

# TCGA-COAD data:
gsva_metab_cptac_coad <- run_gsva_metabolic(gene_exp_tcga_coad,
                                      kcdf="Gaussian")
kmeans_clusters_cptac_coad <- kmeans_gsva_metabolic(gsva_metab_cptac_coad)
kmeans_clusters_cptac_coad$heatmap
kmeans_metab_clust_surv(kmeans_res=kmeans_clusters_cptac_coad$kmean_res,
                        clinical_data=clinical_tcga_coad,
                        sample_id_col="case_id",
                        surv_status_col="overall_free_status",
                        surv_time_col="overall_survival")
