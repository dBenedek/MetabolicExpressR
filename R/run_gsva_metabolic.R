run_gsva_metabolic <- function(gene_exp_data,
                               kcdf,
                               kegg_gs){
  library("GSVA")
  library("parallel")
  
  # Performs GSVA with the KEGG metabolic gene sets on the normalized, filtered
  # gene expression data, returns a data frame with ES scores per pathway
  
  if(!is.data.frame(gene_exp_data)) stop("gene_exp_data must be a dataframe")
  if(!is.character(kcdf)) stop('kcdf must be a character ("Gaussian" or "Poisson")')
  if(!is.list(kegg_gs)) stop("kegg_gs must be a named list (gene set collection) with N elements (pathways)")
  
  # Filter rows with 0 variance:
  genes_var <- apply(gene_exp_data, 1, var)
  gene_exp_data <- as.matrix(round(gene_exp_data[genes_var > 0,]))
  
  # Run GSVA:
  gsva_es <- gsva(expr=gene_exp_data, 
                  gset.idx.list=kegg_gs, 
                  verbose=FALSE, 
                  kcdf=kcdf, 
                  method="gsva", 
                  min.sz=9, 
                  max.sz=300,
                  parallel.sz=1L)
  
  return(gsva_es)
}