#' Gene Set Variation Analysis (GSVA) on gene expression data of tumors with 
#' KEGG metabolic pathways
#'
#' Performs Gene Set Variation Analysis (GSVA) on gene expression matrix of 
#' tumor samples using KEGG metabolic pathways gene set collection
#' @param gene_exp_data N x M gene expression matrix genes as 
#' row names and samples as column names 
#' @param kcdf Character, either "Gaussian" or "Poisson" 
#' (see GSVA documentation)
#' @param kegg_gs Named list with K elements of KEGG metabolic pathways gene set
#' @param n_cores Number of cores used for GSVA calculation (used by BiocParallel::MulticoreParam)
#' @return A K x M matrix with enrichment scores of samples per KEGG 
#' metabolic pathway
#' @examples 
#' run_gsva_metabolic(gene_exp_matrix, "Gaussian", kegg_gs);
#' @import GSVA
#' @import BiocParallel
#' @export
run_gsva_metabolic <- function(gene_exp_data,
                               kcdf,
                               kegg_gs,
                               n_cores){
  
  # Performs GSVA with the KEGG metabolic gene sets on the normalized, filtered
  # gene expression data, returns a data frame with ES scores per pathway
  
  if(!is.data.frame(gene_exp_data)) stop("gene_exp_data must be a dataframe")
  if(!is.character(kcdf)) stop('kcdf must be a character ("Gaussian" or "Poisson")')
  if(!kcdf %in% c("Gaussian", "Poisson")) stop('kcdf should be either "Gaussian" or "Poisson"')
  if(!is.list(kegg_gs)) stop("kegg_gs must be a named list (gene set collection) with N elements (pathways)")
  if(!is.numeric(n_cores)) stop('n_cores must be numeric')
  
  # Check GSVA version:
  gsva_version <- packageVersion("GSVA")
  
  # Filter rows with 0 variance:
  genes_var <- apply(gene_exp_data, 1, var)
  gene_exp_data <- as.matrix(round(gene_exp_data[genes_var > 0,]))
  
  # Run GSVA:
  multicoreParam <- BiocParallel::MulticoreParam(workers = n_cores)
  if (gsva_version >= "1.52.0"){
    gsvaPar <- GSVA::gsvaParam(exprData=gene_exp_data, 
                               geneSets=kegg_gs,
                               kcdf=kcdf,
                               minSize=9,
                               maxSize=300)
    gsva_es <- GSVA::gsva(param=gsvaPar,
                          BPPARAM=multicoreParam,
                          verbose=FALSE)
  } else {
    gsva_es <- GSVA::gsva(expr=gene_exp_data, 
                          gset.idx.list=kegg_gs, 
                          verbose=FALSE, 
                          kcdf=kcdf, 
                          method="gsva", 
                          min.sz=9, 
                          max.sz=300,
                          parallel.sz=n_cores)
  }
  
  return(gsva_es)
}