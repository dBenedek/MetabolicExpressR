#' Gene Set Variation Analysis (GSVA) on gene expression data of tumors with 
#' MSigDB Hallmarks 
#'
#' Performs Gene Set Variation Analysis (GSVA) on gene expression matrix of 
#' tumor samples using MSigDB Hallmarks gene set collection
#' @param gene_exp_data N x M gene expression matrix genes as 
#' row names and samples as column names 
#' @param kcdf Character, either "Gaussian" or "Poisson" 
#' (see GSVA documentation)
#' @param n_cores Number of cores used for GSVA calculation (used by BiocParallel::MulticoreParam)
#' @return A K x M matrix with enrichment scores of samples per Hallmark gene set
#' @examples 
#' run_gsva_hallmarks(gene_exp_matrix, "Gaussian", 4);
#' @import GSVA
#' @import BiocParallel
#' @import msigdbr
#' @export
run_gsva_hallmarks <- function(gene_exp_data,
                               kcdf,
                               n_cores){
  
  # Performs GSVA with the MSigDB gene sets on the normalized, filtered
  # gene expression data, returns a data frame with ES scores per hallmark
  
  if(!is.data.frame(gene_exp_data)) stop("gene_exp_data must be a dataframe")
  if(!is.character(kcdf)) stop('kcdf must be a character ("Gaussian" or "Poisson")')
  if(!kcdf %in% c("Gaussian", "Poisson")) stop('kcdf should be either "Gaussian" or "Poisson"')
  if(!is.numeric(n_cores)) stop('n_cores must be numeric')
  
  # Check GSVA version:
  gsva_version <- packageVersion("GSVA")
  
  # Filter rows with 0 variance:
  genes_var <- apply(gene_exp_data, 1, var)
  gene_exp_data <- as.matrix(round(gene_exp_data[genes_var > 0,]))
  
  # Load gene set data:
  hum_hallmark <- msigdbr(species = "Homo sapiens",
                          category = "H") 
  gs = split(x = hum_hallmark$gene_symbol, 
             f = hum_hallmark$gs_name)
  
  # Run GSVA:
  multicoreParam <- BiocParallel::MulticoreParam(workers = n_cores)
  if (gsva_version >= "1.52.0"){
    gsvaPar <- GSVA::gsvaParam(exprData=gene_exp_data, 
                               geneSets=gs,
                               kcdf=kcdf,
                               minSize=9,
                               maxSize=300)
    gsva_es <- GSVA::gsva(param=gsvaPar,
                          BPPARAM=multicoreParam,
                          verbose=FALSE)
  } else {
    gsva_es <- GSVA::gsva(expr=gene_exp_data, 
                          gset.idx.list=gs, 
                          verbose=FALSE, 
                          kcdf=kcdf, 
                          method="gsva", 
                          min.sz=9, 
                          max.sz=300,
                          parallel.sz=n_cores)
  }
  
  return(gsva_es)
}