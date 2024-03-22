#' Runs PROGENy on the input gene expression data 
#'
#' Runs PROGENy on the input gene expression data and compares PROGENy pathway
#' activtiy between the identified clusters. Wilcoxon Rank Sum test is used if 
#' number of clusters = 2, if number of clusters > 2 then Kruskal-Wallis Rank 
#' Sum test is used.
#' @param gene_exp_data N x M gene expression matrix genes as 
#' row names and samples as column names 
#' @param kmeans_res Data frame holding the cluster membership of each sample
#' @param vst Logical, whether to do VST on the gene expression data or not
#' (default: FALSE)
#' @param top Numeric, he top n genes for generating the model matrix according 
#' to significance (p-value) (see PROGENy documentation for details)
#' @return A data frame with PROGENy pathway activity scores per sample, and
#' a boxplot showing the differential pathway testing results between clusters 
#' @examples 
#' run_progeny(gene_exp_data, kmeans_res, vst = FALSE, top = 100);
#' @import DESeq2
#' @import progeny
#' @import tibble
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import stats 
#' @export
run_progeny <- function(gene_exp_data,
                        kmeans_res,
                        vst = FALSE,
                        top = 100){
  # Runs PROGENy on the input gene expression data and compares pathway activity
  # between the identified clusters
  
  if(!is.data.frame(gene_exp_data)) stop("gene_exp_data must be a dataframe")
  if(!is.data.frame(kmeans_res)) stop("kmeans_res must be a dataframe")
  if(!vst %in% c(TRUE, FALSE)) stop("vst must be either FALSE (default) or TRUE")
  if(!is.null(top)){
    if(!is.numeric(top)) stop("top should be numeric specifying the top n genes for generating the model matrix according to significance")
  }
  
  # Do VST, if specified:
  if (isTRUE(vst)){
    gene_exp_data <- DESeq2::vst(as.matrix(round(gene_exp_cptac_hnscc)))
  }
  
  # Calculate PROGENy scores:
  progeny_scores_matrix <- progeny::progeny(as.matrix(gene_exp_data), 
                                            scale = TRUE,
                                            top = top) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("sample_ID") %>% 
    tidyr::pivot_longer(-sample_ID, names_to="pathway",
                        values_to = "score") %>% 
    dplyr::left_join(kmeans_res,
                     by="sample_ID")
  
  # Do differential testing for pathway acrivity between clusters -
  # Wilcoxon test if clusters == 2, Kruskal-Wallis test, if k > 2:
  if (length(unique(kmeans_res$cluster)) == 2){
    testing_results <- progeny_scores_matrix %>% 
      dplyr::group_by(pathway) %>% 
      dplyr::do(w = stats::wilcox.test(score~cluster, data=., paired = FALSE)) %>% 
      dplyr::summarise(pathway, p.value = w$p.value) %>% 
      dplyr::mutate(p.adj=stats::p.adjust(p.value, method="BH")) %>% 
      dplyr::arrange(p.adj, p.value)
  } else if (length(unique(kmeans_res$cluster)) > 2){
    testing_results <- progeny_scores_matrix %>% 
      dplyr::group_by(pathway) %>% 
      dplyr::do(k = stats::kruskal.test(score~cluster, data=.)) %>% 
      dplyr::summarise(pathway, p.value = k$p.value) %>% 
      dplyr::mutate(p.adj=stats::p.adjust(p.value, method="BH")) %>% 
      dplyr::arrange(p.adj, p.value)
  } else {stop("For running differential testing on data, at least 2 clusters should be provided.")}
  
  # Prepare data for plotting:
  plot_data <- progeny_scores_matrix %>% 
    dplyr::left_join(testing_results, 
                     by="pathway") %>% 
    dplyr::mutate(label=paste0(pathway, ", P.adj=", round(p.adj, 3)),
                  label=ifelse(grepl("adj=0$", label),
                               stringr::str_replace(label, "adj=0$",
                                                    "adj<0.001"),
                               label))
  
  # Plot data:
  progeny_cluster_plot <- ggplot2::ggplot(plot_data, 
                                          aes(x=as.factor(cluster), 
                                              y=score))+ 
    ggplot2::geom_boxplot(outlier.shape = NA, width=.3, alpha=1)+
    ggplot2::geom_jitter(width = 0.25, alpha=.5)+
    facet_wrap(ggplot2::vars(label), scales = "free")+
    ggplot2::theme_bw()+
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(size=16),
          axis.text = ggplot2::element_text(size=14),
          axis.title = ggplot2::element_text(size=16))+
    ggplot2::xlab("Cluster")+
    ggplot2::ylab("PROGENy score")
  
  # Return data:
  return(list(progeny_res=progeny_scores_matrix,
              progeny_boxplot=progeny_cluster_plot))
}