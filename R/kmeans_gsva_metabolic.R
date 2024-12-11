#' K-means clustering on GSVA matrix of KEGG metabolic pathway scores 
#' 
#' Performs k-means clustering with optimal k or user-defined k on the GSVA 
#' matrix of KEGG metabolic pathway enrichment scores of tumors
#' @param gsva_data N x M matrix with N pathways and M tumor samples
#' @param kegg_gs Named list with K elements of KEGG metabolic pathways gene set
#' @param user_def_k Logical, whether to use a user-defined k for k-means 
#' clustering (default: FALSE)
#' @param k Numeric, user-defined k for k-means clustering (default: NULL)
#' @return A data frame with cluster membership of each sample, and a heatmap 
#' plot with GSVA enrichment scores per tumor and pathway (samples are 
#' ordered by cluster)
#' @examples 
#' kmeans_gsva_metabolic(gsva_matrix, kegg_gs);
#' kmeans_gsva_metabolic(gsva_matrix, kegg_gs, user_def_k=TRUE, k=3);
#' @import NbClust
#' @import dplyr
#' @import stats
#' @import ComplexHeatmap
#' @import tibble
#' @import viridis
#' @export
kmeans_gsva_metabolic <- function(gsva_data,
                                  kegg_gs,
                                  user_def_k=FALSE,
                                  k=NULL){
  
  # Runs K-means clustering with the optimal K (or user defined k)
  # on the GSVA KEGG metabolic data
  # and assigns cluster membership to samples
  
  # Check arguments:
  if(!is.matrix(gsva_data)) stop("gsva_data must be an N x M matrix with N pathways and M samples")
  if(!is.list(kegg_gs)) stop("kegg_gs must be a named list (gene set collection) with N elements (pathways)")
  if(!user_def_k %in% c(TRUE, FALSE)) stop("user_def_k must be either FALSE (default) or TRUE")
  if(!is.null(k) & isFALSE(user_def_k)) message("Warning: user_def_k = FALSE, so optimal k will be used")
  if(isTRUE(user_def_k) & is.null(k)) stop("if user_def_k = TRUE then you should specify the parameter k")
  if(!is.null(k)){
    if(!is.numeric(k) | k > ncol(gsva_data)) stop("k should be numeric and smaller than the number of samples")
  }
  
  if (isFALSE(user_def_k)){
  
  # Define optimal number of clusters:
  selected <- c("kl", "ch", "hartigan", "db", 
                "silhouette", "duda", "pseudot2", "beale", 
                "ratkowsky", "ball", "ptbiserial", "gap", 
                "frey", "mcclain", "gamma", "gplus", 
                "tau", "dunn", "sdindex", "sdbw")
  results <- vector("list", length(selected))
  
  for (i in 1:length(selected)) {
    possibleError <- tryCatch(
      {
      results[[i]] <- NbClust::NbClust(t(gsva_data), 
                                       min.nc=2, 
                                       max.nc=10, 
                                       method="kmeans", 
                                       index=selected[i])},
      error = function(e){
        message("An error occured: ", e$message)
      },
      warning = function(w){
        message("A warning occured: ", w$message)
      },
      finally = {
        message("Index tested: ", selected[i])
      })}
  
  # Filter out failed indices:
  selected <- selected[unlist(lapply(1:length(results), 
                                     function(x) "Best.nc" %in% names(results[[x]])))]
  results <- results[unlist(lapply(1:length(results), 
                    function(x) "Best.nc" %in% names(results[[x]])))]
  
  best_nc <- lapply(results, function(x) x$Best.nc)
  names(best_nc) <- selected
  
  # Optimal number of k:
  k <- best_nc[!sapply(best_nc,is.null)] %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    dplyr::arrange(Number_clusters) %>% 
    dplyr::group_by(Number_clusters) %>% 
    dplyr::summarize(count=n()) %>% 
    dplyr::arrange(desc(count)) %>% 
    head(1) %>% 
    dplyr::pull(Number_clusters)
  }
  
  # Perform K-means clustering:
  set.seed(123)
  km_res <- stats::kmeans(t(gsva_data), 
                   k, 
                   nstart = 100,
                   iter.max = 100)
  
  # Plot heatmap:
  pathways <- data.frame(name=names(kegg_gs)) %>% 
    dplyr::mutate(category=str_replace(as.character(name), "^\\d+\\s", ""),
           category=str_replace(as.character(category), "\\s-\\s\\d+.+", "")) %>% 
    dplyr::mutate(pathway=str_extract(name, "-\\s\\d+\\s.+"),
           pathway=str_replace(pathway, "-\\s\\d+\\s", ""),
           pathway=str_replace(pathway, "\\s\\[PATH.*", ""),
           pathway=stringr::str_to_title(pathway)) 
  
  plot_data <- gsva_data %>% 
    scale() %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("pathway") %>% 
    dplyr::mutate(pathway=str_extract(pathway, "-\\s\\d+\\s.+\\[PATH"),
           pathway=str_replace(pathway, "-\\s\\d+\\s", ""),
           pathway=str_replace(pathway, "\\s\\[PATH", ""),
           pathway=stringr::str_to_title(pathway)) %>% 
    dplyr::left_join(pathways, by="pathway") %>% 
    dplyr::arrange(category, desc(pathway)) %>% 
    dplyr::select(-category) %>% 
    tibble::column_to_rownames("pathway")
  
  plot_annot <- as.data.frame(km_res$cluster) %>% 
    dplyr::rename("cluster"="km_res$cluster") %>% 
    tibble::rownames_to_column("sample") %>% 
    dplyr::arrange(cluster) %>% 
    dplyr::mutate(cluster=factor(cluster)) %>% 
    tibble::column_to_rownames("sample")
  
  plot_data <- plot_data %>% 
    dplyr::select(rownames(plot_annot))
  
  row_ha = ComplexHeatmap::rowAnnotation(df = data.frame(rownames(plot_data)) %>% 
                           dplyr::left_join(pathways,
                                     by=c("rownames.plot_data."="pathway")) %>% 
                           dplyr::select(category),
                         col = list(category = c("Amino acid metabolism"="#EF5350",
                                                 "Biosynthesis of other secondary metabolites"="#F48FB1",
                                                 "Carbohydrate metabolism"="#AB47BC",
                                                 "Energy metabolism"="#5C6BC0",
                                                 "Glycan biosynthesis and metabolism"="#29B6F6",  
                                                 "Lipid metabolism"="#26A69A",
                                                 "Metabolism of cofactors and vitamins"="#9CCC65",
                                                 "Metabolism of other amino acids"="#FFEE58",
                                                 "Metabolism of terpenoids and polyketides"="#FFA726",    
                                                 "Not included in regular maps"="#8D6E63", 
                                                 "Nucleotide metabolism"="#78909C",
                                                 "Xenobiotics biodegradation and metabolism"="#ECEFF1")),
                         show_annotation_name = F)
  
  column_ha = ComplexHeatmap::HeatmapAnnotation(cluster = plot_annot$cluster, 
                                show_annotation_name = T,
                                annotation_name_side = "left")
  
  heatmap <- ComplexHeatmap::Heatmap(plot_data,
                     cluster_rows = F,
                     cluster_columns = F,
                     col = rev(viridis::viridis(n=100, option="magma")),
                     show_column_names = F,
                     name="GSVA score",
                     row_names_max_width = ComplexHeatmap::max_text_width(
                       rownames(plot_data), 
                       gp = grid::gpar(fontsize = 12)
                       ),
                     column_split = plot_annot$cluster,
                     column_order = colnames(plot_data),
                     column_gap = unit(2, "mm"),
                     left_annotation = row_ha,
                     row_names_side = "left",
                     heatmap_legend_param = list(
                       ncol=1),
                     width = unit(12, "cm"),
                     border = TRUE,
                     row_split = data.frame(rownames(plot_data)) %>% 
                       dplyr::left_join(pathways,
                                 by=c("rownames.plot_data."="pathway")) %>% 
                       dplyr::select(category),
                     show_row_names = T,
                     show_row_dend = F,
                     row_gap = unit(2, "mm"),
                     row_title = NULL
  )
  
  return(list(kmean_res=data.frame(sample_ID=names(km_res$cluster),
                         cluster=km_res$cluster), 
              heatmap=heatmap))
}