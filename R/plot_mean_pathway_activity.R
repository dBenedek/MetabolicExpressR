#' Plot mean pathway enrichment scores per cluster and pathway
#'
#' Creates a barplot of mean GSVA enrichment scores per cluster and pathway 
#' (pathways are grouped by main KEGG metabolic category)
#' @param gsva_data N x M matrix with N pathways and M tumor samples
#' @param kegg_gs N long named list of KEGG metabolic pathways gene set
#' @param kmeans_res Data frame holding the cluster membership of each sample
#' @return A barplot of mean enrichment scores
#' @examples 
#' plot_mean_pathway_activity(gsva_matrix, kegg_gs, kmeans_res);
#' @export
plot_mean_pathway_activity <- function(gsva_data,
                                       kegg_gs,
                                       kmeans_res){
  # This function plots the mean GSVA activity per pathway and cluster
  
  if(!is.data.frame(kmeans_res)) stop("kmeans_res must be a dataframe")
  if(!is.matrix(gsva_data)) stop("gsva_data must be an N x M matrix with N pathways and M samples")
  if(!is.list(kegg_gs)) stop("kegg_gs must be a named list (gene set collection) with N elements (pathways)")
     
  # pathways category data:
  pathways <- data.frame(name=names(kegg_gs)) %>% 
    mutate(category=str_replace(as.character(name), "^\\d+\\s", ""),
           category=str_replace(as.character(category), "\\s-\\s\\d+.+", "")) %>% 
    mutate(pathway=str_extract(name, "-\\s\\d+\\s.+"),
           pathway=str_replace(pathway, "-\\s\\d+\\s", ""),
           pathway=str_replace(pathway, "\\s\\[PATH.*", "")) %>% 
    arrange(category)
  
  # prepare plot data:
  plot_data <- gsva_data |> 
    t() |> 
    as.data.frame() |> 
    rownames_to_column("sample_ID") |> 
    pivot_longer(-sample_ID,
                 names_to="pathway",
                 values_to="ES") |> 
    left_join(kmeans_res,
              by="sample_ID") |> 
    dplyr::select(-sample_ID) |> 
    group_by(pathway, cluster) |> 
    summarise(across(everything(), mean, na.rm = TRUE)) |> 
    mutate(pathway=str_extract(pathway, "-\\s\\d+\\s.+"),
           pathway=str_replace(pathway, "-\\s\\d+\\s", ""),
           pathway=str_replace(pathway, "\\s\\[PATH.*", ""),
           cluster=as.factor(cluster)) |> 
    left_join(pathways,
              by="pathway") |> 
    arrange(factor(pathway, levels = rev(pathways$pathway)))
  
  # create barplot:
  barplot <- ggplot(plot_data, aes(x=pathway, y=ES, fill=cluster))+
    geom_col(position=position_dodge2(preserve = "single"), 
             width=.6)+
    coord_flip()+
    theme_bw()+
    theme(axis.text.y = element_text(size=12),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(size=14),
          legend.position = "right",
          axis.title.x = element_text(size=16),
          axis.title.y = element_blank(),
          legend.title = element_text(size=14,face="bold"),
          legend.text = element_text(size=14),
          plot.title = element_text(size=18, hjust=.5))+
    ggtitle("Mean ES scores per cluster")+
    xlab("Metabolic pathway")+
    ylab("Mean ES")+
    facet_grid(rows=vars(category), scales="free", space = "free_y")+
    theme(
      strip.background = element_rect(),
      strip.text = element_text(size=0),
      panel.spacing = unit(2, "mm")
    )
  
  return(barplot)
}