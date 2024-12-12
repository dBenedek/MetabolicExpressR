#' Plot mean hallmark enrichment scores per cluster and hallmark
#'
#' Creates a barplot of mean GSVA enrichment scores per cluster and hallmark 
#' (hallmarks are grouped by main category)
#' @param gsva_data N x M matrix with N pathways and M tumor samples
#' @param kmeans_res Data frame holding the cluster membership of each sample
#' @param hallmarks_data Data frame with category information of the Hallmarks
#' gene set collection
#' @return A barplot of mean enrichment scores
#' @examples 
#' plot_mean_hallmark_activity(gsva_matrix, kmeans_res);
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import ggplot2
#' @export
plot_mean_hallmark_activity <- function(gsva_data,
                                        kmeans_res,
                                        hallmarks_data){
  # This function plots the mean GSVA activity per hallmark and cluster
  
  if(!is.data.frame(kmeans_res)) stop("kmeans_res must be a dataframe")
  if(!is.matrix(gsva_data)) stop("gsva_data must be an N x M matrix with N pathways and M samples")
  if(!is.data.frame(hallmarks_data)) stop("hallmarks_data must be a dataframe")
  
  # prepare plot data:
  plot_data <- gsva_data |> 
    t() |> 
    as.data.frame() |> 
    tibble::rownames_to_column("sample_ID") |> 
    tidyr::pivot_longer(-sample_ID,
                        names_to="pathway",
                        values_to="ES") |> 
    dplyr::left_join(kmeans_res,
                     by="sample_ID") |> 
    dplyr::select(-sample_ID) |> 
    dplyr::group_by(pathway, cluster) |> 
    dplyr::summarise(across(everything(), mean, na.rm = TRUE)) |> 
    dplyr::left_join(hallmarks_data,
                     by=c("pathway"="Name")) |> 
    dplyr::arrange(factor(pathway, levels = rev(hallmarks_data$Name))) |> 
    dplyr::mutate(pathway=str_replace(pathway, "HALLMARK_", ""),
                  pathway=str_replace_all(pathway, "_", " "),
                  cluster=as.factor(cluster)) 
  
  # create barplot:
  barplot <- ggplot2::ggplot(plot_data, aes(x=pathway, y=ES, fill=cluster))+
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
    xlab("Hallmark")+
    ylab("Mean ES")+
    facet_grid(rows=vars(Category), scales="free", space = "free_y")+
    theme(
      strip.background = element_rect(),
      strip.text = element_text(size=7),
      panel.spacing = unit(2, "mm")
    )
  
  return(barplot)
}