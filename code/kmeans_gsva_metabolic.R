kmeans_gsva_metabolic <- function(gsva_data){
  
  library("tidyverse")
  library("NbClust")
  library("pheatmap")
  library("fgsea")
  
  # Runs K-means clustering with the optimal K on the GSVA KEGG metabolic data
  # and assigns cluster membership to samples
  
  # Load gene set data:
  kegg_gs <- gmtPathways("/home/benedek_danko/R/genexp_metabolic_subtyping/data/kegg_metabolic_human_20211026.gmt")
  
  # Define optimal number of clusters:
  selected <- c("kl", "ch", "hartigan", "db", 
                "silhouette", "duda", "pseudot2", "beale", 
                "ratkowsky", "ball", "ptbiserial", "gap", 
                "frey", "mcclain", "gamma", "gplus", 
                "tau", "dunn", "hubert", "sdindex", 
                "dindex", "sdbw")
  results <- vector("list", length(selected))
  
  for (i in 1:length(selected)) {
    results[[i]] <- try(NbClust(t(gsva_data), 
                                min.nc=2, 
                                max.nc=10, 
                                method="kmeans", 
                                index=selected[i]))
  }
  
  best_nc <- lapply(results, function(x) x$Best.nc)
  names(best_nc) <- selected
  # Optimal number of k:
  optimal_k <- best_nc[!sapply(best_nc,is.null)] %>% 
    as.data.frame() %>% 
    t() %>% 
    as.data.frame() %>% 
    arrange(Number_clusters) %>% 
    group_by(Number_clusters) %>% 
    summarize(count=n()) %>% 
    arrange(desc(count)) %>% 
    head(1) %>% 
    pull(Number_clusters)
  
  # Perform K-means clustering:
  set.seed(123)
  km_res <- kmeans(t(gsva_data), 
                   optimal_k, 
                   nstart = 100,
                   iter.max = 100)
  
  # Plot heatmap:
  pathways <- data.frame(name=names(kegg_gs)) %>% 
    mutate(category=str_replace(as.character(name), "^\\d+\\s", ""),
           category=str_replace(as.character(category), "\\s-\\s\\d+.+", "")) %>% 
    mutate(pathway=str_extract(name, "-\\s\\d+\\s.+"),
           pathway=str_replace(pathway, "-\\s\\d+\\s", ""),
           pathway=str_replace(pathway, "\\s\\[PATH.*", "")) %>% 
    arrange(category)
  
  plot_data <- gsva_data %>% 
    scale() %>% 
    as.data.frame() %>% 
    rownames_to_column("pathway") %>% 
    mutate(pathway=str_extract(pathway, "-\\s\\d+\\s.+\\[PATH"),
           pathway=str_replace(pathway, "-\\s\\d+\\s", ""),
           pathway=str_replace(pathway, "\\s\\[PATH", "")) %>% 
    column_to_rownames("pathway")
  
  plot_annot <- as.data.frame(km_res$cluster) %>% 
    dplyr::rename("cluster"="km_res$cluster") %>% 
    rownames_to_column("sample") %>% 
    arrange(cluster) %>% 
    mutate(cluster=factor(cluster)) %>% 
    column_to_rownames("sample")
  
  plot_data <- plot_data %>% 
    dplyr::select(rownames(plot_annot))
  
  row_ha = rowAnnotation(df = data.frame(rownames(plot_data)) %>% 
                           left_join(pathways,
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
  
  column_ha = HeatmapAnnotation(cluster = plot_annot$cluster, 
                                show_annotation_name = T,
                                annotation_name_side = "left")
  
  heatmap <- Heatmap(plot_data,
                     cluster_rows = F,
                     cluster_columns = F,
                     col = rev(viridis(n=100, option="magma")),
                     show_column_names = F,
                     name="GSVA score",
                     row_names_max_width = max_text_width(
                       rownames(plot_data), 
                       gp = gpar(fontsize = 12)
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
                       left_join(pathways,
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