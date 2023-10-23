#' Plot Kaplan-Meier analysis of clusters for the selected endpoint
#'
#' Creates a Kaplan-Meier (KM) plot of clusters for the selected clinical 
#' endpoint
#' @param kmeans_res Data frame holding the cluster membership of each sample
#' @param clinical_data Data frame with survival data of the cases analyzed
#' @param sample_id_col Column name of sample IDs in the clinical_data table
#' @param surv_status_col Column name of survival status (such as overall 
#' survival status) in the clinical_data table
#' @param surv_time_col Column name of survival time (such as overall survival time) 
#' in the clinical_data table
#' @return A Kaplan-Meier plot of the clusters for the selected endpoint
#' @examples 
#' kmeans_metab_clust_surv(kmeans_results, clinical_table, "case_id", 
#' "overall_surv_status", "overall_surv_time");
#' @export
kmeans_metab_clust_surv <- function(kmeans_res, 
                                    clinical_data, 
                                    sample_id_col, 
                                    surv_status_col, 
                                    surv_time_col) {
  
  # This function looks for association between the identified metabolic clusters and survival
  
  if(!is.data.frame(kmeans_res)) stop("kmeans_res must be a dataframe")
  if(!is.data.frame(clinical_data)) stop("clinical_data must be a dataframe")
  if(!is.character(clinical_data[[sample_id_col]])) stop("sample_id_col must be of character type")
  if(!is.numeric(clinical_data[[surv_status_col]])) stop("surv_status_col must be of character type")
  if(!is.numeric(clinical_data[[surv_time_col]])) stop("surv_time_col must be of character type")
  
  # Merge data:
  surv_data <-  kmeans_res %>% 
    left_join(clinical_data, by = c("sample_ID" = sample_id_col)) %>% 
    mutate(cluster = as.factor(cluster))
  
  # Run KM analysis:
  formula <- as.formula(paste("Surv(", surv_time_col, 
                              ",", surv_status_col, ") ~ cluster"))
  fit <- surv_fit(formula, data = surv_data)
  
  # Plot:
  survplot <- ggsurvplot(
    fit,
    data = surv_data,  
    pval = TRUE,
    conf.int = TRUE,
    title = str_replace_all(eval(surv_time_col),
                            "_", " ") %>% 
      str_to_title,
    conf.int.style = "step",
    xlab = "Time in days",
    ggtheme = theme_light(),
    risk.table = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.y.text = TRUE,
    ncensor.plot = FALSE,
    surv.median.line = "hv"
  )

  return(survplot)
}
