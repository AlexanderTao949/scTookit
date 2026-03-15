#' Determine Significant Principal Components
#'
#' Calculates the number of PCs needed to explain 90% of cumulative variance.
#'
#' @param seu Seurat object
#' @param embedding Name of reduction embedding to use (default: "pca")
#'
#' @return Integer indicating the number of PCs to retain
#'
#' @examples
#' # Using default PCA embedding:
#' pcs <- pcs_determine(seu)
#'
#' # Using a custom embedding:
#' pcs <- pcs_determine(seu, embedding = "custom_pca")
#'
#' @export
pcs_determine <- function(seu, embedding = "pca"){
  xx <- cumsum(seu[[embedding]]@stdev^2)
  xx <- xx / max(xx)
  ndim <- which(xx > 0.9)[1]
  return(ndim)
}

#' Calculate Silhouette Widths for Cluster Validation
#'
#' Computes silhouette widths to evaluate cluster separation quality at a given resolution.
#' Returns per-cell metrics and cluster-level summary statistics.
#'
#' @param resolution_col Character name of the metadata column containing cluster assignments
#' @param seu Seurat object with dimensionality reduction (default: harmony embeddings)
#'
#' @return Nested list containing:
#'   - `per_cell`: Dataframe with cell-level silhouette metrics
#'   - `per_cluster`: Dataframe with cluster-level summary statistics
#'
#' @importFrom bluster approxSilhouette
#' @importFrom dplyr group_by summarise mutate %>%
#' @importFrom stats quantile
#' @export
#'
#' @examples
#' sil_results <- caculate_silhouette_width("RNA_snn_res.0.8", seu_object)
caculate_silhouette_width <- function(resolution_col, seu) {
  embeddings <- seu@reductions$harmony@cell.embeddings
  clusters <- seu@meta.data[[resolution_col]]

  si <- bluster::approxSilhouette(embeddings, clusters)

  result <- data.frame(
    cell = rownames(si),
    cluster = as.character(si$cluster),
    other = as.character(si$other),
    silhouette_width = si$width,
    resolution = resolution_col
  )

  cluster_stats <- result %>%
    group_by(cluster) %>%
    summarise(
      mean_silhouette = mean(silhouette_width, na.rm = TRUE),
      median_silhouette = median(silhouette_width, na.rm = TRUE),
      min_silhouette = min(silhouette_width, na.rm = TRUE),
      q25 = quantile(silhouette_width, 0.25, na.rm = TRUE),
      q75 = quantile(silhouette_width, 0.75, na.rm = TRUE),
      n_cells = n(),
      prop_negative = sum(silhouette_width < 0) / n()
    ) %>%
    mutate(resolution = resolution_col)

  return(list(
    per_cell = result,
    per_cluster = cluster_stats
  ))
}

#' Visualize Silhouette Widths for Cluster Validation
#'
#' Creates a boxplot of silhouette widths per cluster to evaluate clustering quality at a given resolution.
#'
#' @param silhouette_res_list List object containing silhouette results (output from `caculate_silhouette_width`)
#' @param resolution_col Character string specifying the resolution column name to visualize
#'
#' @return A ggplot object showing silhouette width distributions by cluster
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_hline labs theme_classic theme
#' @importFrom ggplot2 element_text element_line element_blank element_rect
#' @export
#'
#' @examples
#' plot <- silhouette_width_plot(sil_results, "RNA_snn_res.0.8")
#' print(plot)
silhouette_width_plot <- function(silhouette_res_list, resolution_col){
  plot_data <- silhouette_res_list[[resolution_col]][["per_cell"]]
  p <- ggplot(plot_data, aes(x = cluster, y = silhouette_width)) +
    geom_boxplot(color = "black",
                 fill = NA,
                 width = 0.7,
                 alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.5) +
    labs(
      title = resolution_col,
      x = "Cluster",
      y = "Silhouette Width"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white", color = NA)
    )
}
