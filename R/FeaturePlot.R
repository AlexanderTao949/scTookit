#' Threshold-based Feature Visualization for Seurat Objects
#'
#' Visualizes single-cell data by highlighting cells with feature expression above a specified threshold
#' on a dimensional reduction plot. Cells below the threshold are colored in grey.
#'
#' @param seu A Seurat object containing single-cell data and dimensional reductions.
#' @param feature Character vector of feature names (genes or metadata columns) to visualize.
#' @param threshold Numeric threshold value for highlighting cells. If NULL (default),
#'   automatically calculated as mean + 2SD of the feature values.
#' @param reduction Character string specifying the dimensional reduction to use (e.g., "umap", "tsne", "pca").
#' @param dims Integer vector of length 2 specifying which dimensions to plot. Default is c(1, 2).
#' @param color Character vector of colors for the gradient scale. Default uses RColorBrewer "Reds" palette.
#' @param alpha Numeric value between 0 and 1 specifying point transparency. Default is 0.6.
#' @param dot.size Numeric value specifying point size. Default is 0.5.
#' @param title Character string for plot title. If NULL, uses the feature name.
#' @param title.size Numeric value for title font size. Default is 8.
#' @param legend.title Character string for legend title. If NULL, no title is shown.
#' @param legend.position Character string specifying legend position. Default is "bottom".
#' @param legend.title.size Numeric value for legend title font size. Default is 8.
#' @param legend.text.size Numeric value for legend text font size. Default is 7.
#' @param ncol Integer specifying number of columns for arranging multiple plots. Default is 1.
#'
#' @returns A patchwork object containing the arranged ggplot(s). If multiple features are provided,
#'   plots are arranged in a grid according to ncol parameter with shared legend.
#'
#' @importFrom Seurat Embeddings FetchData
#' @importFrom ggplot2 ggplot geom_point aes labs theme_void theme element_text scale_color_gradientn margin
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom scales squish
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with automatic threshold (mean + 2SD)
#' ThreshPlot(seurat_obj, feature = "CD3E", reduction = "umap")
#'
#' # Multiple features with custom threshold
#' ThreshPlot(seurat_obj, feature = c("CD3E", "CD4", "CD8A"),
#'            threshold = 2, reduction = "umap", ncol = 3)
#' }
ThreshPlot <- function(seu,
                       feature,
                       threshold = NULL,
                       reduction = "umap",
                       dims = c(1, 2),
                       color = RColorBrewer::brewer.pal(9, "Reds")[c(2,4,6,7,8,9)],
                       alpha = 0.6,
                       dot.size = 0.5,
                       title = NULL,
                       title.size = 8,
                       legend.title = NULL,
                       legend.position = "bottom",
                       legend.title.size = 8,
                       legend.text.size = 7,
                       ncol = 1){

  plist <- lapply(feature, function(x){

    plotdata <- cbind(seu@meta.data,
                      Embeddings(seu, reduction),
                      FetchData(seu, var = x))

    mean.val <- mean(plotdata[, x])
    sd.val <- sd(plotdata[, x])
    thresh <- if (is.null(threshold)) round(mean.val + 2*sd.val, 2) else threshold
    plotdata$color_val <- ifelse(plotdata[[x]] > thresh, plotdata[[x]], NA)
    plotdata <- plotdata[order(plotdata[[x]], decreasing = FALSE), ]

    x_col <- colnames(Embeddings(seu, reduction))[dims[1]]
    y_col <- colnames(Embeddings(seu, reduction))[dims[2]]

    p <- ggplot(plotdata) + geom_point(
      aes(x = .data[[x_col]],
          y = .data[[y_col]],
          color = color_val),
      size = dot.size,
      stroke = 0,
      shape = 16,
      alpha = alpha
    ) +
      labs(color = legend.title)+
      theme_void() +
      theme(
        aspect.ratio = 1,
        text = element_text(family = "Arial")
      )+
      scale_color_gradientn(colours = color,
                            na.value = "grey80",
                            limits = c(thresh, max(plotdata[[x]], na.rm = TRUE)),
                            oob = scales::squish)
      if(is.null(title)){
        p <- p+labs(title = x)
      }else{
        p <- p+labs(title = title)
      }
      p <- p+theme(
        plot.title = element_text(
          size = title.size,
          face = "bold",
          family = "Arial",
          color = "black",
          hjust = 0.5,
          vjust = 1
        ))

  })

  p <- plist %>%
    patchwork::wrap_plots(ncol = ncol)+
    patchwork::plot_layout(guides = 'collect')&
    ggplot2::theme(
      legend.position = legend.position,
      legend.title = element_text(size = legend.title.size,
                                  face = "bold",
                                  family = "Arial",
                                  color = "black"),
      legend.text = element_text(size = legend.text.size),
      legend.key.size = grid::unit(0.15, "inch"),
      legend.box = "horizontal",
      legend.box.margin = margin(10, 0, 0, 0)
    )


  return(p)

}
