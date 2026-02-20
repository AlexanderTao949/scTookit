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

    meta <- seu@meta.data
    embed <- Embeddings(seu, reduction)
    feat_data <- FetchData(seu, var = x)
    
    meta <- meta[, !colnames(meta) %in% colnames(feat_data), drop = FALSE]
    plotdata <- cbind(meta, embed, feat_data)

    mean.val <- mean(plotdata[, x])
    sd.val <- sd(plotdata[, x])
    thresh <- if (is.null(threshold)) round(mean.val + 2*sd.val, 2) else threshold

    plotdata$color_val <- ifelse(plotdata[[x]] > thresh, plotdata[[x]], NA)
    plotdata <- plotdata[order(plotdata[[x]], decreasing = FALSE), ]

    x_col <- colnames(embed)[dims[1]]
    y_col <- colnames(embed)[dims[2]]

    p <- ggplot(plotdata) + 
      geom_point(
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
      scale_color_gradientn(
        colours = color,
        na.value = "grey80",
        limits = c(thresh, max(plotdata[[x]], na.rm = TRUE)),
        oob = scales::squish
      )
    
    if(is.null(title)){
      p <- p + labs(title = x)
    } else {
      p <- p + labs(title = title)
    }
    
    p <- p + theme(
      plot.title = element_text(
        size = title.size,
        face = "bold",
        family = "Arial",
        color = "black",
        hjust = 0.5,
        vjust = 1
      )
    )
    return(p)
  })

  p <- plist %>%
    patchwork::wrap_plots(ncol = ncol) +
    patchwork::plot_layout(guides = 'collect') &
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


ViolinPlot <- function(seu = NULL, 
                       metacol = NULL, 
                       cellcol = NULL, 
                       cellchoose = NULL, 
                       feature = NULL,
                       cell_filter_thres = 0){
  
  plotdata <- FetchData(seu, var = c(cellcol, metacol, feature))


  plotdata <- plotdata %>% dplyr::filter(major_celltype == cellchoose & TNF > 0)
  
  p <- ggplot(data = plotdata, aes(x = ERScore_group,y = TNF, color = ERScore_group)) + 
    geom_violin(trim = FALSE,position = position_dodge(0.9)) +  # size：加粗小提琴的边框,trim = FALSE：确保小提琴图的尾部不会被裁剪
    geom_signif(comparisons = list(c("High_ERScore", "Low_ERScore")), # 添加显著性，y_position控制标记的y轴位置，textsize控制文本大小，size控制线的粗细
                test = "t.test",colour = "black",
                test.args = list(var.equal = TRUE, alternative = "two.sided"),
                map_signif_level = TRUE, textsize = 6, tip_length = c(0, 0),
                y_position = max(plotdata[,"TNF"]),
                vjust=0
    )+  
    stat_summary(aes(group = ERScore_group,color = ERScore_group), position = position_dodge(0.9),
                 width = 0.2, fun = mean, geom = "errorbar", 
                 fun.min = function(x) quantile(x, 0.25),
                 fun.max = function(x) quantile(x, 0.75)) + 
    stat_summary( aes(group = ERScore_group,color = ERScore_group), position = position_dodge(0.9),
                  width = 0.2,  show.legend = FALSE,
                  fun = mean, geom = "crossbar") + 
    scale_color_manual(values = c("Low_ERScore"="#667a57","High_ERScore"="#dfadab")) + 
    scale_y_continuous(limits = c(0, max(plotdata[,"TNF"])*1.1)) + 
    labs(title = cellchoose,
         x = "",
         y = paste0("TNF expression levels")) + 
    theme_classic() + 
    theme (legend.position = "none",
           axis.line = element_line(color = "black"), 
           axis.title = element_text(size = 10,,face = "bold"),
           axis.text = element_text(size = 9),
           plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  return(p)
}