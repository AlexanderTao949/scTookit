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

#' Visualize average expression and percentage of expressed genes by cluster
#'
#' This function generates a dot plot visualizing marker gene expression patterns
#' across cell clusters. Dot size represents the percentage of cells expressing a gene,
#' while color represents relative expression (scaled per gene across clusters).
#'
#' @param markerlist A named list of gene vectors. Each element name represents a functional
#'   group (e.g., pathways), and the vector contains associated marker genes.
#' @param seu A Seurat object containing single-cell data.
#' @param assay Name of assay to use (default: "RNA").
#' @param slot Data slot to use (default: "counts").
#' @param scale Logical value indicating whether to apply Z-score normalization to gene expression values.
#'   If \code{TRUE} (default), expression values for each gene are scaled across all cell groups
#'   (mean = 0, SD = 1), with colors representing deviation from the mean (red = high, blue = low).
#'   This facilitates cross-gene comparison of relative expression patterns.
#'   If \code{FALSE}, raw expression values are log1p-transformed without standardization,
#'   preserving absolute expression magnitude but making cross-gene comparison less meaningful.
#'   Similar to the \code{scale} parameter in Seurat's DotPlot function.
#' @param zscore_cutoff Numeric value specifying the truncation threshold for Z-scores.
#'   Default is \code{2.5}. Extreme Z-score values (e.g., > 2.5 or < -2.5) are winsorized
#'   (capped at this threshold) to prevent outlier cell groups from compressing the color scale
#'   resolution for the majority of data points. This ensures that subtle expression differences
#'   across most groups remain visually distinguishable.
#'   Only applicable when \code{scale = TRUE}.
#' @param facet A logical value that controls whether faceted display is enabled.
#'   - If `TRUE` (default): enables faceting using the \code{names(markerlist)}.
#'   - If `FALSE`: generates a unified plot without faceting.
#'
#' @return A ggplot object showing marker expression patterns across clusters.
#'
#' @importFrom Seurat Idents WhichCells GetAssayData
#' @importFrom Matrix rowMeans rowSums
#' @importFrom ggplot2 ggplot aes geom_point scale_color_gradientn guide_colorbar
#'   scale_size_continuous facet_grid theme_minimal theme element_text element_rect
#'   element_blank xlab ylab
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats setNames
#'
#' @examples
#' # Example using Seurat object 'pbmc' and marker list
#' Idents(pbmc) <- pbmc$seurat_clusters
#' Idents(pbmc) <- factor(Idents(pbmc), levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"))
#' marker_list <- list(
#'   T_cell = c("CD3D", "CD8A"),
#'   B_cell = c("CD79A", "MS4A1")
#' )
#' MarkDotplot(markerlist = marker_list, seu = pbmc)
#'
#' @export

MarkDotplot <- function(markerlist = NULL, seu = NULL, assay = "RNA", 
                        slot = "data", facet = TRUE, 
                        scale = FALSE,
                        zscore_cutoff = 2.5) {
  
  if (is.null(seu) || is.null(markerlist)) {
    stop("Both seu and markerlist must be provided")
  }
  
  genes_all <- unlist(markerlist)
  genes_valid <- intersect(genes_all, rownames(seu))
  clusters <- levels(Idents(seu))
  
  if (length(genes_valid) == 0) {
    stop("No valid genes found in the Seurat object")
  }
  
  cluster_cells <- lapply(setNames(nm = clusters), function(cl) {
    WhichCells(seu, ident = cl)
  })
  
  counts_data <- GetAssayData(seu, assay = assay, slot = slot)
  
  result_data <- do.call(rbind, lapply(names(markerlist), function(func_name) {
    genes_func <- intersect(markerlist[[func_name]], genes_valid)
    if (length(genes_func) == 0) return(NULL)
    
    do.call(rbind, lapply(clusters, function(cl) {
      cells <- cluster_cells[[cl]]
      mat <- counts_data[genes_func, cells, drop = FALSE]
      avg_expr <- Matrix::rowMeans(mat)
      pct_expr <- Matrix::rowSums(mat > 0) / length(cells) * 100
      
      data.frame(
        Function = func_name,
        Gene = genes_func,
        Group = cl,
        AvgExpr = avg_expr,
        PctExpr = pct_expr,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  if (scale) {
    result_data <- result_data[order(result_data$Gene, result_data$Group), ]

    avg.exp.scaled <- sapply(unique(result_data$Gene), function(x) {
      data.use <- result_data[result_data$Gene == x, "AvgExpr"]
      if (length(data.use) > 1 && sd(data.use) > 0) {
        data.use <- scale(data.use)[,1]  
        data.use[data.use < -zscore_cutoff] <- -zscore_cutoff
        data.use[data.use > zscore_cutoff] <- zscore_cutoff
      } else {
        data.use <- rep(0, length(data.use))  
      }
      return(data.use)
    })
    
    result_data$AvgExpr <- as.vector(t(avg.exp.scaled))
    color_title <- "Average\nexpression"
    color_limits <- c(-zscore_cutoff, zscore_cutoff)
  } else {
    result_split <- split(result_data$AvgExpr, result_data$Gene)
    max_vals <- vapply(result_split, max, numeric(1))
    max_vals[max_vals <= 0] <- 1
    result_data$AvgExpr <- result_data$AvgExpr / max_vals[result_data$Gene]
    color_title <- "Relative\nExpression"
    color_limits <- NULL
  }

  result_data$Function <- factor(result_data$Function, levels = names(markerlist))
  unique_genes <- rev(genes_all[!duplicated(genes_all)])
  result_data$Gene <- factor(result_data$Gene, levels = unique_genes)
  result_data$Group <- factor(result_data$Group, levels = rev(levels(Idents(seu))))
  
  p <- ggplot(result_data,
              aes(y = Group, x = Gene,
                  color = AvgExpr, size = PctExpr)) +
    geom_point() +
    scale_color_gradientn(
      colours = rev(brewer.pal(n = 10, name = "RdBu")),
      limits = color_limits,
      guide = guide_colorbar(
        ticks.colour = "black", 
        frame.colour = "black",
        title = color_title
      ),
      name = color_title
    ) +
    scale_size_continuous(
      name = "Percentage\nExpressed",
      range = c(0, 6)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.border = element_rect(fill = NA, color = "black", size = 1.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black", size = 1.5),
      strip.text = element_text(face = "bold"),
      panel.background = element_rect(fill = "white", color = "black", size = 1.5),
      plot.background = element_rect(fill = "white")
    ) +
    xlab("Genes") +
    ylab("Cell Clusters")
  
  if (facet) {
    p <- p + facet_grid(
      . ~ Function,
      scales = "free_x",
      space = "free_x"
    )
  }
  
  return(p)
}
