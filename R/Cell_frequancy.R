#' Visualize cluster composition across groups
#'
#' This function generates a two-panel visualization showing the distribution of
#' cell clusters across different experimental groups or samples. The left panel
#' displays the proportional composition of each cluster (stacked percentage bar),
#' while the right panel shows the absolute cell counts per cluster.
#'
#' @param data A Seurat object containing single-cell RNA-seq data
#' @param clusters Character string specifying the column name in meta.data
#'   containing cluster identities (default: "seurat_clusters")
#' @param group Character string specifying the column name in meta.data
#'   defining the grouping variable (e.g., "orig.ident", "sample", "condition")
#' @param widths Numeric vector of length 2 specifying relative widths for the
#'   composition plot and cell count plot (default: c(3, 1))
#' @param log Logical indicating whether to display cell counts on log10 scale
#'   in the right panel (default: TRUE)
#' @param legend.title Character string for the legend title (default: "Group")
#' @param color Integer specifying color palette: 1 = D3 Category 20c
#'   (supports up to 20 groups), 2 = ColorBrewer Set2 (supports up to 8 groups)
#' @param xlab Character string for x-axis label of the composition plot
#'   (default: "")
#'
#' @returns A ggplot object combining two plots via patchwork: left panel shows
#'   stacked percentage bars of group composition per cluster, right panel shows
#'   total cell counts per cluster. Returns NULL if input data is invalid.
#'
#' @details
#' Clusters are automatically sorted: numeric cluster IDs (e.g., "0", "1") are
#' sorted in descending order; non-numeric IDs are sorted by cell count
#' (ascending). The function handles up to 20 groups with color scheme 1 or
#' up to 8 groups with color scheme 2, issuing warnings if group numbers exceed
#' these limits.
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme_bw coord_flip scale_fill_manual
#'   scale_fill_brewer scale_x_log10 scale_y_continuous guides guide_legend
#'   theme xlab ylab element_text
#' @importFrom patchwork wrap_plots
#' @importFrom data.table melt
#' @importFrom paletteer paletteer_d
#' @importFrom scales percent
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' ClusterCompPlot(data = seurat_obj,
#'                 clusters = "seurat_clusters",
#'                 group = "orig.ident")
#'
#' # Custom grouping with condition labels and Set2 colors
#' ClusterCompPlot(data = seurat_obj,
#'                 clusters = "cell_type",
#'                 group = "treatment",
#'                 legend.title = "Treatment",
#'                 color = 2,
#'                 xlab = "Cell Clusters")
#' }
CellCompPlot <-  function(data,
                          clusters,
                          group,
                          widths = c(3,1),
                          log =TRUE,
                          legend.title = "Group",
                          color = 1,
                          xlab = ""){

  mytheme = theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                  axis.title = element_text(size = 10,color ="black"),
                  axis.text = element_text(size=10,color = "black"),
                  #axis.line = element_line(color = "black"),
                  #axis.ticks = element_line(color = "black"),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  # panel.grid=element_blank(),
                  # legend.position = "none",
                  legend.text = element_text(size=8),
                  legend.title= element_text(size= 8),
                  # axis.text.x = element_text(angle = 45, hjust=1, vjust=1)
  )

  count_table <- table(data@meta.data[,clusters], data@meta.data[,group])
  count_mtx <- as.data.frame.matrix(count_table)
  count_mtx$cluster <- rownames(count_mtx)
  melt_mtx <- data.table::melt(count_mtx)
  melt_mtx$cluster <- as.factor(melt_mtx$cluster)

  cluster_size <- aggregate(value ~ cluster, data = melt_mtx, FUN = sum)

  if("0" %in% cluster_size$cluster){
    sorted_labels <- paste(sort(as.integer(levels(cluster_size$cluster)),decreasing = T))
  }else{
    sorted_labels <- paste(cluster_size$cluster[order(cluster_size$value)])
  }

  cluster_size$cluster <- factor(cluster_size$cluster,levels = sorted_labels)
  melt_mtx$cluster <- factor(melt_mtx$cluster,levels = sorted_labels)

  colnames(melt_mtx)[2] <- "dataset"

  ################### p1
  if(log){
    p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") +
      theme_bw() + scale_x_log10() + xlab("Cells per cluster, log10 scale") + ylab("") + mytheme
  }else{
    p1 <- ggplot(cluster_size, aes(y= cluster,x = value)) + geom_bar(position="dodge", stat="identity",fill = "grey60") +
      theme_bw() + xlab("Cells per cluster") + ylab("") + mytheme
  }
  ################### p2
  ########################### color 1
  if(color==1){
    if(length(unique(melt_mtx$dataset)) < 21){
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
        geom_bar(position="fill", stat="identity",) + theme_bw() + coord_flip() +
        scale_fill_manual(values = alpha(paletteer::paletteer_d('ggsci::category20c_d3'), 0.9))+
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
    }else{
      warning("The color limit is <21")
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
        geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
    }
  }
  ########################### color 2
  if(color==2){
    if(length(unique(melt_mtx$dataset)) < 9){
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
        geom_bar(position="fill", stat="identity",) + theme_bw() + coord_flip() +
        scale_fill_brewer(palette = "Set2")+
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
    }else{
      warning("The color limit is <9")
      p2 <- ggplot(melt_mtx,aes(x=cluster,y=value,fill=dataset)) +
        geom_bar(position="fill", stat="identity") + theme_bw() + coord_flip() +
        ylab(paste0("Fraction of cells in each ",tolower(legend.title))) + xlab(xlab) +
        theme(legend.position="top") + guides(fill = guide_legend(title = legend.title)) +
        scale_y_continuous(labels = scales::percent,expand = c(0.01, 0.01))+ mytheme
    }
  }
  wrap_plots(ncol = 2,p2,p1,widths = widths)
}
