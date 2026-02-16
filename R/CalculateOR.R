#' Cell-wise Fisher's Exact Test for Contingency Tables
#'
#' Performs Fisher's exact test on each cell of a contingency table to evaluate
#' associations between row and column categories. For every cell (i,j), constructs
#' a 2x2 contingency table comparing the observed count against the row and column
#' marginals, then calculates odds ratios and p-values.
#'
#' @param count.dist A numeric matrix or data.frame representing the contingency
#'   table, where rows correspond to one categorical variable and columns to another.
#'   Row and column names are preserved in the output.
#' @param min.rowSum Numeric threshold for row filtering. Rows with marginal sums
#'   strictly less than this value are excluded from analysis. Default is 0.
#'
#' @return A data.table in long format containing:
#'   \itemize{
#'     \item rid: Row identifier (from input rownames)
#'     \item cid: Column identifier (from input colnames)
#'     \item count: Observed frequency in the cell
#'     \item p.value: Two-sided p-value from Fisher's exact test
#'     \item OR: Estimated odds ratio for the cell against all others
#'   }
#'
#' @importFrom data.table setDT melt as.data.table
#' @importFrom plyr ldply
#' @importFrom stats fisher.test
test.dist.table <- function(count.dist, min.rowSum=0){

  count.dist <- count.dist[rowSums(count.dist) >= min.rowSum, ,drop=F]
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(count.dist)

  setDT(count.dist.tb, keep.rownames=T)
  count.dist.melt.tb <- melt(count.dist.tb, id.vars="rn")
  colnames(count.dist.melt.tb) <- c("rid", "cid", "count")

  count.dist.melt.ext.tb <- as.data.table(ldply(seq_len(nrow(count.dist.melt.tb)), function(i){

    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c   <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col] - this.c

    this.m <- matrix(c(this.c,
                       sum.row[this.row] - this.c,
                       other.col.c,
                       sum(sum.col) - sum.row[this.row] - other.col.c),
                     ncol=2)
    res.test <- fisher.test(this.m)

    data.frame(rid = this.row,
               cid = this.col,
               p.value = res.test$p.value,
               OR  = res.test$estimate)
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb, count.dist.melt.ext.tb,
                                  by=c("rid", "cid"))
  return(count.dist.melt.ext.tb)
}

#' Calculate Directional Odds Ratios for Cell Type Enrichment Between Two Groups
#'
#' Computes cell-wise odds ratios (OR) to quantify enrichment of cell types
#' across two specified groups/conditions within a Seurat object. For each cell
#' type, the function determines which group shows stronger association, retains
#' the higher odds ratio, and assigns a directional sign (positive for the
#' reference group, negative for the alternative group).
#'
#' @param seu A \code{\link[Seurat]{Seurat}} object containing single-cell
#'   RNA-seq data and associated metadata.
#' @param metacol Character string specifying the metadata column name
#'   containing the grouping variable (e.g., "condition", "treatment", "disease_state").
#' @param metacol_values Character vector of length 2 specifying the two groups
#'   to compare. The first element serves as the alternative (positive direction),
#'   and the second as the reference (negative direction) (e.g., \code{metacol_values = c("treat", "control")}).
#'   These values must exist as levels in \code{metacol}.
#' @param cellcol Character string specifying the metadata column name
#'   containing cell type identities or cluster annotations (e.g., "cell_type",
#'   "seurat_clusters", "annotation").
#'
#' @returns A data.table with columns:
#'   \itemize{
#'     \item \strong{celltype}: Factor (ordered by OR), the cell type identifier
#'     \item \strong{group}: Character, the group showing stronger enrichment
#'     \item \strong{OR}: Numeric, signed odds ratio where positive values
#'       indicate enrichment in \code{metacol_values[1]} and negative values
#'       indicate enrichment in \code{metacol_values[2]}
#'   }
#'   Only cell types with finite odds ratios for both groups are retained.
#'   The data is sorted by OR in descending order.
#'
#'
#' This signed odds ratio format is particularly useful for generating
#' diverging bar plots or lollipop charts showing bidirectional enrichment.
#'
#' @importFrom data.table dcast melt
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Calculate directional enrichment between disease and control
#'   or_results <- CalculateOR(
#'     seu = seurat_obj,
#'     metacol = "disease_status",
#'     metacol_values = c("disease", "control"),
#'     cellcol = "cell_type"
#'   )
#'
#'   # Plot results
#'   library(ggplot2)
#'   ggplot(or_results, aes(x = celltype, y = OR, fill = group)) +
#'     geom_col() +
#'     coord_flip() +
#'     labs(y = "Signed Odds Ratio", x = "Cell Type")
#' }
CalculateOR <- function(seu, metacol, metacol_values, cellcol){

  cellInfo.tb  = as.data.frame(seu@meta.data)
  # metacol_values = unique(cellInfo.tb[[metacol]])
  meta.cluster = cellInfo.tb[[cellcol]]
  cellInfo.tb$meta.cluster = as.character(meta.cluster)

  loc = cellInfo.tb[[metacol]]
  if(is.factor(loc)){
    cellInfo.tb$loc = loc
  }else{cellInfo.tb$loc = as.factor(loc)}
  loc.avai.vec <- levels(cellInfo.tb[["loc"]])

  count.dist <- table(cellInfo.tb$meta.cluster, cellInfo.tb$loc)
  count.dist <- unclass(count.dist)[, loc.avai.vec, drop = FALSE]
  freq.dist <- sweep(count.dist, 1, rowSums(count.dist), "/")
  freq.dist.bin <- floor(freq.dist * 100 / 10)
  count.dist.melt.ext.tb <- test.dist.table(count.dist)

  p.dist.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p.value")
  test <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var="OR")

  test <- test[is.finite(test[[metacol_values[1]]]),]
  test <- test[is.finite(test[[metacol_values[2]]]),]
  for (i in 1:length(test$rid)){
    ifelse(test[[metacol_values[1]]][[i]]>test[[metacol_values[2]]][[i]], test[[metacol_values[2]]][[i]]<-0, test[[metacol_values[1]]][[i]]<-0)
  }
  test <- data.table::melt(test, id.vars = "rid")
  test <- test[test$value>0,]
  colnames(test) <- c('celltype', 'group', 'OR')
  for (i in 1:length(test$celltype)){
    if (test$group[i] == metacol_values[2]){
      test$OR[i] = 0 - test$OR[i]
    }
  }
  test <- test[order(test$OR,decreasing = TRUE),]
  test$celltype <- factor(test$celltype, levels = test$celltype)

  return(test)
}


#' Plot Odds Ratio (OR) Values with Bar and Point Geoms
#'
#' Creates a ggplot visualization showing odds ratio values as bars with overlaid
#' points colored by group. The y-axis is symmetric around zero to facilitate
#' comparison of effect sizes in both directions.
#'
#' @param data A data frame containing the variables to be plotted.
#' @param x_var Character string specifying the column name for the x-axis
#'   (categorical variable, e.g., cell types). Default is "celltype".
#' @param y_var Character string specifying the column name for the y-axis
#'   (numeric variable, typically odds ratios). Default is "OR".
#' @param group_var Character string specifying the column name for grouping
#'   (color/fill aesthetic). Default is "group".
#' @param colors Character vector of colors for the groups. Default is
#'   `c("#D8383A", "dodgerblue4", "snow3")`.
#' @param bar_width Numeric value specifying the width of the bars. Default is 0.1.
#' @param point_size Numeric value specifying the size of the points. Default is 6.
#' @param x_text_size Numeric value specifying the font size for x-axis labels.
#'   Default is 6.
#' @param y_text_size Numeric value specifying the font size for y-axis labels.
#'   Default is 8.
#' @param x_angle Numeric value specifying the rotation angle for x-axis labels.
#'   Default is 90.
#' @param hline_color Character string specifying the color of the horizontal
#'   reference line at y=0. Default is "black".
#' @param hline_size Numeric value specifying the line width of the horizontal
#'   reference line. Default is 0.5.
#'
#' @returns A ggplot object that can be further modified or printed.
#'
#' @importFrom rlang sym
#' @importFrom ggplot2 ggplot aes geom_bar geom_point scale_color_manual
#'   scale_fill_manual geom_hline theme_classic ylab scale_y_continuous
#'   theme element_text
#'
#' @export
PlotOR <- function(data,
                   x_var = "celltype",
                   y_var = "OR",
                   group_var = "group",
                   colors = c("#D8383A", "dodgerblue4", "snow3"),
                   bar_width = 0.1,
                   point_size = 6,
                   x_text_size = 8,
                   y_text_size = 10,
                   x_angle = 90,
                   hline_color = "black",
                   hline_size = 0.5) {

  x_sym <- rlang::sym(x_var)
  y_sym <- rlang::sym(y_var)
  group_sym <- rlang::sym(group_var)

  y_max <- max(abs(data[[y_var]]), na.rm = TRUE)

  p <- ggplot(data, aes(x = !!x_sym, y = !!y_sym)) +
    geom_bar(stat = "identity", width = bar_width) +
    geom_point(size = point_size, aes(color = !!group_sym, fill = !!group_sym)) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    geom_hline(yintercept = 0, color = hline_color, linewidth = hline_size) +
    theme_classic() +
    ylab(y_var) +
    scale_y_continuous(limits = c(-y_max, y_max)) +
    theme(
      axis.text.x = element_text(size = x_text_size, angle = x_angle, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = y_text_size)
    )

  return(p)
}
