#' Create a Scatter Plot with Correlation Statistics
#'
#' Generates a scatter plot with optional regression line, confidence interval, 
#' and correlation coefficient annotation using ggpubr and ggplot2.
#'
#' @param data A data frame containing the variables to be plotted.
#' @param x_var Character string specifying the column name for the x-axis variable.
#' @param y_var Character string specifying the column name for the y-axis variable.
#' @param pointsize Numeric value specifying the size of points. Default is 1.5.
#' @param alpha Numeric value between 0 and 1 specifying point transparency. Default is 0.6.
#' @param add Character string specifying the type of line to add (e.g., "reg.line", "none"). Default is "reg.line".
#' @param line_color Character string specifying the color of the regression line. Default is "#1F78B4".
#' @param fill_color Character string specifying the fill color for confidence interval. Default is "#A6CEE3".
#' @param linesize Numeric value specifying the thickness of the regression line. Default is 1.
#' @param conf Logical value indicating whether to display confidence interval. Default is TRUE.
#' @param cor_method Character string specifying correlation method: "pearson", "spearman", or "kendall". Default is "spearman".
#' @param label.x.d Numeric value for horizontal positioning of correlation label (as proportion of x-range from minimum). Default is 0.1.
#' @param label.y.d Numeric value for vertical positioning of correlation label (as proportion of y-range from maximum). Default is 0.1.
#' @param Xlab Character string for x-axis label.
#' @param Ylab Character string for y-axis label.
#' @param title Character string for plot title.
#' @param textsize Numeric value specifying the text size of correlation annotation. Default is 4.
#'
#' @return A ggplot object representing the scatter plot with correlation statistics.
#' 
#' @importFrom ggpubr ggscatter stat_cor
#' @importFrom ggplot2 xlab ylab ggtitle theme element_text
#' 
#' @export 
#' 
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' df <- data.frame(
#'   height = rnorm(100, 170, 10),
#'   weight = rnorm(100, 70, 15)
#' )
#' 
#' # Basic scatter plot with Spearman correlation
#' ScatterPlot(data = df, x_var = "height", y_var = "weight",
#'             Xlab = "Height (cm)", Ylab = "Weight (kg)",
#'             title = "Height vs Weight")
#'
#' # Customized plot with Pearson correlation
#' ScatterPlot(data = df, x_var = "height", y_var = "weight",
#'             pointsize = 2, alpha = 0.8,
#'             line_color = "darkred", fill_color = "lightcoral",
#'             cor_method = "pearson",
#'             Xlab = "Height (cm)", Ylab = "Weight (kg)",
#'             title = "Correlation Analysis")
#' }
ScatterPlot <- function(data, x_var, y_var, pointsize = 1.5, alpha = 0.6,
                        add = "reg.line", line_color = "#1F78B4", fill_color = "#A6CEE3",
                        linesize = 1, conf = TRUE, cor_method = "spearman", label.x.d = .1, label.y.d = .1,
                        Xlab, Ylab, title, textsize = 4){

  if (!requireNamespace("ggpubr", quietly = TRUE)) stop("Package 'ggpubr' is required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required")
                          
  corpval_x <- min(data[[x_var]]) + label.x.d
  corpval_y <- max(data[[y_var]]) - label.y.d
  
  plot <- ggpubr::ggscatter(data, x = x_var, y = y_var,
                            size = pointsize,
                            alpha = alpha,
                            add = add,
                            add.params = list(color = line_color, fill = fill_color, size = linesize),
                            conf.int = conf) +
    ggpubr::stat_cor(method = cor_method, label.x = corpval_x, label.y = corpval_y, label.sep = "\n", size = textsize) +
    ggplot2::xlab(Xlab) +
    ggplot2::ylab(Ylab) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                   aspect.ratio = 1)
  return(plot)
}