#' Generate Pseudo-bulk Expression Matrix from Seurat Object
#'
#' Aggregates single-cell RNA-seq counts into pseudo-bulk expression profiles by sample.
#' Optionally filters for a specific cell type before aggregation.
#'
#' @param seu A Seurat object containing single-cell RNA-seq data.
#' @param group_by Character string specifying the metadata column name to group cells by
#'   (e.g., "orig.ident", "sample_id"). This determines the columns of the output matrix.
#' @param cell_col Optional character string specifying the metadata column for cell type
#'   annotation (e.g., "cell_type", "major_celltype"). Required when \code{cell_choose} is provided.
#' @param cell_choose Optional character string specifying a single cell type to subset
#'   (e.g., "Macrophages", "T_cells"). When provided, \code{cell_col} must also be provided
#'   (cannot be NULL) and \code{cell_choose} must be a single value (length = 1).
#'
#' @returns A data frame with genes as rows and samples as columns, containing summed
#'   pseudo-bulk expression counts.
#'
#' @importFrom Seurat GetAssayData
#' @importFrom dplyr filter pull sym
#' @importFrom rlang !!
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' # Aggregate all cells by sample
#' pseudo_bulk <- getpseudobulk(seu, group_by = "orig.ident")
#'
#' # Aggregate only Macrophages by sample (requires cell_col)
#' pseudo_bulk_mac <- getpseudobulk(
#'   seu,
#'   group_by = "orig.ident",
#'   cell_col = "major_celltype",
#'   cell_choose = "Macrophages"
#' )
getpseudobulk <- function(seu, group_by, cell_col = NULL, cell_choose = NULL){

  # Validate parameter constraints
  if (!is.null(cell_choose)) {
    if (is.null(cell_col)) {
      stop("cell_col must be provided when cell_choose is specified")
    }
    if (length(cell_choose) != 1) {
      stop("cell_choose must be a single value")
    }
  }

  counts <- Seurat::GetAssayData(seu, assay = "RNA", slot = "counts")
  metadata <- seu@meta.data

  if(!is.null(cell_col) && !is.null(cell_choose)){
    barcodes_choose <- metadata %>% dplyr::filter(!!sym(cell_col) == cell_choose) %>% rownames()
    counts <- counts[, barcodes_choose]
    metadata <- metadata[barcodes_choose, ]
  }

  samples <- metadata %>% dplyr::pull(!!sym(group_by)) %>% as.factor()
  names(samples) <- rownames(metadata)

  mat.summary <- do.call(cbind, lapply(levels(samples), function(x){
    cells <- names(samples)[samples == x]
    pseudobulk <- rowSums(counts[, cells])
    return(pseudobulk)
  })) %>% as.data.frame()
  colnames(mat.summary) <- levels(samples)

  return(mat.summary)
}
