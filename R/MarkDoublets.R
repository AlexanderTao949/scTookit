#' Calculate pN Parameter for Doublet Detection
#'
#' Computes the pN parameter used in DoubletDecon based on input cell count. The function
#' uses predefined thresholds to determine pN for typical cell counts and extrapolates
#' for counts beyond the maximum threshold. Output is scaled to a proportion (0-1).
#'
#' @param cell_count Numeric. Total number of cells in the dataset. Must be > 0.
#'
#' @return Numeric pN value between 0 and 1 for doublet detection algorithms.
#' @importFrom utils tail
calculate_pN <- function(cell_count) {
  thresholds <- c(500, 1000, 2000, 3000, 4000, 5000,
                  6000, 7000, 8000, 9000)

  pN_values <- c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9,
                 4.6, 5.4, 6.1, 6.9)

  if (cell_count <= 0) {
    stop("Error: Cell count must be greater than 0.")
  }

  idx <- findInterval(cell_count, thresholds) + 1

  if (cell_count > max(thresholds)) {

    base_rate <- tail(pN_values, 1)
    last_threshold <- tail(thresholds, 1)

    extra_cells <- cell_count - last_threshold
    extra_units <- extra_cells / 1000

    last_increments <- diff(tail(pN_values, 4))
    avg_increment <- mean(last_increments)

    pN_percent <- base_rate + (extra_units * avg_increment)
  } else {
    pN_percent <- pN_values[idx]
  }
  pN <- pN_percent / 100

  return(pN)
}

#' Detect doublet cells using Scrublet
#'
#' @description Performs doublet detection by wrapping the Python Scrublet package via reticulate.
#'   Note: Users must manually initialize the Python environment containing the scrublet library
#'   before running this function, e.g., using \code{reticulate::use_condaenv()}.
#'
#' @param seu A Seurat object containing RNA count data
#' @param PN Expected doublet rate. If NULL, automatically calculated using calculate_pN()
#'
#' @return A seurat object whose metadata contain doublet scores and predictions for each barcode.
#'
#' @export
#'
#' @importFrom reticulate import py_run_string
#' @importFrom tibble rownames_to_column
#' @importFrom Matrix t
#' @importFrom dplyr %>%
MarkDoublets <- function(seu, PN = NULL){

  counts_matrix <- Matrix::t(seu@assays$RNA@counts)
  metadata <- seu@meta.data %>% rownames_to_column("barcode")
  gene_names <- colnames(counts_matrix)

  if(is.null(PN)){
    PN <- calculate_pN(nrow(counts_matrix))
  }

  scipy <- reticulate::import("scipy.sparse", convert = FALSE)
  prepare_counts_matrix <- function(seu_counts_t) {
    if (inherits(seu_counts_t, "dgCMatrix")) {
      return(scipy$csc_matrix(seu_counts_t))
    } else {
      return(as.matrix(seu_counts_t))
    }
  }
  counts_matrix <- prepare_counts_matrix(counts_matrix)

  reticulate::py_run_string("
import pandas as pd
import anndata
import scrublet as scr

def detect_doublets(counts_matrix, metadata, gene_names, PN):
    # Create AnnData object
    adata = anndata.AnnData(X=counts_matrix)
    adata.obs = metadata
    if 'barcode' not in adata.obs:
        raise ValueError(\"Metadata must contain 'barcode' column\")
    adata.obs.index = adata.obs['barcode'].astype(str)

    # Check for unique barcodes
    if adata.obs.index.has_duplicates:
        raise ValueError(\"Barcodes must be unique for index assignment\")

    # Set gene names
    adata.var.index = gene_names

    # Run Scrublet
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=PN)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()

    # Create results DataFrame
    df = pd.DataFrame({
        'doublet_score': scrub.doublet_scores_obs_,
        'predicted_doublet': scrub.predicted_doublets_
    }, index=adata.obs.index)

    return df
")

  py <- reticulate::import_main()
  result_df <- py$detect_doublets(
    counts_matrix = counts_matrix,
    metadata = metadata,
    gene_names = gene_names,
    PN = PN
  )

  seu$doublet_score <- result_df$doublet_score
  seu$predicted_doublet <- result_df$predicted_doublet
  return(seu)
}
