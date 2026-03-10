#' Run PCA, Clustering, and t-SNE on a Seurat Object
#'
#' Performs PCA, neighbor finding, clustering, and t-SNE dimensional
#' reduction on a Seurat object. Designed for spatial transcriptomics
#' workflows but compatible with any Seurat object.
#'
#' @param object A Seurat object.
#' @param assay Assay to use for PCA (default = "SCT").
#' @param dims Dimensions to use for PCA and downstream steps (default = 1:30).
#' @param resolution Clustering resolution (default = 1).
#' @param return_plots Logical; whether to return t-SNE and spatial plots (default = TRUE).
#'
#' @return If return_plots = FALSE, returns a processed Seurat object.
#' If return_plots = TRUE, returns a list containing:
#' \itemize{
#'   \item object: Processed Seurat object
#'   \item tsne_plot: t-SNE DimPlot
#'   \item spatial_plot: SpatialDimPlot
#'   \item combined_plot: Combined t-SNE + spatial plot
#' }
#'
#' @export
run_tsne_clustering <- function(
    object,
    assay = "SCT",
    dims = 1:30,
    resolution = 1,
    return_plots = TRUE
) {

  object <- Seurat::RunPCA(object, assay = assay, verbose = FALSE)

  object <- Seurat::FindNeighbors(
    object,
    reduction = "pca",
    dims = dims
  )

  object <- Seurat::FindClusters(
    object,
    resolution = resolution,
    verbose = FALSE
  )

  object <- Seurat::RunTSNE(
    object,
    reduction = "pca",
    dims = dims
  )

  return(object)
}
