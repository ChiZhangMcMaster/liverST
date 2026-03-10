#' Run PCA, Clustering, and UMAP on a Seurat Object
#'
#' Performs PCA, neighbor finding, clustering, and UMAP dimensional
#' reduction on a Seurat object.
#'
#' @param object A Seurat object.
#' @param assay Assay to use for PCA (default = "SCT").
#' @param dims Dimensions to use for PCA and downstream steps (default = 1:30).
#' @param resolution Clustering resolution (default = 1).
#'
#' @return Processed Seurat object
#' @export
run_umap_clustering <- function(
    object,
    assay = "SCT",
    dims = 1:30,
    resolution = 1
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

  object <- Seurat::RunUMAP(
    object,
    reduction = "pca",
    dims = dims
  )

  return(object)
}

