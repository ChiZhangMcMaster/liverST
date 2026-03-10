#' Integrate Multiple SCT-Normalized Seurat Objects Using Harmony
#'
#' @description
#' Performs merging and Harmony-based integration across multiple SCT-normalized
#' Seurat objects, followed by clustering and UMAP visualization.
#'
#' @param object.list A list of SCT-normalized Seurat objects.
#' @param project.name Project name for merged object.
#' @param npcs Number of PCs to compute (default = 30).
#' @param dims Dimensions to use for downstream analysis (default = 1:30).
#' @param resolution Clustering resolution (default = 0.8).
#' @param verbose Logical; whether to print progress messages.
#'
#' @return Integrated Seurat object with Harmony reduction and UMAP embedding.
#' @export
integrate_ST_samples <- function(object.list,
                                 project.name = "STprotocol",
                                 npcs = 30,
                                 dims = 1:30,
                                 resolution = 0.8,
                                 nfeatures = 3000,
                                 verbose = TRUE) {

  if (length(object.list) < 2) {
    stop("At least two Seurat objects are required for integration.")
  }

  # Merge
  merged.obj <- merge(
    x = object.list[[1]],
    y = object.list[-1],
    project = project.name
  )

  Seurat::DefaultAssay(merged.obj) <- "SCT"

  # Select integration features
  var.features <- Seurat::SelectIntegrationFeatures(
    object.list = object.list,
    nfeatures = nfeatures
  )
  Seurat::VariableFeatures(merged.obj) <- var.features

  # PCA
  merged.obj <- Seurat::RunPCA(merged.obj, npcs = npcs, verbose = FALSE)

  # Harmony integration
  merged.obj <- Seurat::IntegrateLayers(
    object = merged.obj,
    method = Seurat::HarmonyIntegration,
    orig.reduction = "pca",
    normalization.method = "SCT",
    new.reduction = "harmony",
    verbose = verbose
  )

  # Neighbors
  merged.obj <- Seurat::FindNeighbors(
    merged.obj,
    reduction = "harmony",
    dims = dims
  )

  # Clustering
  merged.obj <- Seurat::FindClusters(
    merged.obj,
    resolution = resolution,
    cluster.name = "cluster",
    verbose = FALSE
  )

  # Explicitly set identities to cluster
  Seurat::Idents(merged.obj) <- merged.obj$cluster

  # UMAP
  merged.obj <- Seurat::RunUMAP(
    merged.obj,
    reduction = "harmony",
    dims = dims,
    reduction.name = "umap.harmony"
  )

  # Prepare SCT assay for FindMarkers / FindAllMarkers
  merged.obj <- Seurat::PrepSCTFindMarkers(merged.obj, assay = "SCT", verbose = verbose)

  if (verbose) message("Integration complete: object is ready for DE testing.")

  return(merged.obj)
}
