#' Plot spatial regions from a Seurat object with optional rotation
#'
#' @param seurat_obj A Seurat object with spatial data
#' @param group_by Metadata column to color by (default = "region")
#' @param image Name of spatial image
#' @param alpha Transparency vector (default = c(0.5, 0.5))
#' @param pt_size_factor Point size scaling factor (default = 1.6)
#' @param rotate_degree Degrees to rotate plot (0, 90, 180, 270)
#'
#' @return A ggplot object
#' @export
plot_spatial_regions <- function(
    seurat_obj,
    group_by = "region",
    image,
    alpha = c(0.5, 0.5),
    pt_size_factor = 1.6,
    rotate_degree = 0
) {

  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.")
  }

  if (!group_by %in% colnames(seurat_obj@meta.data)) {
    stop("group_by column not found in Seurat metadata.")
  }

  # Base plot
  p <- Seurat::SpatialDimPlot(
    object = seurat_obj,
    group.by = group_by,
    images = image,
    alpha = alpha,
    pt.size.factor = pt_size_factor
  )

  return(p)
}

