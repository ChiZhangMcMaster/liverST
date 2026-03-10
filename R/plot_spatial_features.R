#' Plot Spatial Feature Expression
#'
#' @param obj A Seurat object
#' @param features Character vector of gene names
#' @param file Optional file name for SVG export
#' @param width SVG width (default = 5)
#' @param height SVG height (default = 6)
#'
#' @return A ggplot object
#' @export
plot_spatial_features <- function(
    obj,
    features,
    file = NULL,
    width = 5,
    height = 6
) {

  p <- Seurat::SpatialFeaturePlot(
    object = obj,
    features = features
  )

  if (!is.null(file)) {
    grDevices::graph2svg(
      x = NULL,
      file = file,
      font = "Arial",
      cairo = TRUE,
      width = width,
      height = height,
      bg = "white"
    )
  }

  return(p)
}
