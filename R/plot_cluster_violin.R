#' Plot Violin Plot for Selected Genes
#'
#' @param obj A Seurat object
#' @param features Character vector of gene names
#' @param ncol Number of columns (default = 2)
#' @param file Optional file name for SVG export (without extension)
#' @param width SVG width (default = 6.5)
#' @param height SVG height (default = 5)
#'
#' @return A ggplot object
#' @export
plot_cluster_violin <- function(
    obj,
    features,
    ncol = 2,
    file = NULL,
    width = 6.5,
    height = 5
) {

  p <- Seurat::VlnPlot(
    object = obj,
    features = features,
    ncol = ncol
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
