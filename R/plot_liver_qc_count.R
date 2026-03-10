#' Plot spatial count QC for liver Seurat objects
#'
#' @param seurat_obj A Seurat object
#' @param file Character or NULL. Output filename (without extension) for SVG export
#' @param width Numeric. Width of saved figure
#' @param height Numeric. Height of saved figure
#'
#' @return A patchwork plot object
#' @export
plot_liver_qc_count <- function(
    seurat_obj,
    file = NULL,
    width = 6,
    height = 4.5
) {

  feature <- "nCount_Spatial"

  p1 <- Seurat::VlnPlot(
    seurat_obj,
    features = feature,
    pt.size = 0.1
  ) +
    Seurat::NoLegend()

  p2 <- Seurat::SpatialFeaturePlot(
    seurat_obj,
    features = feature,
    slot = "counts"
  ) +
    ggplot2::theme(legend.position = "right")

  p <- patchwork::wrap_plots(p1, p2)

  if (!is.null(file)) {
    grDevices::svg(
      filename = paste0(file, ".svg"),
      width = width,
      height = height,
      bg = "transparent"
    )
    print(p)
    grDevices::dev.off()
  }

  return(p)
}
