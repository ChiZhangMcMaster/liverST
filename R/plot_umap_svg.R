#' Plot and save UMAP DimPlot
#'
#' @param obj Seurat object
#' @param reduction Dimensionality reduction to use (default "umap.harmony")
#' @param label Logical; whether to label clusters
#' @param group.by Metadata column(s) to group by (default: "ident")
#' @param file File name for SVG output (without extension)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param font Font family for SVG (default "Arial")
#' @param bg Background color (default "transparent")
#' @export
plot_umap_svg <- function(obj,
                          reduction = "umap.harmony",
                          label = TRUE,
                          group.by = "ident",
                          file = "DimPlot",
                          width = 10,
                          height = 7,
                          font = "Arial",
                          bg = "transparent") {

  p <- Seurat::DimPlot(
    obj,
    reduction = reduction,
    label = label,
    group.by = group.by
  )

  # Save as SVG
  export::graph2svg(
    x = p,
    file = paste0(file, ".svg"),
    font = font,
    cairo = TRUE,
    width = width,
    height = height,
    bg = bg
  )

  return(p)
}
