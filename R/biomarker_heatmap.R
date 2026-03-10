#' Generate Heatmap from Marker Table
#'
#' @param obj A Seurat object
#' @param marker_table A data frame containing marker genes
#'
#' @return A ggplot heatmap object
#' @export
biomarker_heatmap <- function(
    obj,
    marker_table
) {

  Seurat::DoHeatmap(
    object = obj,
    features = marker_table$gene
  ) +
    Seurat::NoLegend()
}
