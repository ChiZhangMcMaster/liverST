#' Rank, Filter, and Export Cluster Biomarkers
#'
#' @param markers Data frame returned from FindAllMarkers
#' @param logfc.threshold Minimum avg_log2FC threshold (default = 0.5)
#' @param padj.threshold Adjusted p-value cutoff (default = 0.05)
#' @param min.pct.1 Minimum pct.1 expression (default = 0)
#' @param max.pct.2 Maximum pct.2 expression (default = 1)
#' @param min.pct.diff Minimum difference between pct.1 and pct.2 (optional)
#' @param top_n Number of top markers per cluster (default = 10)
#' @param file File name for CSV output (default = "cluster_markers.csv")
#'
#' @return Filtered and ranked marker data frame (also exported to CSV)
#' @export
rank_cluster_biomarkers <- function(
    markers,
    logfc.threshold = 0.5,
    padj.threshold = 0.05,
    min.pct.1 = 0,
    max.pct.2 = 1,
    min.pct.diff = NULL,
    top_n = 10,
    file = "cluster_markers.csv"
) {

  if (!is.data.frame(markers)) {
    stop("markers must be a data frame returned from FindAllMarkers.")
  }

  # Adjusted p-value filter
  markers <- dplyr::filter(markers, p_val_adj < padj.threshold)

  # Log2FC filter
  markers <- dplyr::filter(markers, avg_log2FC > logfc.threshold)

  # Expression percentage filters
  markers <- dplyr::filter(markers, pct.1 >= min.pct.1)
  markers <- dplyr::filter(markers, pct.2 <= max.pct.2)

  # Optional pct difference filter
  if (!is.null(min.pct.diff)) {
    markers <- dplyr::filter(markers, (pct.1 - pct.2) >= min.pct.diff)
  }

  # Rank within cluster
  markers <- dplyr::group_by(markers, cluster)
  markers <- dplyr::slice_max(markers, order_by = avg_log2FC, n = top_n)
  top_markers <- dplyr::ungroup(markers)

  # Export CSV
  utils::write.csv(top_markers, file = file, row.names = FALSE)
  message("Marker table exported to: ", file)

  return(top_markers)
}
