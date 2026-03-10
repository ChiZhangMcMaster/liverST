#' Percentile-based feature filtering for Seurat objects
#'
#' @param obj A Seurat object
#' @param feature_col Metadata column to filter on
#' @param lower_pct Lower percentile (default = 0.01)
#' @param upper_pct Upper percentile (default = 0.99)
#'
#' @return A filtered Seurat object
#' @export
quantile_filter <- function(
    obj,
    feature_col = "nFeature_Spatial",
    lower_pct = 0.01,
    upper_pct = 0.99
) {

  if (!inherits(obj, "Seurat")) {
    stop("obj must be a Seurat object.")
  }

  if (!feature_col %in% colnames(obj@meta.data)) {
    stop(paste0("Column '", feature_col, "' not found in metadata."))
  }

  values <- obj@meta.data[[feature_col]]

  q_low <- stats::quantile(values, lower_pct, na.rm = TRUE)
  q_high <- stats::quantile(values, upper_pct, na.rm = TRUE)

  keep_cells <- colnames(obj)[
    values > q_low & values < q_high
  ]

  obj_filtered <- subset(obj, cells = keep_cells)

  kept <- length(keep_cells)
  total <- ncol(obj)
  pct <- round(100 * kept / total, 2)

  message(
    paste0(
      "Filtering based on ", feature_col,
      " (", lower_pct * 100, "–", upper_pct * 100, " percentile): ",
      kept, "/", total, " retained (", pct, "%)"
    )
  )

  return(obj_filtered)
}
