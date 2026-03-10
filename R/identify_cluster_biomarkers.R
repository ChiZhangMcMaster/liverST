#' Identify Cluster Biomarkers
#'
#' @param obj A Seurat object
#' @param only.pos Logical; return only positive markers (default = TRUE)
#' @param assay Assay to use (default = DefaultAssay(obj))
#'
#' @return Data frame of all marker results
#' @export
identify_cluster_biomarkers <- function(obj, only.pos = TRUE) {

  markers <- Seurat::FindAllMarkers(
    object = obj,
    only.pos = only.pos
  )

  return(markers)
}
