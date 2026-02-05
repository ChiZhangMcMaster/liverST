#' Load 10X Visium spatial transcriptomics data
#'
#' A wrapper around Seurat::Load10X_Spatial for simplified loading.
#'
#' @param data_dir Path to the Visium data directory
#' @param filename HDF5 file name
#' @param assay Assay name (default: "Spatial")
#' @param slice Slice name
#' @param filter.matrix Logical; filter matrix
#' @param to.upper Logical; convert gene names to upper case
#'
#' @return A Seurat object
#' @export
load_spatial <- function(
    data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice,
    filter.matrix = TRUE,
    to.upper = FALSE
) {
  Seurat::Load10X_Spatial(
    data.dir = data_dir,
    filename = filename,
    assay = assay,
    slice = slice,
    filter.matrix = filter.matrix,
    to.upper = to.upper
  )
}
