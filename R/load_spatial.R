#' Load and standardize 10X Visium spatial transcriptomics data
#'
#' Wrapper around Seurat::Load10X_Spatial with automatic
#' metadata standardization for multi-sample workflows.
#'
#' @param data_dir Path to the Visium data directory
#' @param sample_id Character. Unique sample identifier
#' @param filename HDF5 file name
#' @param assay Assay name (default: "Spatial")
#' @param filter.matrix Logical; filter matrix
#' @param to.upper Logical; convert gene names to upper case
#'
#' @return A standardized Seurat object
#' @export
load_spatial <- function(
    data_dir,
    sample_id,
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    filter.matrix = TRUE,
    to.upper = FALSE
) {

  if (!is.character(sample_id) || length(sample_id) != 1) {
    stop("sample_id must be a single character string.")
  }

  # 1. Load data (slice = sample_id)
  obj <- Seurat::Load10X_Spatial(
    data.dir = data_dir,
    filename = filename,
    assay = assay,
    slice = sample_id,
    filter.matrix = filter.matrix,
    to.upper = to.upper
  )

  # 2. Add metadata
  obj$sample_id <- sample_id
  obj$orig.ident <- sample_id

  # 3. Set identities
  Seurat::Idents(obj) <- obj$sample_id

  # 4. Prefix cell names (critical for integration safety)
  obj <- Seurat::RenameCells(
    obj,
    add.cell.id = sample_id
  )

  # 5. Set project name
  obj@project.name <- sample_id

  return(obj)
}
