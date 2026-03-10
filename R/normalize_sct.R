#' Apply SCTransform normalization to liver Seurat objects
#'
#' This replaces NormalizeData(), ScaleData(), and FindVariableFeatures().
#' Transformed data will be available in the SCT assay, which is set as the default after running SCTransform.
#'
#' @param seurat_objs A Seurat object or a named list of Seurat objects
#' @param assay Character. Assay to normalize (default "Spatial")
#' @param method Character. Method for sctransform (default "glmGamPoi")
#' @param verbose Logical. Whether to print progress messages
#'
#' @return A Seurat object if input is a single object, or a list of Seurat objects if input is a list
#' @export
normalize_sct <- function(
    seurat_objs,
    assay = "Spatial",
    method = "glmGamPoi",
    verbose = TRUE
) {

  # Function to apply SCTransform to a single object
  .apply_sct <- function(obj) {
    Seurat::SCTransform(
      obj,
      method = method,
      assay = assay,
      verbose = verbose
    )
  }

  # Check if input is a list or a single Seurat object
  if (inherits(seurat_objs, "Seurat")) {
    return(.apply_sct(seurat_objs))
  } else if (is.list(seurat_objs)) {
    lapply(seurat_objs, .apply_sct)
  } else {
    stop("Input must be a Seurat object or a list of Seurat objects.")
  }
}
