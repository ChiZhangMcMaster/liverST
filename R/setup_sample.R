#' Standardize Seurat metadata for downstream analyses
#'
#' @param seurat_obj A Seurat object
#' @param condition Character. Condition/sample name to assign
#'
#' @return A modified Seurat object
#' @export
setup_sample <- function(
    seurat_obj,
    condition
) {

  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object.")
  }

  if (!is.character(condition) || length(condition) != 1) {
    stop("condition must be a single character string.")
  }

  # 1. Set orig.ident
  seurat_obj$orig.ident <- condition

  # 2. Set identities to condition
  Seurat::Idents(seurat_obj) <- seurat_obj$orig.ident

  # 3. Set project name (optional but clean)
  seurat_obj@project.name <- condition

  return(seurat_obj)
}

