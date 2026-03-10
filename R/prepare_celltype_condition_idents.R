#' Prepare Seurat Object for Condition-Specific DE with Control/Treatment Grouping
#'
#' Creates new identity and grouping columns to enable DE testing
#' between conditions within specific cell types, with combined Control/Treatment groups.
#'
#' @param obj A Seurat object.
#' @param celltype_col Optional metadata column containing cell type labels.
#'   If NULL, current Idents(obj) are used.
#' @param condition_col Metadata column containing condition/sample info (default: "orig.ident").
#' @param control_conditions Character vector of condition names to treat as Control.
#' @param new_ident_name Name of the combined celltype + condition column (default: "celltype_condition").
#' @param sep Separator used when combining cell type and condition/group (default: "_").
#'
#' @return Seurat object with new metadata columns:
#'   'celltype_original', 'celltype_condition', 'group', 'celltype_group'.
#'   Also sets Idents(obj) = celltype_condition.
#' @export
prepare_celltype_condition_idents <- function(
    obj,
    celltype_col = NULL,
    condition_col = "orig.ident",
    control_conditions = NULL,
    new_ident_name = "celltype_group",
    sep = "_"
) {

  if (!inherits(obj, "Seurat")) stop("obj must be a Seurat object.")

  md <- obj@meta.data

  if (!condition_col %in% colnames(md)) {
    stop(paste0("Column '", condition_col, "' not found in metadata."))
  }

  # Determine cell type
  celltype_vector <- Seurat::Idents(obj)
  if (!is.null(celltype_col)) {
    if (!celltype_col %in% colnames(md)) {
      stop(paste0("Column '", celltype_col, "' not found in metadata."))
    }
    celltype_vector <- md[[celltype_col]]
  }

  # Original condition vector
  condition_vector <- md[[condition_col]]

  # Create 'group' column: map controls to "Control", others to "Treatment"
  if (!is.null(control_conditions)) {
    md$group <- ifelse(condition_vector %in% control_conditions, "Control", "Treat")
  } else {
    md$group <- condition_vector
  }

  # Add celltype_original to metadata
  md$celltype_original <- celltype_vector

  # Create celltype_condition column
  md[[new_ident_name]] <- paste(celltype_vector, condition_vector, sep = sep)

  # Create celltype_group column
  md$celltype_group <- paste(celltype_vector, md$group, sep = sep)

  # Update metadata in Seurat object
  obj@meta.data <- md

  # Set new Idents
  Seurat::Idents(obj) <- md[[new_ident_name]]

  return(obj)
}
