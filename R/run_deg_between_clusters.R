#' Run differential expression between two identities in a Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param ident_1 Identity/group 1 (e.g., "1: Hepatocytes2")
#' @param ident_2 Identity/group 2 (e.g., "0: Hepatocytes1")
#' @param group_col Metadata column to use as identities (default = NULL,
#'   meaning current Idents are used)
#' @param assay Assay to use (default = NULL, uses active assay)
#' @param slot Slot to use for testing (default = "data")
#' @param logfc.threshold Minimum log fold change threshold (default = 0.25)
#' @param min.pct Minimum percent expression (default = 0.1)
#' @param test.use Statistical test to use (default = "wilcox")
#'
#' @return A data.frame of differential expression results
#' @export
run_deg_between_clusters <- function(
    seurat_obj,
    ident_1,
    ident_2,
    group_col = NULL,
    assay = NULL,
    slot = "data",
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox"
) {

  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  # Optionally set identities from metadata column
  if (!is.null(group_col)) {
    if (!group_col %in% colnames(seurat_obj@meta.data)) {
      stop("group_col not found in metadata.")
    }
    Seurat::Idents(seurat_obj) <- seurat_obj[[group_col, drop = TRUE]]
  }

  # Optionally set assay
  if (!is.null(assay)) {
    Seurat::DefaultAssay(seurat_obj) <- assay
  }

  message("Running differential expression: ",
          ident_1, " vs ", ident_2)

  deg_results <- Seurat::FindMarkers(
    object = seurat_obj,
    ident.1 = ident_1,
    ident.2 = ident_2,
    slot = slot,
    logfc.threshold = logfc.threshold,
    min.pct = min.pct,
    test.use = test.use
  )

  return(deg_results)
}
