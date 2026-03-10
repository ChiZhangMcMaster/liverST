#' Perform pseudobulk aggregation from a Seurat object (raw counts)
#'
#' @param seurat_obj Seurat object
#' @param condition_col Metadata column for condition (default = "group")
#' @param replicate_col Metadata column for replicate/sample (default = "orig.ident")
#' @param celltype_col Metadata column for cell type (default = "celltype")
#' @param assay Assay to use for pseudobulk (default = "RNA" raw counts)
#' @param return_seurat Logical; return Seurat object (default = TRUE)
#'
#' @return Pseudobulk Seurat object with `celltype_condition_replicate` as active identity
perform_pseudobulk <- function(
    seurat_obj,
    condition_col = "group",
    replicate_col = "orig.ident",
    celltype_col = "celltype",
    assay = "Spatial",
    return_seurat = TRUE
) {

  obj <- seurat_obj

  # ---- Check metadata ----
  required_cols <- c(condition_col, replicate_col, celltype_col)
  missing_cols <- required_cols[!required_cols %in% colnames(obj@meta.data)]
  if(length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))

  # ---- Clean celltype names ----
  obj$celltype_clean <- gsub("^[0-9]+_", "", gsub(" ", "_", obj[[celltype_col]][,1]))

  # ---- Create pseudobulk ID ----
  obj$pseudobulk_id <- paste(obj$celltype_clean,
                             obj[[condition_col]][,1],
                             obj[[replicate_col]][,1],
                             sep = "_")

  # ---- Save metadata mapping before aggregation ----
  meta_mapping <- unique(obj@meta.data[, c("pseudobulk_id", celltype_col, condition_col, replicate_col)])
  rownames(meta_mapping) <- meta_mapping$pseudobulk_id

  # ---- Aggregate counts ----
  pseudo_obj <- Seurat::AggregateExpression(
    obj,
    assays = assay,
    slot = "counts",
    return.seurat = TRUE,
    group.by = "pseudobulk_id"
  )

  # ---- Clean Seurat "g" prefix if present in metadata ----
  pseudo_obj$pseudobulk_id <- sub("^g", "", pseudo_obj$pseudobulk_id)

  # ---- Attach metadata mapping ----
  pseudo_obj <- Seurat::AddMetaData(pseudo_obj, meta_mapping)

  # ---- Set identity ----
  Seurat::Idents(pseudo_obj) <- "pseudobulk_id"

  cat("Pseudobulk samples created:\n")
  print(table(pseudo_obj$pseudobulk_id))

  return(pseudo_obj)
}
