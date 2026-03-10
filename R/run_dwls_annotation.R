#' Run DWLS deconvolution using a scRNA reference and return annotated Seurat object
#'
#' @param spatial_obj Seurat spatial object (to annotate)
#' @param reference_obj Seurat scRNA-seq reference object
#' @param ref_celltype_col Column name in reference containing cell types
#' @param spatial_cluster_col Column in spatial object used for clusters
#' @param spatial_assay Assay name in spatial object (default = "Spatial")
#' @param ref_assay Assay name in reference object (default = "SCT")
#' @param downsample Number of cells per cell type for reference (default = 200)
#' @param n_cell Number of cells for DWLS (default = 10)
#'
#' @return Seurat object with DWLS deconvolution results
#' @export
run_dwls_annotation <- function(
    spatial_obj,
    reference_obj,
    ref_celltype_col = "celltype",
    spatial_cluster_col = "seurat_clusters",
    spatial_assay = "Spatial",
    ref_assay = "SCT",
    downsample = 200,
    n_cell = 10
) {

  message("Preparing reference object...")

  # Set reference identities
  Seurat::Idents(reference_obj) <- reference_obj[[ref_celltype_col, drop = TRUE]]

  # Downsample reference
  reference_obj <- subset(reference_obj, downsample = downsample)

  # Convert reference to Giotto
  giotto_ref <- Giotto::seuratToGiottoV5(
    sobject = reference_obj,
    spatial_assay = ref_assay
  )

  giotto_ref <- Giotto::normalizeGiotto(giotto_ref)

  message("Finding marker genes...")
  markers <- Giotto::findMarkers_one_vs_all(
    gobject = giotto_ref,
    method = "scran",
    expression_values = "normalized",
    cluster_column = ref_celltype_col,
    min_feats = 3
  )

  top_markers <- markers[, head(.SD, 10), by = "cluster"]$feats

  message("Building signature matrix...")
  sign_matrix <- Giotto::makeSignMatrixDWLSfromMatrix(
    matrix = Giotto::getExpression(
      giotto_ref,
      values = "normalized",
      output = "matrix"
    ),
    cell_type = Giotto::pDataDT(giotto_ref)[[ref_celltype_col]],
    sign_gene = top_markers
  )

  message("Preparing spatial object...")

  spatial_obj$cell_type <- spatial_obj[[spatial_cluster_col, drop = TRUE]]

  giotto_spatial <- Giotto::seuratToGiottoV5(
    sobject = spatial_obj,
    spatial_assay = spatial_assay
  )

  message("Running DWLS deconvolution...")
  giotto_spatial <- Giotto::runDWLSDeconv(
    giotto_spatial,
    sign_matrix = sign_matrix,
    cluster_column = "cell_type",
    n_cell = n_cell,
    return_gobject = TRUE
  )

  message("Converting back to Seurat...")
  result_seurat <- Giotto::giottoToSeuratV5(giotto_spatial)

  message("DWLS annotation complete.")
  return(result_seurat)
}
