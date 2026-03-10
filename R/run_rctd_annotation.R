#' Run RCTD deconvolution using a scRNA-seq reference
#'
#' @param spatial_obj Seurat spatial object
#' @param reference_obj Seurat scRNA reference object
#' @param ref_celltype_col Column name in reference metadata containing cell types
#' @param ref_assay Assay name in reference (default = "RNA")
#' @param spatial_assay Assay name in spatial object (default = "SCT")
#' @param min_cells Minimum number of cells per reference cell type (default = 25)
#' @param max_cores Number of cores for RCTD (default = 4)
#' @param doublet_mode RCTD doublet mode ("doublet" or "full"; default = "doublet")
#'
#' @return Seurat object with RCTD results added to metadata
#' @export
run_rctd_annotation <- function(
    spatial_obj,
    reference_obj,
    ref_celltype_col = "celltype",
    ref_assay = "RNA",
    spatial_assay = "SCT",
    min_cells = 25,
    max_cores = 4,
    doublet_mode = "doublet"
) {

  # ---- Dependency check ----
  if (!requireNamespace("spacexr", quietly = TRUE)) {
    stop(
      "Package 'spacexr' is required for RCTD annotation.\n",
      "Install with:\n",
      "  remotes::install_github('dmcable/spacexr')"
    )
  }

  message("Preparing reference object...")

  # Extract reference data
  Seurat::Idents(reference_obj) <- reference_obj[[ref_celltype_col, drop = TRUE]]

  counts <- Seurat::GetAssayData(reference_obj, assay = ref_assay, slot = "counts")
  cluster <- as.factor(reference_obj[[ref_celltype_col, drop = TRUE]])
  nUMI <- reference_obj$nCount_RNA

  levels(cluster) <- gsub("/", "-", levels(cluster))
  cluster <- droplevels(cluster)

  reference <- spacexr::Reference(counts, cluster, nUMI)

  # Filter rare cell types
  celltype_counts <- table(reference@cell_types)
  valid_celltypes <- names(celltype_counts[celltype_counts > min_cells])

  selected_cells <- names(reference@cell_types)[
    reference@cell_types %in% valid_celltypes
  ]

  reference_subset <- methods::new(
    "Reference",
    counts = reference@counts[, selected_cells],
    nUMI = reference@nUMI[selected_cells],
    cell_types = droplevels(reference@cell_types[selected_cells])
  )

  message("Preparing spatial query...")

  counts_spatial <- Seurat::GetAssayData(
    spatial_obj,
    assay = spatial_assay,
    slot = "counts"
  )

  coords <- Seurat::GetTissueCoordinates(spatial_obj)
  coords <- coords[colnames(counts_spatial), 1:2]

  query <- spacexr::SpatialRNA(
    coords,
    counts_spatial,
    colSums(counts_spatial)
  )

  message("Running RCTD...")

  rctd_obj <- spacexr::create.RCTD(
    query,
    reference_subset,
    max_cores = max_cores
  )

  rctd_obj <- spacexr::run.RCTD(
    rctd_obj,
    doublet_mode = doublet_mode
  )

  message("Adding results to Seurat object...")

  spatial_obj <- Seurat::AddMetaData(
    spatial_obj,
    metadata = rctd_obj@results$results_df
  )

  message("RCTD annotation complete.")

  return(spatial_obj)
}
