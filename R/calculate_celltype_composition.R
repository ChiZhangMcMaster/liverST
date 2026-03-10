#' Prepare average percent cell type composition across groups
#'
#' Computes the average number of cells per sample within each group,
#' converts to percent within each cell type, and orders cell types
#' by variability across groups.
#'
#' @param seurat_obj A Seurat object
#' @param selected_celltypes Optional character vector of cell types to keep (default = NULL)
#' @param sample_col Metadata column indicating sample/replicate ID (default = "orig.ident")
#' @param celltype_col Metadata column indicating cell type (default = "celltype")
#' @param group_col Metadata column indicating experimental group (default = "group")
#'
#' @return A data frame containing group, cell type, average cell count,
#'         and percent composition (ordered by variability)
calculate_celltype_composition <- function(
    seurat_obj,
    selected_celltypes = NULL,
    sample_col = "orig.ident",
    celltype_col = "celltype",
    group_col = "group"
) {

  md <- seurat_obj@meta.data |> as.data.frame()

  # Remove numeric-only cluster labels
  md <- md |>
    dplyr::filter(!grepl("^[0-9]+$", .data[[celltype_col]]))

  # Optional subset
  if (!is.null(selected_celltypes)) {
    md <- md |>
      dplyr::filter(.data[[celltype_col]] %in% selected_celltypes)
  }

  # Count cells per sample × cell type
  df_counts <- md |>
    dplyr::group_by(.data[[sample_col]],
                    .data[[group_col]],
                    .data[[celltype_col]]) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

  # Compute total cells per sample
  df_totals <- md |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::summarise(total_cells = dplyr::n(), .groups = "drop")

  # Join totals
  df_counts <- df_counts |>
    dplyr::left_join(df_totals, by = sample_col)

  # Compute fraction within each sample
  df_counts <- df_counts |>
    dplyr::mutate(Percent = n_cells / total_cells * 100)

  # Average percentages across samples within group
  df_avg <- df_counts |>
    dplyr::group_by(.data[[group_col]],
                    .data[[celltype_col]]) |>
    dplyr::summarise(Percent = mean(Percent), .groups = "drop")

  # ----------------------------
  # Row-normalize so each cell type sums to 100%
  # ----------------------------
  df_avg <- df_avg |>
    dplyr::group_by(.data[[celltype_col]]) |>
    dplyr::mutate(NormPercent = Percent / sum(Percent) * 100) |>
    dplyr::ungroup()

  # Order cell types by variability across groups (optional)
  ordering <- df_avg |>
    dplyr::group_by(.data[[celltype_col]]) |>
    dplyr::summarise(range = max(NormPercent) - min(NormPercent), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(range)) |>
    dplyr::pull(.data[[celltype_col]])

  df_avg[[celltype_col]] <- factor(df_avg[[celltype_col]],
                                   levels = ordering)

  return(df_avg)
}
