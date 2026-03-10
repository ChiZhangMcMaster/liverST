add_celltype_stats <- function(
    seurat_obj,
    selected_celltypes = NULL,
    sample_col = "orig.ident",
    celltype_col = "celltype",
    group_col = "group"
) {

  md <- seurat_obj@meta.data

  if (!is.null(selected_celltypes)) {
    md <- md[md[[celltype_col]] %in% selected_celltypes, ]
  }

  df_counts <- md |>
    dplyr::group_by(
      .data[[sample_col]],
      .data[[group_col]],
      .data[[celltype_col]]
    ) |>
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop")

  df_prop <- df_counts |>
    dplyr::group_by(.data[[sample_col]]) |>
    dplyr::mutate(prop = n_cells / sum(n_cells) * 100) |>
    dplyr::ungroup()

  groups <- unique(df_prop[[group_col]])

  if (length(groups) != 2) {
    stop("Statistical testing currently implemented for 2 groups only.")
  }

  df_stats <- df_prop |>
    dplyr::group_by(.data[[celltype_col]]) |>
    dplyr::summarise(
      p_value = stats::wilcox.test(
        prop ~ .data[[group_col]]
      )$p.value,
      .groups = "drop"
    )

  df_stats$p_adj <- stats::p.adjust(df_stats$p_value, method = "BH")

  return(df_stats)
}
