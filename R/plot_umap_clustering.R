#' Plot UMAP and Spatial Clusters
#'
#' Generates UMAP and spatial cluster plots and optionally combines them.
#'
#' @param object A processed Seurat object.
#' @param label Logical; whether to label clusters (default = TRUE).
#' @param label.size Size of spatial labels (default = 3).
#' @param combine Logical; whether to combine plots using patchwork (default = TRUE).
#'
#' @return A ggplot object (combined or individual plots)
#' @export
plot_umap_clustering <- function(
    object,
    label = TRUE,
    label.size = 3,
    combine = TRUE
) {

  # Ensure consistent cluster factor levels
  clusters <- Seurat::Idents(object)
  all_clusters <- sort(unique(as.numeric(as.character(clusters))))
  object$seurat_clusters <- factor(clusters, levels = all_clusters)

  # UMAP
  p1 <- Seurat::DimPlot(
    object,
    reduction = "umap",
    group.by = "seurat_clusters",
    label = label
  )

  # Spatial reference to extract colors
  p_ref <- Seurat::SpatialDimPlot(
    object,
    group.by = "seurat_clusters",
    label = FALSE,
    combine = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p_ref)$data[[1]]
  cluster_colors_df <- unique(plot_data[, c("group", "fill")])

  default_colors <- scales::hue_pal()(length(all_clusters))
  cluster_colors <- setNames(default_colors, as.character(all_clusters))
  cluster_colors[cluster_colors_df$group] <- cluster_colors_df$fill

  # Spatial plot
  p2 <- Seurat::SpatialDimPlot(
    object,
    group.by = "seurat_clusters",
    label = label,
    label.size = label.size
  ) +
    ggplot2::scale_fill_manual(values = cluster_colors)

  if (combine) {
    return(patchwork::wrap_plots(p1, p2))
  } else {
    return(list(umap_plot = p1, spatial_plot = p2))
  }
}


