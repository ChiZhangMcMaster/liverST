#' Plot cell type annotation on UMAP and SpatialDimPlot with consistent colors
#'
#' Generates UMAP and SpatialDimPlot for Seurat objects using cell type annotations,
#' ensures consistent colors between plots, and optionally combines them.
#'
#' @param object Seurat object with cell type annotations in Idents(object) or "celltype" column.
#' @param label Logical; whether to label clusters (default = TRUE).
#' @param label.size Numeric; size of cluster labels in SpatialDimPlot (default = 3).
#' @param combine Logical; whether to combine UMAP and Spatial plots using patchwork (default = TRUE).
#'
#' @return Either a combined ggplot object (if combine = TRUE) or a list of plots (umap_plot, spatial_plot).
#' @export
plot_celltype_annotation <- function(
    object,
    label = TRUE,
    label.size = 3,
    combine = TRUE
) {

  # Use celltype identities if present
  if ("celltype" %in% colnames(object@meta.data)) {
    Seurat::Idents(object) <- object$celltype
  }

  # Extract celltype identities
  celltype <- Seurat::Idents(object)

  # Extract numeric prefix (before ":") for ordering
  numeric_prefix <- as.numeric(sub(":.*", "", celltype))

  # Create unique levels in numeric order
  all_clusters <- unique(celltype[order(numeric_prefix)])

  object$celltype <- factor(celltype, levels = all_clusters)

  # UMAP plot
  p1 <- Seurat::DimPlot(
    object,
    group.by = "celltype",
    label = FALSE
  )

  # Spatial reference to extract colors
  p_ref <- Seurat::SpatialDimPlot(
    object,
    group.by = "celltype",
    label = FALSE,
    combine = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p_ref)$data[[1]]
  cluster_colors_df <- unique(plot_data[, c("group", "fill")])

  # Assign consistent colors
  default_colors <- scales::hue_pal()(length(all_clusters))
  cluster_colors <- setNames(default_colors, all_clusters)
  cluster_colors[cluster_colors_df$group] <- cluster_colors_df$fill

  # Spatial plot
  p2 <- Seurat::SpatialDimPlot(
    object,
    group.by = "celltype",
    label = label,
    label.size = label.size
  ) +
    ggplot2::scale_fill_manual(values = cluster_colors)

  # Return combined or separate plots
  if (combine) {
    return(patchwork::wrap_plots(p1, p2))
  } else {
    return(list(umap_plot = p1, spatial_plot = p2))
  }
}
