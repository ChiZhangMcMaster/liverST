#' Plot and save SpatialDimPlot with consistent cluster colors (preserves cluster 0)
#'
#' @param obj Seurat object
#' @param label Logical; whether to label regions
#' @param label.size Label size (default 3)
#' @param file File name for SVG output (without extension)
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param font Font family for SVG (default "Arial")
#' @param bg Background color (default "transparent")
#' @export
plot_spatial_svg <- function(obj,
                             label = TRUE,
                             label.size = 3,
                             file = "SpatialDimPlot",
                             width = 10,
                             height = 7,
                             font = "Arial",
                             bg = "transparent") {

  # --- Step 1: Create an unlabeled plot to extract cluster colors ---
  p_tmp <- Seurat::SpatialDimPlot(
    obj,
    label = FALSE,
    label.size = label.size
  )

  # Extract ggplot data
  plot_data <- ggplot2::ggplot_build(p_tmp)$data[[1]]
  cluster_colors_df <- unique(plot_data[, c("group", "fill")])

  # Convert cluster names to numeric, preserving 0
  cluster_colors_df$group <- as.numeric(as.character(cluster_colors_df$group))

  # Ensure cluster 0 is first and sort remaining clusters
  cluster_levels <- sort(unique(cluster_colors_df$group))
  cluster_colors_df <- cluster_colors_df[match(cluster_levels, cluster_colors_df$group), ]

  # Named vector for consistent colors
  cluster_colors_refined <- setNames(cluster_colors_df$fill, cluster_levels)

  # --- Step 2: Create final plot with consistent colors ---
  p <- Seurat::SpatialDimPlot(
    obj,
    label = label,
    label.size = label.size
  ) & ggplot2::scale_fill_manual(values = cluster_colors_refined)

  # --- Step 3: Save as SVG ---
  export::graph2svg(
    x = p,
    file = paste0(file, ".svg"),
    font = font,
    cairo = TRUE,
    width = width,
    height = height,
    bg = bg
  )

  return(p)
}
