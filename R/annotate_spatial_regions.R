#' Interactively annotate spatial regions in a Seurat object
#'
#' @param seurat_obj A Seurat object with spatial data
#' @param image Name of spatial image (e.g. "condition2a")
#' @param region_names Character vector of region names to annotate
#' @param assay Assay to use for expression data (default = "SCT")
#'
#' @return Seurat object with new metadata column "region"
#' @export
annotate_spatial_regions <- function(
    seurat_obj,
    image,
    region_names = c("Region1"),
    assay = "SCT"
) {

  if (!requireNamespace("Seurat", quietly = TRUE) ||
      !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install Seurat and ggplot2.")
  }

  # Initialize region column
  seurat_obj$region <- "other"

  # Extract coordinates
  coords <- Seurat::GetTissueCoordinates(
    seurat_obj,
    scale = NULL,
    cols = c("imagerow", "imagecol")
  )

  # Rename columns for ggplot consistency
  colnames(coords) <- c("x", "y")

  for (region in region_names) {

    message("Draw region for: ", region)

    p <- ggplot2::ggplot(coords, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point() +
      ggplot2::theme_void()

    selected <- Seurat::CellSelector(p)

    if (length(selected) == 0) {
      warning("No spots selected for region: ", region)
      next
    }

    seurat_obj$region[selected] <- region
  }

  return(seurat_obj)
}
