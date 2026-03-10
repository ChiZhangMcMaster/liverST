#' Plot CellChat heatmaps or pathway signaling role networks
#'
#' @param cellchat_list Named list of CellChat objects
#' @param save_dir Directory to save plots
#' @param measure Interaction metric: "count" or "weight"
#' @param pathway Optional signaling pathway. If NULL, plots interaction heatmap
#' @param skip_existing Logical; skip plots that already exist
#' @param file_format Output format: "tif", "png", or "svg"
#' @param width Plot width (pixels for raster, inches for svg)
#' @param height Plot height (pixels for raster, inches for svg)
#' @param res Resolution for raster formats (dpi)
#' @param color.heatmap Color palette for heatmap
#' @return Invisibly returns input list
#' @importFrom CellChat netVisual_heatmap netAnalysis_signalingRole_network
#' @importFrom grDevices tiff png svg dev.off
#' @export
cellchat_heatmap <- function(
    cellchat_list,
    save_dir,
    measure = c("count", "weight"),
    pathway = NULL,
    skip_existing = FALSE,
    file_format = c("tif", "png", "svg"),
    width = NULL,
    height = NULL,
    res = 200,
    color.heatmap = "Blues"
) {

  measure <- match.arg(measure)
  file_format <- match.arg(file_format)

  if (is.null(width))  width  <- ifelse(file_format == "svg", 8, 1600)
  if (is.null(height)) height <- ifelse(file_format == "svg", 8, 1600)

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  if (!is.list(cellchat_list) || is.null(names(cellchat_list))) {
    stop("cellchat_list must be a named list.")
  }

  for (nm in names(cellchat_list)) {

    message("\n----------------------------------")
    message("Plotting: ", nm)
    message("----------------------------------")

    obj <- cellchat_list[[nm]]

    outfile <- file.path(
      save_dir,
      paste0(
        "cellchat_", nm,
        if (!is.null(pathway)) paste0("_", pathway),
        if (is.null(pathway)) paste0("_heatmap_", measure),
        ".", file_format
      )
    )

    if (skip_existing && file.exists(outfile)) {
      message("✔ Plot already exists. Skipping.")
      next
    }

    # open graphics device
    if (file_format == "tif") {
      grDevices::tiff(outfile, width = width, height = height, res = res)
    } else if (file_format == "png") {
      grDevices::png(outfile, width = width, height = height, res = res)
    } else {
      grDevices::svg(outfile, width = width, height = height)
    }

    if (!is.null(pathway)) {

      p <- CellChat::netAnalysis_signalingRole_network(
        obj,
        signaling = pathway,
        width = 8,
        height = 2.5,
        font.size = 10
      )

      if (!is.null(p)) print(p)

    } else {

      p <- CellChat::netVisual_heatmap(
        obj,
        measure = measure,
        color.heatmap = color.heatmap
      )

      if (!is.null(p)) print(p)

    }

    grDevices::dev.off()
    message("✔ Saved: ", outfile)
  }

  invisible(cellchat_list)
}
