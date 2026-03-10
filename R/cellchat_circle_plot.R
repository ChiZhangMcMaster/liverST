#' Plot aggregated or pathway-specific CellChat communication networks
#'
#' @param cellchat_list Named list of CellChat objects
#' @param save_dir Directory to save plots
#' @param skip_existing Logical; skip plots that already exist
#' @param file_format Output format: "tif", "png", or "svg"
#' @param pathway Optional signaling pathway. If NULL, plots aggregated networks
#' @param width Plot width (pixels for raster, inches for svg)
#' @param height Plot height (pixels for raster, inches for svg)
#' @param res Resolution for raster formats (dpi)
#' @return Invisibly returns input list
#' @importFrom CellChat netVisual_circle netVisual_aggregate
#' @importFrom graphics par
#' @importFrom grDevices tiff png svg dev.off
#' @export
cellchat_circle_plot <- function(
    cellchat_list,
    save_dir,
    skip_existing = FALSE,
    file_format = c("tif", "png", "svg"),
    pathway = NULL,
    width = NULL,
    height = NULL,
    res = 200
) {

  file_format <- match.arg(file_format)

  if (is.null(width))  width  <- ifelse(file_format == "svg", 12, 2400)
  if (is.null(height)) height <- ifelse(file_format == "svg", 6, 1200)

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  if (!is.list(cellchat_list) || is.null(names(cellchat_list))) {
    stop("cellchat_list must be a named list.")
  }

  scale_vertex <- function(x, min_size = 3, max_size = 10) {
    if (length(x) == 1) return(max_size)
    scaled <- (x - min(x)) / (max(x) - min(x) + 1e-6)
    min_size + scaled * (max_size - min_size)
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
        "_circle_network.", file_format
      )
    )

    if (skip_existing && file.exists(outfile)) {
      message("✔ Plot already exists. Skipping.")
      next
    }

    # open device
    if (file_format == "tif") {
      grDevices::tiff(outfile, width = width, height = height, res = res)
    } else if (file_format == "png") {
      grDevices::png(outfile, width = width, height = height, res = res)
    } else {
      grDevices::svg(outfile, width = width, height = height)
    }

    graphics::par(
      mfrow = if (is.null(pathway)) c(1,2) else c(1,1),
      mar = c(6,6,4,2),
      xpd = TRUE
    )

    if (!is.null(pathway)) {

      CellChat::netVisual_aggregate(
        obj,
        signaling = pathway,
        layout = "circle"
      )

      graphics::title(main = nm)

    } else {

      if (is.null(obj@net$count) || is.null(obj@net$weight)) {
        message("⚠ Network not aggregated yet. Skipping.")
        grDevices::dev.off()
        next
      }

      count_mat  <- obj@net$count
      weight_mat <- obj@net$weight

      vertex_count  <- scale_vertex(rowSums(count_mat))
      vertex_weight <- scale_vertex(rowSums(weight_mat))

      CellChat::netVisual_circle(
        count_mat,
        vertex.weight = vertex_count,
        weight.scale = TRUE,
        label.edge = FALSE,
        title.name = paste0(nm, ": Number of interactions")
      )

      CellChat::netVisual_circle(
        weight_mat,
        vertex.weight = vertex_weight,
        weight.scale = TRUE,
        label.edge = FALSE,
        title.name = paste0(nm, ": Interaction weights/strength")
      )
    }

    grDevices::dev.off()
    message("✔ Saved: ", outfile)
  }

  invisible(cellchat_list)
}
