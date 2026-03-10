#' Plot spatial CellChat signaling networks
#'
#' Visualizes pathway-specific signaling on spatial coordinates using
#' `netVisual_aggregate(..., layout = "spatial")`.
#'
#' @param cellchat_list Named list of CellChat objects
#' @param pathway Signaling pathway(s) to visualize
#' @param save_dir Directory to save plots
#' @param sample.use Spatial sample name inside the CellChat object
#' @param skip_existing Logical; skip plots that already exist
#' @param file_format Output format: "tif", "png", or "svg"
#' @param width Plot width (pixels for raster, inches for svg)
#' @param height Plot height (pixels for raster, inches for svg)
#' @param res Resolution for raster formats (dpi)
#'
#' @return Invisibly returns input list
#'
#' @importFrom CellChat netVisual_aggregate
#' @importFrom graphics par
#' @importFrom grDevices tiff png svg dev.off
#' @export
cellchat_spatial_plot <- function(
    cellchat_list,
    pathway,
    file_format = "tif",
    save_dir = "cellchat_spatial_plots",
    width = 8,
    height = 8,
    res = 300
) {

  if (!dir.exists(save_dir))
    dir.create(save_dir, recursive = TRUE)

  if (is.null(names(cellchat_list)))
    names(cellchat_list) <- paste0("sample_", seq_along(cellchat_list))

  for (nm in names(cellchat_list)) {

    obj <- cellchat_list[[nm]]

    message("\n----------------------------------")
    message("Plotting spatial network: ", nm)
    message("----------------------------------")

    #################################
    # verify pathway exists
    #################################

    if (!pathway %in% obj@netP$pathways) {
      message("Skipping ", nm, ": pathway not present")
      next
    }

    #################################
    # output file
    #################################

    file_path <- file.path(
      save_dir,
      paste0(nm, "_", pathway, "_spatial.", file_format)
    )

    #################################
    # open graphics device
    #################################

    if (file_format == "pdf") {
      pdf(file_path, width = width, height = height)
    } else if (file_format == "png") {
      png(file_path, width = width, height = height, units = "in", res = res)
    } else if (file_format %in% c("tif","tiff")) {
      tiff(file_path, width = width, height = height, units = "in", res = res)
    } else {
      stop("Unsupported format")
    }

    par(mfrow = c(1,1))

    #################################
    # plot spatial signaling
    #################################

    p <- CellChat::netVisual_aggregate(
      obj,
      signaling = pathway,
      layout = "spatial",
      edge.width.max = 2,
      alpha.image = 0.2,
      vertex.size.max = 3,
      vertex.label.cex = 0
    )

    print(p)

    dev.off()

    message("Saved: ", file_path)
  }

  message("\nFinished spatial plotting.")
}
