#' Plot ligand-receptor contributions for CellChat objects (tiff/png/svg)
#'
#' Computes the contribution of each ligand-receptor pair to a signaling pathway
#' using `netAnalysis_contribution()` and generates a plot for each sample.
#'
#' @param cellchat_list Named list of CellChat objects
#' @param pathway Signaling pathway(s) to visualize
#' @param save_dir Directory to save plots
#' @param file_format Output format: "tif", "png", or "svg"
#' @param width Plot width (pixels for raster, inches for svg)
#' @param height Plot height (pixels for raster, inches for svg)
#' @param res Resolution for raster formats (dpi)
#' @param skip_existing Logical; skip plots that already exist
#'
#' @return Invisibly returns a named list of contribution tables
#' @importFrom CellChat netAnalysis_contribution
#' @importFrom grDevices tiff png svg dev.off
#' @importFrom graphics par title
#' @export
cellchat_LR_contribution_plot <- function(
    cellchat_list,
    pathway,
    save_dir = "cellchat_LR_plots",
    file_format = c("tif", "png", "svg"),
    width = NULL,
    height = NULL,
    res = 300,
    skip_existing = TRUE
) {

  file_format <- match.arg(file_format)

  # default sizes like your heatmap function
  if (is.null(width))  width  <- ifelse(file_format == "svg", 8, 1600)
  if (is.null(height)) height <- ifelse(file_format == "svg", 8, 1600)

  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  if (is.null(names(cellchat_list))) names(cellchat_list) <- paste0("sample_", seq_along(cellchat_list))

  contrib_list <- list()

  for (nm in names(cellchat_list)) {
    obj <- cellchat_list[[nm]]

    message("\n----------------------------------")
    message("Processing sample: ", nm)
    message("----------------------------------")

    # check pathway exists
    if (is.null(obj@netP$pathways) || !pathway %in% obj@netP$pathways) {
      message("Skipping ", nm, ": pathway not present")
      next
    }

    # output file
    file_path <- file.path(save_dir, paste0(nm, "_", pathway, "_LR_contribution.", file_format))
    if (skip_existing && file.exists(file_path)) {
      message("Skipping ", nm, ": file already exists")
      next
    }

    # compute contribution
    contrib <- CellChat::netAnalysis_contribution(obj, signaling = pathway)
    contrib_list[[nm]] <- contrib

    # open graphics device
    if (file_format == "svg") {
      svg(file_path, width = width, height = height)
    } else if (file_format == "png") {
      png(file_path, width = width, height = height, units = "px", res = res)
    } else if (file_format %in% c("tif","tiff")) {
      tiff(file_path, width = width, height = height, units = "px", res = res)
    } else {
      stop("Unsupported file format: ", file_format)
    }

    par(mar = c(5, 12, 4, 2)) # more left margin for LR names

    # simple horizontal barplot
    contrib_values <- contrib$contribution  # numeric values
    names(contrib_values) <- rownames(contrib) # LR pair names

    barplot(
      sort(contrib_values, decreasing = TRUE),
      horiz = TRUE,
      las = 1,
      col = "steelblue",
      main = paste(nm, pathway),
      xlab = "Contribution",
      cex.names = 0.8
    )

    dev.off()
    message("Saved plot: ", file_path)
  }

  message("\nFinished plotting LR contributions.")
  invisible(contrib_list)
}
