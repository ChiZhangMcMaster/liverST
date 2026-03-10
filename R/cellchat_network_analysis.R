#' Finalize CellChat Network Analysis
#'
#' Performs downstream CellChat analysis including communication filtering,
#' pathway-level communication inference, network aggregation, and centrality
#' computation. Each CellChat object in the input list is processed sequentially
#' and saved immediately after completion to allow crash recovery.
#'
#' @param cellchat_list Named list of CellChat objects.
#' @param save_dir Directory where processed objects will be saved.
#' @param min_cells Minimum number of cells required for a cell group to be
#' included in communication filtering. Default is 10.
#' @param slot_name Slot used for centrality calculation. Default `"netP"`.
#' @param skip_existing Logical; if TRUE, already processed samples are loaded
#' from disk instead of recomputed. Default TRUE.
#'
#' @return Named list of processed CellChat objects.
#'
#' @importFrom CellChat filterCommunication
#' @importFrom CellChat computeCommunProbPathway
#' @importFrom CellChat aggregateNet
#' @importFrom CellChat netAnalysis_computeCentrality
#' @export
cellchat_network_analysis <- function(
    cellchat_list,
    save_dir,
    min_cells = 10,
    slot_name = "netP",
    skip_existing = TRUE
) {

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  if (!is.list(cellchat_list) || is.null(names(cellchat_list))) {
    stop("cellchat_list must be a named list.")
  }

  results <- list()
  failed <- character(0)

  for (nm in names(cellchat_list)) {

    message("\n----------------------------------")
    message("Processing: ", nm)
    message("----------------------------------")

    outfile <- file.path(
      save_dir,
      paste0("cellchat_", nm, "_final.rds")
    )

    if (skip_existing && file.exists(outfile)) {
      message("✔ Found existing file. Loading.")
      results[[nm]] <- readRDS(outfile)
      next
    }

    obj <- cellchat_list[[nm]]

    tryCatch({

      #########################################
      # 4. Filter communications
      #########################################
      obj <- CellChat::filterCommunication(obj, min.cells = min_cells)

      #########################################
      # 5. Compute pathway-level communication
      #########################################
      obj <- CellChat::computeCommunProbPathway(obj)

      #########################################
      # 6. Aggregate network
      #########################################
      obj <- CellChat::aggregateNet(obj)

      #########################################
      # 7. Network centrality
      #########################################
      obj <- CellChat::netAnalysis_computeCentrality(
        obj,
        slot.name = slot_name
      )

      saveRDS(obj, outfile)
      message("✔ Saved: ", outfile)

      results[[nm]] <- obj

    }, error = function(e) {

      message("✖ Failed for sample: ", nm)
      message("  → ", e$message)

      failed <<- c(failed, nm)

    })
  }

  if (length(failed) > 0) {
    message("\n⚠ Failed samples:")
    message(paste(failed, collapse = ", "))
  } else {
    message("\n✔ All samples completed.")
  }

  return(results)
}
