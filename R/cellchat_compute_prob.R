#' Compute Communication Probabilities (Crash-Safe & Resume-Safe)
#'
#' Runs `computeCommunProb()` for each CellChat object in a named list.
#' Each object is saved immediately after successful computation.
#' If a sample fails, processing continues for the remaining samples.
#' Existing output files can be skipped for resume-safe execution.
#'
#' @param cellchat_list A named list of CellChat objects.
#' @param save_dir Directory where processed objects will be saved.
#' @param type Method for averaging expression. Default "truncatedMean".
#' @param trim Trimming fraction for truncatedMean. Default 0.1.
#' @param distance.use Logical; whether to use spatial distance. Default FALSE.
#' @param reload Logical; reload object from disk after saving. Default TRUE.
#' @param skip_existing Logical; skip samples if output file already exists. Default TRUE.
#'
#' @return A named list of successfully processed CellChat objects.
#'
#' @importFrom CellChat computeCommunProb
#' @export
cellchat_compute_prob <- function(
    cellchat_list,
    save_dir,
    type = "truncatedMean",
    trim = 0.1,
    distance.use = FALSE,
    reload = TRUE,
    skip_existing = TRUE
) {

  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  if (!is.list(cellchat_list) || is.null(names(cellchat_list))) {
    stop("cellchat_list must be a named list of CellChat objects.")
  }

  results <- list()
  failed  <- character(0)

  for (nm in names(cellchat_list)) {

    message("\n======================================")
    message("Processing sample: ", nm)
    message("======================================")

    outfile <- file.path(
      save_dir,
      paste0("cellchat_", nm, "_after_CP.rds")
    )

    # Resume-safe skip
    if (skip_existing && file.exists(outfile)) {
      message("✔ File already exists. Loading from disk.")
      results[[nm]] <- readRDS(outfile)
      next
    }

    obj <- cellchat_list[[nm]]

    # Crash-safe wrapper
    tryCatch({

      obj <- CellChat::computeCommunProb(
        obj,
        type = type,
        trim = trim,
        distance.use = distance.use
      )

      saveRDS(obj, outfile)
      message("✔ Successfully saved: ", outfile)

      if (reload) {
        obj <- readRDS(outfile)
      }

      results[[nm]] <- obj

    }, error = function(e) {

      message("✖ ERROR in sample: ", nm)
      message("  → ", e$message)
      failed <<- c(failed, nm)

    })
  }

  if (length(failed) > 0) {
    message("\n⚠ The following samples failed:")
    message(paste(failed, collapse = ", "))
  } else {
    message("\n✔ All samples processed successfully.")
  }

  return(results)
}
