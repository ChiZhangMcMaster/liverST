#' Annotate clusters using GPT and embed cell types into Seurat object
#'
#' @param obj Seurat object
#' @param markers_df Marker dataframe (e.g. output from FindAllMarkers)
#' @param tissuename Character. Tissue description (default = "mouse liver")
#' @param model Character. GPT model name (default = "gpt-4")
#' @param output_csv Character. Output CSV filename
#'
#' @return Seurat object with renamed identities and celltype column
#'
#' @importFrom Seurat RenameIdents Idents
annotate_clusters_gpt <- function(obj,
                                  markers_df,
                                  tissuename = "mouse liver",
                                  model = "gpt-4",
                                  output_csv = "Annotation.csv") {

  message("Preparing marker list for GPT...")
  markers_list <- split(markers_df$gene, markers_df$cluster)

  message("Running GPT cell type annotation...")
  res <- GPTCelltype::gptcelltype(
    markers_list,
    tissuename = tissuename,
    model = model
  )

  # Convert named vector to data.frame if needed
  if (is.atomic(res) && !is.list(res)) {
    res <- data.frame(
      cluster = names(res),
      cell_type = as.character(res),
      stringsAsFactors = FALSE
    )
  }

  # Add numeric suffix only if there are multiple clusters with the same cell type
  res$cell_type <- ave(res$cell_type, res$cell_type, FUN = function(x) {
    if(length(x) == 1) {
      x  # keep original if only one
    } else {
      paste0(x, " ", seq_along(x))  # add suffix if multiple
    }
  })

  # Prepend cluster numbers
  res$cell_type <- paste0(res$cluster, ": ", res$cell_type)

  # Save to CSV
  write.csv(res, output_csv, row.names = FALSE)

  message("Renaming cluster identities...")
  celltype_vector <- res$cell_type
  names(celltype_vector) <- as.character(res$cluster)

  obj <- Seurat::RenameIdents(obj, celltype_vector)

  # Add celltype column into metadata
  obj$celltype <- Seurat::Idents(obj)

  message("Annotation complete.")
  return(obj)
}
