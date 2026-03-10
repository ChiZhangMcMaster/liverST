#' Prepare CellChat Objects for Downstream Inference
#'
#' Assigns the Secreted Signaling ligand–receptor database
#' and performs preprocessing (subset data, identify overexpressed
#' genes and interactions) for a list of CellChat objects.
#'
#' @param cellchat_list A named list of CellChat objects.
#' @param species Character string: "mouse" (default) or "human".
#'
#' @return A named list of processed CellChat objects.
#'
#' @importFrom CellChat subsetDB subsetData identifyOverExpressedGenes identifyOverExpressedInteractions
#' @export
prepare_cellchat_calculations <- function(
    cellchat_list,
    species = "mouse"
) {

  if (!is.list(cellchat_list) || is.null(names(cellchat_list))) {
    stop("cellchat_list must be a named list of CellChat objects.")
  }

  #########################################
  # 1. Select database
  #########################################

  if (species == "mouse") {
    CellChatDB <- CellChat::CellChatDB.mouse
  } else if (species == "human") {
    CellChatDB <- CellChat::CellChatDB.human
  } else {
    stop("species must be 'mouse' or 'human'")
  }

  CellChatDB.use <- CellChat::subsetDB(
    CellChatDB,
    search = "Secreted Signaling",
    key = "annotation"
  )

  #########################################
  # 2. Process each object
  #########################################

  processed_list <- lapply(names(cellchat_list), function(nm) {

    obj <- cellchat_list[[nm]]

    # Assign database
    obj@DB <- CellChatDB.use

    # Remove lowly expressed genes
    obj <- CellChat::subsetData(obj)

    # Identify overexpressed genes
    obj <- CellChat::identifyOverExpressedGenes(obj)

    # Identify overexpressed interactions
    obj <- CellChat::identifyOverExpressedInteractions(obj)

    # Ensure spatial coordinates are matrix
    if (!is.null(obj@images$coordinates)) {
      obj@images$coordinates <- as.matrix(obj@images$coordinates)
    }

    return(obj)
  })

  names(processed_list) <- names(cellchat_list)

  return(processed_list)
}
