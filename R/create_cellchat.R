#' Create CellChat Objects from Integrated Spatial Seurat Object
#'
#' Splits an integrated Seurat object by `orig.ident`,
#' extracts normalized SCT data, loads spatial coordinates,
#' converts them to micron space using spot size, and constructs
#' a CellChat object for each sample.
#'
#' @param seu_integrated A Seurat object containing multiple spatial samples.
#' @param spatial_paths Named character vector of paths to spatial folders
#'   (containing `scalefactors_json.json`). Names must match `orig.ident`.
#' @param ident_col Metadata column used for splitting. Default `"orig.ident"`.
#' @param assay Assay to extract data from. Default `"SCT"`.
#' @param layer Layer to extract (Seurat v5). Default `"data"`.
#' @param spot_size Theoretical Visium spot size in microns. Default 65.
#'
#' @return A named list of CellChat objects.
#'
#' @importFrom Seurat SplitObject GetAssayData Idents GetTissueCoordinates
#' @importFrom jsonlite fromJSON
#' @importFrom CellChat createCellChat
#' @export
create_cellchat <- function(
    seu_integrated,
    spatial_paths,
    ident_col = "orig.ident",
    assay = "SCT",
    layer = "data",
    spot_size = 65
) {

  if (!ident_col %in% colnames(seu_integrated@meta.data)) {
    stop(paste0("Column '", ident_col, "' not found in meta.data."))
  }

  seu_list <- Seurat::SplitObject(seu_integrated, split.by = ident_col)
  sample_names <- names(seu_list)

  if (!all(sample_names %in% names(spatial_paths))) {
    stop("spatial_paths names must match orig.ident values.")
  }

  cellchat_list <- lapply(sample_names, function(s) {

    seu <- seu_list[[s]]

    # 1. Normalized data
    data.input <- Seurat::GetAssayData(object = seu, assay = assay, layer = layer)

    # 2. Metadata
    meta <- data.frame(
      labels = Seurat::Idents(seu),
      samples = s
    )

    meta$labels <- factor(meta$labels, levels = levels(Seurat::Idents(seu)))
    meta$samples <- factor(meta$samples, levels = s)

    # 3. Spatial coordinates
    spatial.locs <- Seurat::GetTissueCoordinates(
      object = seu,
      scale = NULL,
      cols = c("imagerow", "imagecol")
    )

    if (ncol(spatial.locs) > 2) spatial.locs <- spatial.locs[, 1:2]

    # 4. Load scale factors
    sf <- jsonlite::fromJSON(file.path(spatial_paths[[s]], "scalefactors_json.json"))

    if (is.null(sf$spot_diameter_fullres)) {
      stop(paste0("spot_diameter_fullres not found for sample ", s))
    }

    conversion <- spot_size / sf$spot_diameter_fullres
    spatial.locs <- spatial.locs * conversion

    scale.factors <- list(
      spot.diameter = spot_size,
      spot = spot_size
    )

    stopifnot(ncol(data.input) == nrow(spatial.locs))

    # 5. Create CellChat object
    cc <- CellChat::createCellChat(
      object = data.input,
      meta = meta,
      group.by = "labels",
      datatype = "spatial",
      coordinates = spatial.locs,
      scale.factors = scale.factors
    )

    # 6. Clean labels
    clean <- function(x) sub("_(Control|Treat)[0-9]*$", "", x)

    levels(cc@idents) <- clean(levels(cc@idents))
    cc@meta$labels <- factor(clean(cc@meta$labels), levels = levels(cc@idents))

    # 7. Order clusters by numeric prefix before
    ids <- levels(cc@idents)

    cluster_num <- suppressWarnings(as.numeric(sub("^([0-9]+).*", "\\1", ids)))

    ord <- order(cluster_num)

    ordered_levels <- ids[ord]

    cc@idents <- factor(cc@idents, levels = ordered_levels)
    cc@meta$labels <- factor(cc@meta$labels, levels = ordered_levels)

    return(cc)

  })

  names(cellchat_list) <- sample_names
  return(cellchat_list)
}
