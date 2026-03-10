#' Run pseudobulk differential expression analysis for a single cell type
#'
#' Automatically handles Seurat renaming of cell types (leading 'g') and uses all replicates.
#'
#' @param pseudo_obj Seurat pseudobulk object (output of perform_pseudobulk)
#' @param celltype_name Character. Name of the cell type (e.g., "Macrophages")
#' @param assay_counts Character. Assay to use for counts (default = "RNA")
#' @param lfc_threshold Numeric. Log2FC cutoff for labeling DEGs (default = 0.5)
#' @param padj_threshold Numeric. Adjusted p-value cutoff (default = 0.05)
#' @param control_pattern Pattern to detect control replicates (default = "Control")
#' @param treat_pattern Pattern to detect treated replicates (default = "Treat")
#'
#' @return Data frame with DE results (avg_log2FC, p_val, p_val_adj, Group)
run_pseudobulk_DE <- function(
    pseudo_obj,
    celltype_name,
    assay_counts = "Spatial",
    lfc_threshold = 0.5,
    padj_threshold = 0.05,
    condition_col = "group"
) {

  # ---- Check required metadata ----
  if(!"pseudobulk_id" %in% colnames(pseudo_obj@meta.data))
    stop("pseudobulk_id column not found.")

  if(!condition_col %in% colnames(pseudo_obj@meta.data))
    stop("Condition column not found.")

  meta <- pseudo_obj@meta.data

  # ---- Subset samples for this cell type (match anywhere in pseudobulk_id) ----
  keep_samples <- grepl(celltype_name, meta$pseudobulk_id)

  if(sum(keep_samples) == 0)
    stop("No pseudobulk samples found for cell type: ", celltype_name)

  sub_obj <- subset(pseudo_obj, cells = rownames(meta)[keep_samples])

  meta_sub <- sub_obj@meta.data
  counts <- Seurat::GetAssayData(sub_obj, assay = assay_counts, slot = "counts")

  # ---- Replicate sanity check ----
  replicate_table <- table(meta_sub[[condition_col]])
  print(replicate_table)

  if(any(replicate_table < 2))
    warning("Fewer than 2 replicates in at least one condition; DE results may be underpowered.")

  # ---- Ensure factor ----
  meta_sub[[condition_col]] <- factor(meta_sub[[condition_col]])

  # ---- Run DESeq2 ----
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta_sub,
    design = stats::as.formula(paste("~", condition_col))
  )

  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds)
  res <- as.data.frame(res)
  res$gene <- rownames(res)

  # ---- Standardize column naming ----
  res$avg_log2FC <- res$log2FoldChange

  # ---- Remove NA genes ----
  res <- res[!is.na(res$padj) & !is.na(res$avg_log2FC), ]

  # ---- DEG labeling ----
  res$Group <- "normal"
  res$Group[res$padj < padj_threshold & res$avg_log2FC > lfc_threshold] <- "up"
  res$Group[res$padj < padj_threshold & res$avg_log2FC < -lfc_threshold] <- "down"

  # ---- Final tidy output ----
  res <- res[, c("gene", "avg_log2FC", "pvalue", "padj", "Group")]
  colnames(res) <- c("gene", "avg_log2FC", "p_val", "p_val_adj", "Group")

  return(res)
}
