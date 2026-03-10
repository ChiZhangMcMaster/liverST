#' Generate Basic Volcano Plot with Gene Labels
#'
#' @param deg_table DEG results data frame
#' @param logfc_col Column name for log2 fold change
#' @param padj_col Column name for adjusted p-value
#' @param logfc_threshold Log2FC cutoff for coloring
#' @param padj_threshold Adjusted p-value cutoff for coloring
#' @param label_logfc_threshold Log2FC cutoff for labeling genes
#' @param label_logq_threshold -log10(FDR) cutoff for labeling genes
#' @param top_n Number of top genes (by |log2FC|) to label
#' @param save_path Optional file path to save SVG (without extension)
#'
#' @return ggplot object
#' @export
plot_basic_volcano <- function(
    deg_table,
    logfc_col = "avg_log2FC",
    padj_col = "p_val_adj",
    logfc_threshold = 0.5,
    padj_threshold = 0.01,
    label_logfc_threshold = 2,
    label_logq_threshold = 3,
    top_n = 20,
    save_path = NULL
) {

  if (!requireNamespace("ggpubr", quietly = TRUE) ||
      !requireNamespace("ggthemes", quietly = TRUE)) {
    stop("Please install ggpubr and ggthemes.")
  }

  df <- as.data.frame(deg_table)

  if (!all(c(logfc_col, padj_col) %in% colnames(df))) {
    stop("Specified columns not found in deg_table.")
  }

  df$log2FC <- df[[logfc_col]]
  df$FDR <- df[[padj_col]]
  df$logQ <- -log10(df$FDR)
  df$ID <- rownames(df)

  # ---- Classify genes ----
  df$Group <- "normal"
  df$Group[df$FDR < padj_threshold & df$log2FC > logfc_threshold] <- "up"
  df$Group[df$FDR < padj_threshold & df$log2FC < -logfc_threshold] <- "down"

  # ---- Select genes for labeling ----
  df_sig <- df[
    df$logQ > label_logq_threshold &
      abs(df$log2FC) > label_logfc_threshold,
  ]

  if (nrow(df_sig) > 0) {
    df_sig <- df_sig[order(abs(df_sig$log2FC), decreasing = TRUE), ]
    df_sig <- head(df_sig, top_n)
    df$Label <- ifelse(df$ID %in% df_sig$ID, df$ID, "")
  } else {
    df$Label <- ""
  }

  # ---- Plot ----
  p <- ggpubr::ggscatter(
    df,
    x = "log2FC",
    y = "logQ",
    color = "Group",
    palette = c("down" = "#00BA38",
                "normal" = "#BBBBBB",
                "up" = "#F8766D"),
    size = 2,
    label = df$Label,
    font.label = 8,
    repel = TRUE,
    xlab = "log2FC",
    ylab = "-log10(FDR)"
  ) +
    ggthemes::theme_base() +
    ggplot2::geom_hline(
      yintercept = label_logq_threshold,
      linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      xintercept = c(-logfc_threshold, logfc_threshold),
      linetype = "dashed"
    )

  # ---- Optional save ----
  if (!is.null(save_path)) {
    ggplot2::ggsave(
      filename = paste0(save_path, ".svg"),
      plot = p,
      width = 6.5,
      height = 3.5,
      bg = "transparent"
    )
  }

  return(p)
}
