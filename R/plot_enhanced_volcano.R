#' Create enhanced volcano plot from DEG results
#'
#' @param deg_table Data frame of differential expression results
#' @param logfc_col Column name for log2 fold change (default = "avg_log2FC")
#' @param padj_col Column name for adjusted p-value (default = "p_val_adj")
#' @param logfc_threshold Log2FC cutoff (default = 0.5)
#' @param padj_threshold Adjusted p-value cutoff (default = 0.05)
#' @param top_n Number of top genes to label (default = 20)
#' @param title Plot title (default = NULL)
#'
#' @return ggplot object
#' @export
plot_enhanced_volcano <- function(
    deg_table,
    logfc_col = "avg_log2FC",
    padj_col = "p_val_adj",
    logfc_threshold = 0.5,
    padj_threshold = 0.05,
    top_n = 20,
    title = NULL
) {

  if (!all(c(logfc_col, padj_col) %in% colnames(deg_table))) {
    stop("Specified logfc_col or padj_col not found in deg_table.")
  }

  df <- deg_table

  # Standardize column names internally (do not modify original)
  df$log2FoldChange <- df[[logfc_col]]
  df$padj <- df[[padj_col]]

  # Count significant genes
  num_up <- sum(df$padj < padj_threshold & df$log2FoldChange > logfc_threshold, na.rm = TRUE)
  num_down <- sum(df$padj < padj_threshold & df$log2FoldChange < -logfc_threshold, na.rm = TRUE)

  # Assign colors
  df$color <- dplyr::case_when(
    df$padj < padj_threshold & df$log2FoldChange > logfc_threshold ~ "Significant Up",
    df$padj < padj_threshold & df$log2FoldChange < -logfc_threshold ~ "Significant Down",
    TRUE ~ "Not Significant"
  )

  # Select top genes
  top_genes <- df[!is.na(df$padj) &
                    df$padj < padj_threshold &
                    abs(df$log2FoldChange) > logfc_threshold, ]

  if (nrow(top_genes) > 0) {
    top_genes$score <- -log10(top_genes$padj) * abs(top_genes$log2FoldChange)
    top_genes <- top_genes[order(-top_genes$score), ]
    top_genes <- head(top_genes, top_n)
    top_genes$gene <- rownames(top_genes)
  }

  volcano_plot <- ggplot2::ggplot(df,
                                  ggplot2::aes(x = log2FoldChange, y = -log10(padj))
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = color),
      size = 1.5, alpha = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = c(
        "Significant Up" = "red",
        "Significant Down" = "blue",
        "Not Significant" = "grey50"
      )
    ) +
    ggplot2::geom_vline(
      xintercept = c(-logfc_threshold, logfc_threshold),
      color = "grey", linewidth = 0.3
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_threshold),
      color = "grey", linewidth = 0.3
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      caption = paste("Total significant genes:", num_up + num_down)
    ) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text = ggplot2::element_text(color = "black"),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.line = ggplot2::element_line(linewidth = 0.4, color = "black"),
      plot.caption = ggplot2::element_text(face = "italic", size = 10)
    ) +
    ggplot2::annotate(
      "text", x = Inf, y = Inf,
      label = paste("Up:", num_up),
      hjust = 1.1, vjust = 2,
      color = "red", size = 4, fontface = "bold"
    ) +
    ggplot2::annotate(
      "text", x = -Inf, y = Inf,
      label = paste("Down:", num_down),
      hjust = -0.1, vjust = 2,
      color = "blue", size = 4, fontface = "bold"
    )

  # Add gene labels if available
  if (exists("top_genes") && nrow(top_genes) > 0) {
    volcano_plot <- volcano_plot +
      ggrepel::geom_text_repel(
        data = top_genes,
        ggplot2::aes(label = gene),
        size = 3,
        max.overlaps = Inf,
        min.segment.length = 0.2,
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50",
        segment.size = 0.3,
        segment.alpha = 0.6,
        force = 2,
        force_pull = 2,
        max.time = 2,
        max.iter = 10000
      )
  }

  return(volcano_plot)
}
