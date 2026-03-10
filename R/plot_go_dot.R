#' Perform GO enrichment analysis and return dotplot
#'
#' @param deg_table DEG results data frame
#' @param logfc_col Column name for log2 fold change
#' @param padj_col Column name for adjusted p-value
#' @param logfc_threshold Log2FC cutoff
#' @param padj_threshold Adjusted p-value cutoff
#' @param orgdb OrgDb annotation package
#' @param ont Ontology: "BP", "MF", or "CC"
#' @param show_category Number of GO categories to display
#'
#' @return ggplot object
#' @export
plot_go_dot <- function(
    deg_table,
    logfc_col = "avg_log2FC",
    padj_col = "p_val_adj",
    logfc_threshold = 0.5,
    padj_threshold = 0.05,
    orgdb = org.Mm.eg.db::org.Mm.eg.db,
    ont = "ALL",
    show_category = 15
) {

  required_pkgs <- c("clusterProfiler", "enrichplot", "AnnotationDbi")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop("Please install: ", paste(missing_pkgs, collapse = ", "))
  }

  if (!all(c(logfc_col, padj_col) %in% colnames(deg_table))) {
    stop("Specified logfc_col or padj_col not found.")
  }

  # ---- Prepare data ----
  df <- deg_table
  df$log2FoldChange <- df[[logfc_col]]
  df$padj <- df[[padj_col]]
  df$SYMBOL <- rownames(df)

  df_sig <- df[
    !is.na(df$padj) &
      df$padj < padj_threshold &
      abs(df$log2FoldChange) > logfc_threshold,
  ]

  if (nrow(df_sig) < 5) {
    stop("Too few significant genes for enrichment.")
  }

  # ---- Map to ENTREZ ----
  gene.df <- clusterProfiler::bitr(
    df_sig$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = orgdb
  )

  if (nrow(gene.df) == 0) {
    stop("No SYMBOL mapped to ENTREZ.")
  }

  df_sig <- merge(df_sig, gene.df, by = "SYMBOL")

  # ---- Enrichment ----
  ego <- clusterProfiler::enrichGO(
    gene = df_sig$ENTREZID,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = ont,
    pvalueCutoff = padj_threshold,
    qvalueCutoff = padj_threshold
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    stop("No enriched GO terms found.")
  }

  # ---- Remove redundancy ----
  ego <- clusterProfiler::simplify(
    ego,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )

  # ---- Dotplot ----
  p <- enrichplot::dotplot(
    ego,
    showCategory = show_category,
    orderBy = "GeneRatio"
  ) +
    ggplot2::labs(
      title = paste0("GO ", ont, " Enrichment"),
      subtitle = paste0(
        "|log2FC| > ", logfc_threshold,
        " & FDR < ", padj_threshold
      )
    ) +
    ggplot2::theme_bw()

  return(p)
}
