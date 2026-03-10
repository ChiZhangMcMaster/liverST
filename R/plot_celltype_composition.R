#' Plot stacked percent cell type composition (row-normalized)
#'
#' Each cell type is re-normalized so group percentages sum to 100%.
#' Cell types are ordered top-to-bottom by largest difference across groups.
#'
#' @param df_pct Data frame returned from prepare_celltype_percent()
#' @param celltype_col Column name for cell type (default = "celltype")
#' @param percent_col Column name for percent values (default = "Percent")
#' @param group_col Column name for group (default = "group")
#' @param reference_group Optional reference group (default = NULL)
#' @param reference_line Numeric value for vertical reference line (default = 50)
#'
#' @return A ggplot object
plot_celltype_composition <- function(
    df_pct,
    celltype_col = "celltype",
    group_col = "group",
    experimental_group = NULL,
    reference_line = 50
) {

  df_plot <- df_pct

  # ----------------------------
  # Determine experimental group
  # ----------------------------
  if (is.null(experimental_group)) {
    experimental_group <- unique(df_plot[[group_col]])[1]
  }

  # Experimental group first → left in bars, top in legend
  df_plot[[group_col]] <- factor(
    df_plot[[group_col]],
    levels = c(
      experimental_group,
      setdiff(unique(df_plot[[group_col]]), experimental_group)
    )
  )

  # ----------------------------
  # Order cell types by experimental − other difference (using NormPercent)
  # ----------------------------
  exp_df <- df_plot[df_plot[[group_col]] == experimental_group,
                    c(celltype_col, "NormPercent")]
  names(exp_df)[2] <- "exp_percent"

  df_plot <- merge(df_plot, exp_df, by = celltype_col)
  df_plot$diff <- df_plot$exp_percent - df_plot$NormPercent

  # Compute mean difference per cell type (experimental − other)
  diff_order <- aggregate(diff ~ get(celltype_col), df_plot, mean)
  names(diff_order)[1] <- celltype_col

  # Sort from largest positive to largest negative (pure value)
  diff_order <- diff_order[order(-diff_order$diff), ]

  # Factor levels: top = largest positive diff, bottom = largest negative diff
  df_plot[[celltype_col]] <- factor(df_plot[[celltype_col]],
                                    levels = rev(diff_order[[celltype_col]]))

  # ----------------------------
  # Reverse default ggplot2 colors
  # ----------------------------
  n_groups <- length(levels(df_plot[[group_col]]))
  default_colors <- scales::hue_pal()(n_groups)  # ggplot2 default palette
  reversed_colors <- rev(default_colors)

  # ----------------------------
  # Plot
  # ----------------------------
  ggplot2::ggplot(
    df_plot,
    ggplot2::aes(
      x = NormPercent,
      y = .data[[celltype_col]],
      fill = .data[[group_col]]
    )
  ) +
    ggplot2::geom_col(
      width = 0.75,
      position = ggplot2::position_stack(reverse = TRUE)
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", NormPercent)),
      position = ggplot2::position_stack(vjust = 0.5, reverse = TRUE),
      size = 3.5
    ) +
    ggplot2::geom_vline(
      xintercept = reference_line,
      linetype = "dashed",
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = reversed_colors) +
    ggplot2::labs(
      x = "Relative distribution across groups (row-normalized)",
      y = "Cell Type"
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank()
    )
}
