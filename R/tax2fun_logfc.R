#' Effect-size bar plots for taxon-oriented results
#'
#' Create one bar plot per taxon, ranking functions by absolute `logFC`. Positive and negative effects use the Up and Down colours. Recorded effect types label the effect axis.
#'
#' @param pathway_sdaa_result A **named list** where each element is a
#'   `data.frame` of differential results for one taxon (or function), containing
#'   at least the columns `feature`, `logFC`, `p_adjust`, `group1`, `group2`,
#'   and optionally `pathway_name`.
#'
#' @param wrap_width Integer. Character width used to line-wrap labels (via
#'   `stringr::str_wrap`). Default `45`.
#'
#' @return A named `list` of `ggplot` objects, one per entry of
#'   `pathway_sdaa_result`. A plot is returned for each named entry; after filtering, some plots may contain no bars.
#'
#' @export
#'
#' @param top_n Positive integer giving the maximum number of rows per plot, ranked by absolute effect. Significant rows do not bypass this final limit.
#'
#' @param padj_cutoff Adjusted p-value cutoff used in the internal row-selection step. The final plot is still capped at `top_n` and is not restricted to significant rows.
#'
#' @details The column `logFC` may contain a model coefficient rather than a log fold change; its meaning follows `effect_type`. Input must contain `feature`, `logFC`, `p_adjust`, `group1` and `group2`. Restrict input to one contrast before plotting. Missing p-values are treated as 1 and missing effects are omitted.
#'
#' @seealso [fun2tax_logfc()], [pathway_sdaa()]
#'
#' @md
tax2fun_logfc <- function(pathway_sdaa_result,
                               top_n = 15,
                               padj_cutoff = 0.05,
                               wrap_width = 45) {
  taxa = names(pathway_sdaa_result)
  label_col = "pathway_name"
  out <- lapply(taxa, function(tn) {
    if (!tn %in% names(pathway_sdaa_result)) return(NULL)
    plot_pathway_logfc_one(
      df = pathway_sdaa_result[[tn]],
      title = tn,
      top_n = top_n,
      padj_cutoff = padj_cutoff,
      label_col = label_col,
      wrap_width = wrap_width
    )
  })
  names(out) <- taxa
  Filter(Negate(is.null), out)
}

plot_pathway_logfc_one <- function(df,
                                   title = NULL,
                                   padj_cutoff = 0.05,
                                   top_n = 15,
                                   label_col = c("pathway_name", "feature"),
                                   wrap_width = 45) {
  label_col <- label_col[1]

  if (!label_col %in% names(df)) label_col <- "feature"
  df2 <- df %>%
    mutate(
      p_adjust = ifelse(is.na(p_adjust), 1, p_adjust),
      label = .data[[label_col]],
      label = ifelse(is.na(label) | label == "", .data[["feature"]], label),
      label = stringr::str_wrap(label, width = wrap_width)
    ) %>%
    arrange(desc(abs(logFC))) %>%
    # keep the strongest effects, but always keep significant ones
    filter(row_number() <= top_n | p_adjust <= padj_cutoff) %>%
    slice_head(n = top_n) %>%
    filter(!is.na(logFC)) %>%
    mutate(label = stats::reorder(label, logFC))

  effect_types <- if ("effect_type" %in% names(df)) {
    unique(df$effect_type[!is.na(df$effect_type) & nzchar(df$effect_type)])
  } else if ("method" %in% names(df) && identical(unique(as.character(df$method)), "Maaslin2")) {
    "Maaslin2 AST coefficient"
  }
  effect_label <- if (length(effect_types) == 1L) effect_types else "logFC"
  g1 <- ggplot(df2, aes(x = label, y = logFC, fill = logFC > 0)) +
    geom_col(width = 0.8) +
    coord_flip() +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    scale_fill_manual(values = c("FALSE" = "#B4CBE0", "TRUE" = "#F5BB7C"), guide = "none") +
    labs(
      x = NULL, y = effect_label,
      title = if (is.null(title)) "" else title,
      subtitle = paste0(unique(df$group1), " vs ", unique(df$group2))
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title.position = "plot",
      plot.margin = margin(10, 18, 10, 36),   # add left margin to avoid clipping
      axis.text.y = element_text(hjust = 1)
    )
  g1
}


# ---- Merged from fun2tax_complete.R ----

#' Plot signed taxon effects for each function
#'
#' Create one effect-size bar plot per function from function-oriented results. Taxon IDs in `feature` provide bar labels; no statistical analysis is rerun.
#'
#' @return A named list of ggplot objects.
#'
#' @export
#'
#' @param pathway_sdaa_result Named function-oriented list containing `feature`, `logFC`, `p_adjust`, `group1` and `group2` in each result table.
#'
#' @param ... Arguments passed to `tax2fun_logfc()`: `top_n`, `padj_cutoff` and `wrap_width`.
#'
#' @seealso [tax2fun_logfc()]
#'
#' @md
fun2tax_logfc <- function(pathway_sdaa_result, ...) {
  .microfun_result_direction(pathway_sdaa_result, "function")
  result <- lapply(pathway_sdaa_result, function(z) {
    z$pathway_name <- as.character(z$feature)
    z
  })
  tax2fun_logfc(result, ...)
}
