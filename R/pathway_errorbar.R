utils::globalVariables(c("group", "name", "value", "feature", "negative_log10_p", "group_nonsense", "nonsense", "pathway_class", "p_adjust", "log_2_fold_change", "transform_sample_counts", "column_to_rownames", "txtProgressBar", "setTxtProgressBar", "utils", "p_label"))
#' Relative-abundance summaries for differential functional features
#'
#' Plot group mean relative abundances with standard-deviation error bars, optional two-group mean differences and adjusted p-value labels.
#'
#' @param order A parameter controlling the ordering of the rows in the figure. The options are: "p_values" (order by p-values), "name" (order by pathway name), "group" (order by the group with the highest mean relative abundance), or "pathway_class" (order by the pathway category).
#'
#' @param colors A vector of colors to be used to represent the groups in the figure. Each color corresponds to a group.
#'
#' @importFrom stats sd
#'
#' @export
#'
#' @param abundance Numeric feature-by-sample abundance matrix or data frame with feature row names. Sample totals must be positive; values are converted to relative abundances by column.
#'
#' @param daa_results_df Annotated differential-result data frame containing `feature`, `method`, `p_adjust` and the annotation column selected by `x_lab`. Select a single method and contrast first, so each selected feature occurs once.
#'
#' @param Group Group vector in exactly the same order as the abundance columns. Use one non-missing group assignment per sample.
#'
#' @param ko_to_kegg Logical. If `TRUE`, use KEGG pathway labels and class annotations; requires `pathway_class`.
#'
#' @param p_values_threshold Retain features with finite adjusted p-values strictly below this cutoff.
#'
#' @param select Optional vector of feature IDs to retain after significance filtering. If `NULL`, use all significant features. More than 30 retained features raises an error; supply `select` to reduce the set.
#'
#' @param p_value_bar Logical. Show the two-group mean-difference panel when `TRUE`. Automatically disabled unless exactly two groups are observed.
#'
#' @param x_lab Name of the annotation column in `daa_results_df`, not an axis title. Defaults to `"pathway_name"` for KEGG or `"description"` otherwise.
#'
#' @return A patchwork plot combining abundance summaries, adjusted p-value labels and optional KEGG class and two-group panels. The `microfun_min_height` attribute suggests a minimum figure height in inches.
#'
#' @details Features are matched to abundance rows by ID. Missing annotations are dropped. Adjusted p-values are displayed in scientific notation. The two-group panel summarizes differences in observed mean relative abundance; it does not display a backend-specific model coefficient.
#'
#' @seealso [pathway_daa()], [pathway_annotation()]
#'
#' @md
pathway_errorbar <-
  function(abundance,
           daa_results_df,
           Group,
           ko_to_kegg = FALSE,
           p_values_threshold = 0.05,
           order = "group",
           select = NULL,
           p_value_bar = TRUE,
           colors = NULL,
           x_lab = NULL) {
    # Validate inputs at the start of the function.
    if(!is.matrix(abundance) && !is.data.frame(abundance)) {
      stop("'abundance' must be a matrix or data frame")
    }

    if(!is.data.frame(daa_results_df)) {
      stop("'daa_results_df' must be a data frame")
    }

    # Check required columns.
    required_cols <- c("feature", "method", "group1", "group2", "p_adjust")
    missing_cols <- setdiff(required_cols, colnames(daa_results_df))
    if(length(missing_cols) > 0) {
      stop("Missing required columns in daa_results_df: ",
           paste(missing_cols, collapse = ", "))
    }

    # Validate input data.
    if (length(Group) != ncol(abundance)) {
      stop("Length of Group must match number of columns in abundance matrix")
    }

    # Check the number of significant features.
    sig_features <- sum(daa_results_df$p_adjust < 0.05)

    # Identify pathways with missing annotation
    label_column <- if (is.null(x_lab)) {
      if (ko_to_kegg) "pathway_name" else "description"
    } else x_lab
    missing_pathways <- daa_results_df$feature[is.na(daa_results_df[[label_column]])]

    # Inform the user about the missing annotations
    if(length(missing_pathways) > 0) {
      message("The following pathways are missing annotations and have been excluded: ",
              paste(missing_pathways, collapse = ", "))
      message("You can use the 'pathway_annotation' function to add annotations for these pathways.")
    }

    # Get the names of all columns in the data frame
    column_names <- colnames(daa_results_df)

    # Check if there are more than two 'group' columns
    group_columns <- grepl("^group", column_names)
    if (sum(group_columns) > 2) {
      p_value_bar <- FALSE
      message(
        "There are more than two 'group' columns in the 'daa_results_df' data frame. As a result, it is not possible to compute the log2 fold values. The 'p_value_bar' has been automatically set to FALSE."
      )
    }

    # Exclude rows with missing pathway annotation
    x_lab <- label_column
    if (!x_lab %in% names(daa_results_df))
      stop("The x_lab column '", x_lab, "' is missing. Provide it or annotate the DAA results.")
    if (ko_to_kegg && !"pathway_class" %in% names(daa_results_df))
      stop("KEGG results require a 'pathway_class' column.")
    daa_results_df <- daa_results_df[!is.na(daa_results_df[[x_lab]]), , drop = FALSE]

    if (nlevels(factor(daa_results_df$method)) != 1) {
      stop(
        "The 'method' column in the 'daa_results_df' data frame contains more than one method. Please filter it to contain only one method."
      )
    }

    if (nlevels(factor(daa_results_df$group1)) != 1 || nlevels(factor(daa_results_df$group2)) != 1) {
      stop(
        "The 'group1' or 'group2' column in the 'daa_results_df' data frame contains more than one group. Please filter each to contain only one group."
      )
    }

    if (is.null(colors))
      colors <- c("#B4CBE0", "#F5BB7C", "#6faa3e", "#e8a825", "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")
    Group <- droplevels(as.factor(Group))
    if (is.null(names(colors)))
      colors <- .microfun_extend_colours(colors, nlevels(Group))
    if (nlevels(Group) != 2L)
      p_value_bar <- FALSE

    errorbar_abundance_mat <- as.matrix(abundance)

    daa_results_filtered_df <-
      daa_results_df[is.finite(daa_results_df$p_adjust) & daa_results_df$p_adjust < p_values_threshold, , drop = FALSE]

    if (!is.null(select)) {
      daa_results_filtered_sub_df <-
        daa_results_filtered_df[daa_results_filtered_df$feature %in% select, ]
    } else {
      daa_results_filtered_sub_df <- daa_results_filtered_df
    }

    if (nrow(daa_results_filtered_sub_df) > 30) {
      message(
        paste0(
          "The number of features with statistical significance exceeds 30, leading to suboptimal visualization. ",
          "Please use 'select' to reduce the number of features.\n",
          "Currently, you have these features: ",
          paste(paste0('"', daa_results_filtered_sub_df$feature, '"'), collapse = ", "), ".\n",
          "You can find the statistically significant features with the following command:\n",
          "daa_results_df %>% filter(p_adjust < 0.05) %>% select(c(\"feature\",\"p_adjust\"))"
        )
      )
      stop()
    }

    if (nrow(daa_results_filtered_sub_df) == 0){
      stop(
        "Visualization with 'pathway_errorbar' cannot be performed because there are no features with statistical significance. ",
        "For possible solutions, please check the FAQ section of the tutorial."
      )
    }
    # Convert to relative abundance
    relative_abundance_mat <- sweep(errorbar_abundance_mat, 2L, colSums(errorbar_abundance_mat), "/")

    # Subset to only include the features present in daa_results_filtered_sub_df$feature
    daa_results_filtered_sub_df$feature <- as.character(daa_results_filtered_sub_df$feature)
    if (!all(daa_results_filtered_sub_df$feature %in% rownames(relative_abundance_mat)) ||
        anyDuplicated(daa_results_filtered_sub_df$feature))
      stop("Each selected feature must occur once in the DAA results and match an abundance row.")
    sub_relative_abundance_mat <- relative_abundance_mat[daa_results_filtered_sub_df$feature, , drop = FALSE]

    # Create a matrix for the error bars
    error_bar_matrix <- cbind(
      sample = colnames(sub_relative_abundance_mat),
      group = Group,
      t(sub_relative_abundance_mat)
    )
    error_bar_df <- as.data.frame(error_bar_matrix)

    error_bar_df$group <- factor(Group,levels = levels(as.factor(Group)))

      error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df,-c(sample, group))

    error_bar_pivot_longer_tibble <-
      mutate(error_bar_pivot_longer_df, group = as.factor(group))

    error_bar_pivot_longer_tibble$sample <-
      factor(error_bar_pivot_longer_tibble$sample)

    error_bar_pivot_longer_tibble$name <-
      factor(error_bar_pivot_longer_tibble$name)

    error_bar_pivot_longer_tibble$value <-
      as.numeric(error_bar_pivot_longer_tibble$value)

    error_bar_pivot_longer_tibble_summarised <-
      error_bar_pivot_longer_tibble %>% group_by(name, group) %>%
      summarise(mean = mean(value), sd = stats::sd(value), .groups = "drop")

    switch(
      order,
      "p_values" = {
        order <- order(daa_results_filtered_sub_df$p_adjust)
      },
      "name" = {
        order <- order(daa_results_filtered_sub_df$feature)
      },
      "group" = {
        daa_results_filtered_sub_df$pro <- 1
        for (i in levels(error_bar_pivot_longer_tibble_summarised$name)) {
          error_bar_pivot_longer_tibble_summarised_sub <-
            error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name ==
                                                       i,]
          pro_group <- as.character(error_bar_pivot_longer_tibble_summarised_sub$group[
            order(error_bar_pivot_longer_tibble_summarised_sub$mean, decreasing = TRUE, na.last = TRUE)[1L]])
          daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature ==
                                        i,]$pro <- pro_group
        }
        order <-
          order(daa_results_filtered_sub_df$pro,
                daa_results_filtered_sub_df$p_adjust)
      },
      "pathway_class" = {
        if (!"pathway_class" %in% colnames(daa_results_filtered_sub_df)) {
          stop(
            "The 'pathway_class' column is missing in the 'daa_results_filtered_sub_df' data frame. ",
            "Please use the 'pathway_annotation' function to annotate the 'pathway_daa' results."
          )
        }
        order <- order(
          daa_results_filtered_sub_df$pathway_class,
          daa_results_filtered_sub_df$p_adjust
        )
      },
      {
        order <- order
      }
    )

    daa_results_filtered_sub_df <-
      daa_results_filtered_sub_df[order,]

    error_bar_pivot_longer_tibble_summarised_ordered <-
      data.frame(
        name = NULL,
        group = NULL,
        mean = NULL,
        sd = NULL
      )

    for (i in daa_results_filtered_sub_df$feature) {
      error_bar_pivot_longer_tibble_summarised_ordered <-
        rbind(
          error_bar_pivot_longer_tibble_summarised_ordered,
          error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name ==
                                                     i,]
        )
    }

    if (ko_to_kegg == FALSE){
      error_bar_pivot_longer_tibble_summarised_ordered[, x_lab] <-
        rep(daa_results_filtered_sub_df[[x_lab]], each = length(levels(
          factor(error_bar_pivot_longer_tibble_summarised_ordered$group)
        )))
    }

    if (ko_to_kegg == TRUE) {
      error_bar_pivot_longer_tibble_summarised_ordered$pathway_class <-
        rep(daa_results_filtered_sub_df$pathway_class,
            each = length(levels(
              factor(error_bar_pivot_longer_tibble_summarised_ordered$group)
            )))
    }

    error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, levels = rev(daa_results_filtered_sub_df$feature))

    bar_errorbar <-
      ggplot2::ggplot(error_bar_pivot_longer_tibble_summarised_ordered,
             ggplot2::aes(mean, name, fill = group)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(xmax = mean + sd, xmin = 0),
        position = ggplot2::position_dodge(width = 0.8),
        width = 0.5,
        linewidth = 0.5,
        color = "black"
      ) +
      ggplot2::geom_bar(stat = "identity",
               position = ggplot2::position_dodge(width = 0.8),
               width = 0.8) +
      GGally::geom_stripped_cols(width = 10) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::scale_color_manual(values = colors) +
      ggprism::theme_prism(base_size = 12) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(labels = rev(daa_results_filtered_sub_df[[x_lab]])) +
      ggplot2::labs(x = "Relative Abundance", y = NULL) +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(linewidth = 0.5),
        axis.ticks.x = ggplot2::element_line(linewidth = 0.5),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = 10, color = "black"),
        axis.text.x = ggplot2::element_text(margin = ggplot2::margin(r = 0)),
        axis.text.y = ggplot2::element_text(
          size = 10,
          color = "black",
          margin = ggplot2::margin(b = 6)
        ),
        axis.title.x =  ggplot2::element_text(
          size = 10,
          color = "black",
          hjust = 0.5
        ),
        legend.position = "top",
        legend.key.size = ggplot2::unit(0.1, "cm"),
        legend.direction = "vertical",
        legend.justification = "left",
        legend.text = ggplot2::element_text(size = 8, face = "bold"),
        legend.box.just = "right",
        plot.margin = ggplot2::margin(0, 0.5, 0.5, 0, unit = "cm")
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.margin = ggplot2::unit(c(1, 1, 1, 1), "cm")
      )

    if (ko_to_kegg == TRUE) {
      pathway_class_group_mat <-
        daa_results_filtered_sub_df$pathway_class %>%
        table() %>%
        data.frame() %>% tibble::column_to_rownames(".")
      pathway_class_group <- data.frame(.= unique(daa_results_filtered_sub_df$pathway_class),Freq = pathway_class_group_mat[unique(daa_results_filtered_sub_df$pathway_class),])
      start <-
        c(1, head(rev(pathway_class_group$Freq), -1L)) %>%
        cumsum()
      end <- cumsum(rev(pathway_class_group$Freq))
      ymin <- start - 1 / 2
      ymax <- end + 1 / 2
      nPoints <- length(start)
      pCol <- .microfun_extend_colours(c("#D51F26","#272E6A","#208A42","#89288F","#F47D2B",
                 "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                 "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                 "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D"), nPoints, "function")
      pFill <- pCol
      for (i in 1:nPoints)  {
        bar_errorbar <- bar_errorbar +
          ggplot2::annotation_custom(
            grob = grid::rectGrob(
              gp = grid::gpar(
                col = pCol[i],
                fill = pFill[i],
                lty = NULL,
                lwd = NULL,
                alpha = 0.2
              )
            ),
            xmin = ggplot2::unit(-2, 'native'),
            xmax = ggplot2::unit(0, 'native'),
            ymin = ggplot2::unit(ymin[i], 'native'),
            ymax = ggplot2::unit(ymax[i], 'native')
          )
      }
    }
    daa_results_filtered_sub_df <-
      cbind(
        daa_results_filtered_sub_df,
        negative_log10_p = -log10(daa_results_filtered_sub_df$p_adjust),
        group_nonsense = "nonsense",
        log_2_fold_change = NA
      )

    for (i in daa_results_filtered_sub_df$feature){
      mean <- error_bar_pivot_longer_tibble_summarised_ordered[error_bar_pivot_longer_tibble_summarised_ordered$name %in% i,]$mean
      daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature==i,]$log_2_fold_change <- log2(mean[1]/mean[2])
    }
    daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature,levels = rev(daa_results_filtered_sub_df$feature))
    p_values_bar <- daa_results_filtered_sub_df %>%
      ggplot2::ggplot(ggplot2::aes(feature, log_2_fold_change, fill = log_2_fold_change > 0)) +
      ggplot2::geom_bar(stat = "identity",
               position = ggplot2::position_dodge(width = 0.8),
               width = 0.8) +
      ggplot2::labs(y = "log2 fold change", x = NULL) +
      GGally::geom_stripped_cols() +
      ggplot2::scale_fill_manual(values = c("FALSE" = "#B4CBE0", "TRUE" = "#F5BB7C"), guide = "none") +
      ggplot2::geom_hline(ggplot2::aes(yintercept = 0),
                 linetype = 'dashed',
                 color = 'black') +
      ggprism::theme_prism(base_size = 12) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::theme(
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.line.x = ggplot2::element_line(linewidth = 0.5),
        axis.ticks.x = ggplot2::element_line(linewidth = 0.5),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        axis.text = ggplot2::element_text(size = 10, color = "black"),
        axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          size = 10,
          color = "black",
          margin = ggplot2::margin(b = 6)
        ),
        axis.title.x =  ggplot2::element_text(
          size = 11,
          color = "black",
          hjust = 0.5
        ),
        legend.position = "none"
      ) +
      ggplot2::coord_flip()

    if (ko_to_kegg == TRUE) {
      pathway_class_y <- (ymax + ymin) / 2 - 0.5
      pathway_class_plot_df <-
        as.data.frame(
          cbind(
            nonsense = "nonsense",
            pathway_class_y = pathway_class_y,
            pathway_class = rev(unique(
              daa_results_filtered_sub_df$pathway_class
            ))
          )
        )
      pathway_class_plot_df$pathway_class_y <-
        as.numeric(pathway_class_plot_df$pathway_class_y)
      pathway_class_annotation <-
        pathway_class_plot_df %>% ggplot2::ggplot(ggplot2::aes(nonsense, pathway_class_y)) + ggplot2::geom_text(
          ggplot2::aes(nonsense, pathway_class_y, label = pathway_class),
          size = 3.5,
          color = "black",
          fontface = "bold",
          family = "sans"
        ) +
        ggplot2::scale_y_discrete(position = "right") +
        ggprism::theme_prism(base_size = 12) +
        ggplot2::theme(
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.major.x = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_blank(),
          plot.margin = ggplot2::unit(c(0, 0.2, 0, 0), "cm"),
          axis.title.y =  ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          legend.position = "none"
        )
    }
    daa_results_filtered_sub_df$p_label <-
      formatC(daa_results_filtered_sub_df$p_adjust, format = "e", digits = 3L)
    p_label_size <- 3.5
    p_label_width <- grid::grobWidth(grid::textGrob(
      daa_results_filtered_sub_df$p_label,
      gp = grid::gpar(fontfamily = "sans", fontface = "bold", fontsize = p_label_size * ggplot2::.pt)
    )) + grid::unit(2, "mm")
    panel_widths <- function(widths) {
      grid::unit.c(grid::unit(widths, "null"), p_label_width)
    }
    p_annotation <- daa_results_filtered_sub_df %>%
      ggplot2::ggplot(ggplot2::aes(group_nonsense, feature)) +
      ggplot2::geom_text(
        ggplot2::aes(label = p_label),
        size = p_label_size,
        color = "black",
        fontface = "bold",
        family = "sans"
      ) +
      ggplot2::labs(y = "p-value (adjusted)") +
      ggplot2::scale_y_discrete(position = "right",
        limits = levels(daa_results_filtered_sub_df$feature), drop = FALSE) +
      ggprism::theme_prism(base_size = 12) +
      ggplot2::theme(
        axis.ticks = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0, 0.2, 0, 0), "cm"),
        axis.title.y =  ggplot2::element_text(
          size = 11,
          color = "black",
          vjust = 0
        ),
        axis.title.x = ggplot2::element_blank(),
        legend.position = "none"
      )
    if (p_value_bar == TRUE) {
      if (ko_to_kegg == TRUE) {
        combination_bar_plot <-
          patchwork::wrap_plots(pathway_class_annotation, bar_errorbar, p_values_bar, p_annotation,
            ncol = 4, widths = panel_widths(c(1, 1.2, 0.5)))
      }
      else{
        combination_bar_plot <-
          patchwork::wrap_plots(bar_errorbar, p_values_bar, p_annotation,
            ncol = 3, widths = panel_widths(c(2.3, 0.7)))
      }
    }else{
      if (ko_to_kegg == TRUE) {
        combination_bar_plot <-
          patchwork::wrap_plots(pathway_class_annotation, bar_errorbar, p_annotation,
            ncol = 3, widths = panel_widths(c(1, 1.2)))
      }
      else{
        combination_bar_plot <-
          patchwork::wrap_plots(bar_errorbar, p_annotation, ncol = 2, widths = panel_widths(2.5))
      }
    }
    attr(combination_bar_plot, "microfun_min_height") <- (max(24, 4.3 * nrow(daa_results_filtered_sub_df)) + 50) / 25.4
    return(combination_bar_plot)
  }
