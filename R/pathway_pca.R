#' PCA of functional abundance profiles with marginal densities
#'
#' Compute PCA on the transposed abundance matrix with centering and feature scaling. Plot the first two components, 95 percent group ellipses and marginal density plots.
#'
#' @param group A character string specifying the column name in metadata that contains
#'        group information for samples (e.g., "treatment", "condition", "group").
#'
#' @importFrom magrittr %>%
#'
#' @importFrom dplyr select
#'
#' @importFrom stats prcomp
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual stat_ellipse labs theme_classic theme element_line element_text element_blank geom_vline geom_hline geom_density scale_fill_manual coord_flip
#'
#' @importFrom aplot insert_top insert_right
#'
#' @importFrom ggplotify as.ggplot
#'
#' @export
#'
#' @param abundance Numeric matrix or data frame with at least two features in rows and three samples in columns. Supply unique sample column names, finite values, no missing values and no constant feature rows.
#'
#' @param metadata Data frame containing a unique `sample_name` column and the grouping column. Its sample IDs must match the abundance columns exactly; row order may differ.
#'
#' @param colors Optional character vector with one valid colour per observed group. Defaults to the built-in group colours.
#'
#' @return A ggplot object created by converting the composite PCA and marginal-density layout with `ggplotify::as.ggplot()`.
#'
#' @details No relative-abundance or logarithmic transformation is applied before PCA. Group ellipses and density estimates need enough observations and variation within each group; small or degenerate groups may produce plotting warnings.
#'
#' @seealso [pathway_heatmap()]
#'
#' @md
pathway_pca <- function(abundance,
                        metadata,
                        group,
                        colors = NULL) {
  # Input validation
  # Check if inputs are missing
  if (missing(abundance)) {
    stop("Abundance matrix is required")
  }
  if (missing(metadata)) {
    stop("Metadata is required")
  }
  if (missing(group)) {
    stop("Group variable name is required")
  }

  # Check abundance matrix
  if (!is.matrix(abundance) && !is.data.frame(abundance)) {
    stop("Abundance must be a matrix or data frame")
  }
  if (any(is.na(abundance))) {
    stop("Abundance matrix contains missing values (NA)")
  }
  if (!all(apply(abundance, 2, is.numeric))) {
    stop("Abundance matrix must contain only numeric values")
  }
  if (nrow(abundance) < 2) {
    stop("Abundance matrix must contain at least 2 pathways")
  }
  if (ncol(abundance) < 3) {
    stop("Abundance matrix must contain at least 3 samples")
  }

  # Check metadata
  if (!is.data.frame(metadata)) {
    stop("Metadata must be a data frame")
  }
  if (!"sample_name" %in% colnames(metadata)) {
    stop("Metadata must contain a 'sample_name' column")
  }
  if (!group %in% colnames(metadata)) {
    stop(sprintf("Group column '%s' not found in metadata", group))
  }

  # Check sample names match
  if (anyDuplicated(metadata$sample_name) || anyNA(metadata$sample_name) ||
      length(metadata$sample_name) != ncol(abundance)) {
    stop("Number of unique samples in metadata does not match abundance matrix")
  }
  if (is.null(colnames(abundance)) || anyDuplicated(colnames(abundance)) ||
      !all(colnames(abundance) %in% metadata$sample_name)) {
    stop("Some sample names in abundance matrix are not found in metadata")
  }

  # Check group variable
  group_values <- metadata[[group]]
  if (length(unique(group_values)) < 1) {
    stop("Group variable must have at least one level")
  }
  if (!is.factor(group_values)) {
    warning("Converting group variable to factor")
    metadata[[group]] <- factor(group_values)
  }

  # Check colors if provided
  if (!is.null(colors)) {
    if (!is.vector(colors) || !is.character(colors)) {
      stop("Colors must be provided as a character vector")
    }
    if (!all(sapply(colors, function(x) tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)))) {
      stop("Invalid color names provided")
    }
    levels <- length(levels(factor(metadata[[group]])))
    if (length(colors) != levels) {
      stop(sprintf("Number of colors (%d) does not match number of groups (%d)",
                  length(colors), levels))
    }
  }

  # due to NSE notes in R CMD check
  PC1 = PC2 = Group = NULL

  # Perform PCA on the abundance data, keeping the first two principal components
  pca_fit <- stats::prcomp(t(abundance), center = TRUE, scale = TRUE)
  pca_axis <- pca_fit$x[, 1:2, drop = FALSE]

  # Calculate the proportion of total variation explained by each PC
  pca_proportion <- pca_fit$sdev[1:2]^2 / sum(pca_fit$sdev^2) * 100

  # Combine the PCA results with the metadata information
  pca <- as.data.frame(pca_axis)
  pca$Group <- metadata[[group]][match(rownames(pca), metadata$sample_name)]

  levels <- length(levels(factor(pca$Group)))

  # Set default colors if colors are not provided
  if (is.null(colors)) {
    colors <- c("#B4CBE0", "#F5BB7C","#208A42","#89288F","#F47D2B",
                "#FEE500","#8A9FD1","#C06CAB","#E6C2DC","#90D5E4",
                "#89C75F","#F37B7D","#9983BD","#D24B27","#3BBCA8",
                "#6E4B9E","#0C727C", "#7E1416","#D8A767","#3D3D3D")[1:levels]
  }

  # Ensure the number of colors matches the number of levels in Group
  if (length(colors) != levels) {
    stop("The length of colors vector must match the number of levels in Group")
  }

  # Create a ggplot object for the PCA scatter plot
  Fig1a.taxa.pca <- ggplot2::ggplot(pca,ggplot2::aes(PC1,PC2))+
    ggplot2::geom_point(size=4,ggplot2::aes(color=Group),show.legend = T)+
    ggplot2::scale_color_manual(values=colors)+
    ggplot2::stat_ellipse(ggplot2::aes(color = Group),fill="white",geom = "polygon",
                 level=0.95,alpha = 0.01,show.legend = F)+
    ggplot2::labs(x=paste0("PC1(",round(pca_proportion[1],1),"%)"),y=paste0("PC2(",round(pca_proportion[2],1),"%)"),color = group)+
    ggplot2::theme_classic()+
    ggplot2::theme(axis.line=ggplot2::element_line(colour = "black"),
          axis.title=ggplot2::element_text(color="black",face = "bold"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(color="black",size=10,face = "bold"),
          legend.text = ggplot2::element_text(size = 16, face = "bold"),
          legend.title = ggplot2::element_text(size = 16, face = "bold"))+
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black")

  # Create a ggplot object for the density plot of PC1
  # Plot the density of PC1 colored by Group
  Fig1a.taxa.pc1.density <-
    ggplot2::ggplot(pca) +
    ggplot2::geom_density(ggplot2::aes(x=PC1, group=Group, fill=Group), # Plot the density of PC1
                 color="black", alpha=1, position = 'identity',show.legend = F) +
    ggplot2::scale_fill_manual(values=colors) + # Manually set the fill color for each group
    ggplot2::theme_classic() + # Set the classic theme
    ggplot2::scale_y_discrete(expand = c(0,0.001)) + # Scale the y-axis with a small expansion to improve appearance
    ggplot2::labs(x=NULL, y=NULL) + # Remove x and y labels
    ggplot2::theme(axis.text.x=ggplot2::element_blank(), # Remove x axis text
          axis.ticks.x = ggplot2::element_blank()) # Remove x axis ticks

  # Plot the density of PC2 colored by Group
  Fig1a.taxa.pc2.density <-
    ggplot2::ggplot(pca) +
    ggplot2::geom_density(ggplot2::aes(x=PC2, group=Group, fill=Group), # Plot the density of PC2
                 color="black", alpha=1, position = 'identity',show.legend = F) +
    ggplot2::scale_fill_manual(values=colors) + # Manually set the fill color for each group
    ggplot2::theme_classic() + # Set the classic theme
    ggplot2::scale_y_discrete(expand = c(0,0.001)) + # Scale the y-axis with a small expansion to improve appearance
    ggplot2::labs(x=NULL, y=NULL) + # Remove x and y labels
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), # Remove y axis text
          axis.ticks.y = ggplot2::element_blank()) + # Remove y axis ticks
    ggplot2::coord_flip() # Flip the plot to change it from a horizontal to a vertical orientation

  # Combine the two plots into a single plot with the PC1 plot on top and the PC2 plot on the right
  Fig1a.taxa.pca %>%
    aplot::insert_top(Fig1a.taxa.pc1.density, height = 0.3) %>% # Insert the PC1 plot on top with a height of 0.3
    aplot::insert_right(Fig1a.taxa.pc2.density, width=0.3) %>% # Insert the PC2 plot on the right with a width of 0.3
    ggplotify::as.ggplot() # Convert the combined plot to a ggplot object
}
