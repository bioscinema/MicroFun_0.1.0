#' Sankey diagram of taxon -> pathway-class contributions (PICRUSt2 stratified)
#'
#' @description Builds an interactive Sankey diagram (via **networkD3**) showing flows from
#' taxa (left) to pathway superclasses (right), using per-taxon-per-pathway
#' differential results (e.g., PICRUSt2 stratified output summarized by MaAsLin2).
#' Links are colored by effect direction (Up/Down). Node labels for left-side
#' genera can be rendered **to the left of the rectangles** for readability.
#'
#' @param pathway_sdaa_result A named \code{list} where each name is a taxon and each
#'   element is a \code{data.frame} with at least the columns:
#'   \itemize{
#'     \item \code{pathway_class} Character; pathway superclass label (e.g., KEGG level-1).
#'     \item \code{logFC} Numeric; effect sign/size used to color links (>\,0 = Up).
#'     \item \code{p_adjust} Numeric; adjusted p-value (used only for filtering in user code).
#'   }
#'   Rows with \code{NA} in \code{pathway_class} are dropped.
#'
#' @param taxa_list Optional character vector of taxon names to include. Matching is
#'   case-insensitive. If supplied and none match, the function errors.
#'
#' @param node_neutral Optional single colour overriding both node palettes.
#'   NULL (default) uses the light taxa and dark function palettes.
#'
#' @param taxon_groups Optional named character vector mapping each taxon name
#'   to a taxonomy category such as Phylum. Without it, taxa are coloured by name.
#'
#' @param margin A named list with \code{top}, \code{right}, \code{bottom}, \code{left}
#'   margins (pixels) passed to \code{networkD3::sankeyNetwork()}. Increase \code{left}
#'   when \code{label_left_sources = TRUE} to avoid clipping.
#'
#' @importFrom networkD3 sankeyNetwork
#'
#' @importFrom htmlwidgets JS onRender
#'
#' @importFrom jsonlite toJSON
#'
#' @param tree_mode Side tree: 'none' (networkD3), 'taxonomy', or 'phylogeny'.
#'
#' @param physeq Optional phyloseq object with original ASV taxonomy/tree.
#'
#' @param taxonomy Taxonomy data.frame/matrix; ASV row names or ASV/FeatureID column.
#'
#' @param phylo_tree A phylo object or named list of independent phylo trees.
#'
#' @param tree_file Newick file path(s). Each file must contain one tree.
#'
#' @param picrust2_dir PICRUSt2 output folder. Prefers reduced domain trees, then out.tre.
#'
#' @param taxon_level Rank inferred from recorded analysis settings; otherwise Genus. Use 'ASV' for ASV IDs.
#'
#' @param taxonomy_ranks Ordered ranks used for the taxonomic side tree; NULL uses available ancestors.
#'
#' @param tree_tip_level 'ASV' matches tree tips to taxonomy row IDs; 'taxon' requires already taxon-labelled tips.
#'
#' @param show_branch_length Use available phylogenetic branch lengths; otherwise topology only.
#'
#' @param phylogeny_order Taxon ordering in phylogeny mode: 'median' (default)
#'   uses the median displayed ASV tip position; 'first' uses the first tip position.
#'   Median ties use the first tip position. Tree topology determines tip order.
#'   Used only in phylogeny mode.
#'
#' @param missing_taxa 'error' (default) or 'drop' for taxa without matched tree tips.
#'
#' @param height,width Optional tree-mode canvas dimensions in pixels.
#'
#' @export
#'
#' @param node_color Named colours for `Up` and `Down` links. Positive effects are Up; zero or negative effects are Down. Node category colours are set by the palettes.
#'
#' @param tree_options Named list of tree appearance options: `title`, `subtitle`, `font_size`, `font_family`, `tree_width`, `node_width`, `node_gap`, `show_tip_labels`, `tree_colour`, `link_opacity`, `aggregate_links`, `collapse_clades`, `bundle_mapping` and `show_footer`. Default `node_gap = NULL` chooses spacing from label lengths; a positive value sets explicit spacing.
#'
#' @param label_left_sources Logical. In standard mode, place left-node labels to the left and right-node labels to the right. This hook does not wrap labels and does not control tree-mode labels.
#'
#' @param taxon_palette Taxon category colours, optionally named by category. Unnamed palettes are extended when more colours are required.
#'
#' @param function_palette Function category colours, optionally named by category. Unnamed palettes are extended when more colours are required.
#'
#' @param height, width Optional minimum drawing dimensions in pixels for tree modes. Currently ignored in standard `tree_mode = "none"` mode.
#'
#' @return An htmlwidget: a networkD3 widget for `tree_mode = "none"`, or a MicroFun tree Sankey widget for taxonomy/phylogeny mode. Tree widgets retain original edges, tree data and rendering settings for RDS export and redrawing.
#'
#' @details Each retained association contributes one unit of link width; these Sankey widths are association counts, not measured abundance. Links use the sign of `logFC`. No new significance cutoff is applied; filter result tables before plotting if needed. Tree modes integrate a taxonomic or phylogenetic side tree. Named palettes preserve category assignments. Use `save_sankey()` for supported exports.
#'
#' @seealso [save_sankey()], [redraw_sankey()]
#'
#' @md
tax2fun_sankey <- function(pathway_sdaa_result, taxa_list = NULL, node_color = c(Up = "#F5BB7C", Down = "#B4CBE0"),
  node_neutral = NULL, label_left_sources = TRUE, margin = list(top = 10, right = 140, bottom = 10,
    left = 160), taxon_palette = microfun_palette("taxa"), function_palette = microfun_palette("function"),
  taxon_groups = NULL, tree_mode = c("none", "taxonomy", "phylogeny"), physeq = NULL, taxonomy = NULL,
  phylo_tree = NULL, tree_file = NULL, picrust2_dir = NULL, taxon_level = NULL, taxonomy_ranks = NULL,
  tree_tip_level = c("ASV", "taxon"), show_branch_length = TRUE, missing_taxa = c("error", "drop"),
  tree_options = list(), height = NULL, width = NULL, phylogeny_order = c("median", "first")) {
  taxon_level <- .microfun_result_level(pathway_sdaa_result, taxon_level)
  .microfun_result_direction(pathway_sdaa_result, "taxa")
  tree_mode <- match.arg(tree_mode)

  # Build the association table.
  sankey_df <- do.call(rbind, lapply(names(pathway_sdaa_result), function(taxon) {
    df <- pathway_sdaa_result[[taxon]]
    df <- df[!is.na(df$pathway_class), ]
    df$pathway_class <- sub(";.*", "", df$pathway_class)
    if (nrow(df) == 0)
      return(NULL)
    data.frame(source = rep(taxon, nrow(df)), target = df$pathway_class, logFC = df$logFC, qval = df$p_adjust,
      stringsAsFactors = FALSE)
  }))
  if (is.null(sankey_df) || nrow(sankey_df) == 0)
    stop("No taxon-function data.")

  if (!is.null(taxa_list)) {
    taxa_list <- unique(as.character(taxa_list))
    keep_mask <- tolower(sankey_df$source) %in% tolower(taxa_list)
    if (!any(keep_mask))
      stop("None of the taxa in 'taxa_list' are present in the data.")
    sankey_df <- sankey_df[keep_mask, , drop = FALSE]
  }

  sankey_df$group <- ifelse(sankey_df$logFC > 0, "Up", "Down")
  sankey_df$value <- 1
  style <- .microfun_sankey_style(sankey_df, TRUE, taxon_palette, function_palette, taxon_groups, node_color,
    node_neutral)
  if (tree_mode != "none")
    return(.microfun_tree_sankey(style, TRUE, tree_mode, physeq, taxonomy, phylo_tree, tree_file,
      picrust2_dir, taxon_level, taxonomy_ranks, tree_tip_level, show_branch_length, missing_taxa,
      tree_options, height, width, node_color, phylogeny_order))
  nodes <- style$nodes
  sankey_df <- style$edges
  color_scale <- style$scale

  # base widget
  p <- networkD3::sankeyNetwork(Links = sankey_df, Nodes = nodes, Source = "source_id", Target = "target_id",
    Value = "value", NodeID = "name", LinkGroup = "group", NodeGroup = "node_group", colourScale = color_scale,
    fontSize = 18, nodeWidth = 28, sinksRight = TRUE, margin = margin)

  p$x$microfun_colours <- nodes[c("name", "side", "kind", "category", "colour")]

  # move labels: left nodes to left; right nodes to right
  if (label_left_sources) {
    p <- htmlwidgets::onRender(p, "
      function(el, x){
        var svg = d3.select(el).select('svg');
        svg.selectAll('.node text')
          .attr('text-anchor', function(d){
            var isLeft = (d.side === 'left') || (d.x0 === 0 || d.x === 0);
            return isLeft ? 'end' : 'start';
          })
          .attr('x', function(d){
            var w = (d.x1 !== undefined && d.x0 !== undefined) ? (d.x1 - d.x0)
                    : (d.dx !== undefined ? d.dx : 24);
            var isLeft = (d.side === 'left') || (d.x0 === 0 || d.x === 0);
            return isLeft ? -6 : (w + 6);
          });
      }
    ")
  }

  p$x$microfun_taxon_level <- taxon_level
  p <- htmlwidgets::onRender(p, "function(el,x){
    d3.select(el).selectAll('.link').style('stroke-opacity', 0.48);
    d3.select(el).selectAll('.node text').style('font-family','Arial, Helvetica, sans-serif');
  }")
  p
}
