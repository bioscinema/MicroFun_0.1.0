#' Sankey diagram: functions (KO/EC/MetaCyc) to taxa/features
#'
#' Build an interactive Sankey diagram linking functional features (KO, EC, or
#' MetaCyc IDs) — optionally grouped to higher-level categories — to taxa or
#' features reported in differential analysis results (e.g., from
#' \code{pathway_sdaa()}).
#'
#' @description
#' The input \code{pathway_sdaa_result} must be a **named list** where names are
#' function IDs (KO/EC/MetaCyc), and each element is a data frame containing at
#' least \code{feature} (target/taxon), \code{logFC}, and one multiple-testing
#' column (e.g., \code{p_adjust}, \code{padj}, \code{qval}, \code{FDR}). The
#' function:
#' \enumerate{
#'   \item normalizes function IDs (e.g., strips \code{'ko:'}, \code{'ec:'}),
#'   \item (optionally) maps IDs to higher-level groups using reference tables
#'         shipped in \pkg{MicroFun} (\code{inst/extdata/*_reference.RData}),
#'   \item assembles edges \emph{source to target} (group or function → taxon),
#'   \item uses one unit of width per retained association,
#'   \item colors links by direction (Up/Down), and
#'   \item returns an interactive \pkg{networkD3} sankey widget.
#' }
#'
#' @param pathway_sdaa_result Named list. Names are function IDs (KO/EC/MetaCyc).
#'   Each element is a data frame with columns including:
#'   \itemize{
#'     \item \code{feature}: target/taxon (e.g., genus) label,
#'     \item \code{logFC}: signed effect estimate,
#'     \item a q/p-adjust column (any of: \code{p_adjust}, \code{padj},
#'           \code{qval}, \code{q_value}, \code{qvalue}, \code{p.adj}, \code{FDR}, \code{p_adj}).
#'   }
#' @param pathway Character. Functional namespace of the names in
#'   \code{pathway_sdaa_result}; one of \code{'KO'}, \code{'EC'}, \code{'MetaCyc'}.
#' @param group_functions Logical. If \code{TRUE}, map function IDs to higher
#'   groups before linking to taxa (KO to \code{PathwayL1}; EC/MetaCyc to \code{description}).
#'   Requires reference files in \pkg{MicroFun}'s \code{inst/extdata}.
#' @param taxa_list Optional character vector. If provided, keep only edges whose
#'   \emph{source} (group label or function ID) matches one of these entries.
#' @param node_color Named character vector of length 2 with colors for link
#'   groups: names must be \code{c('Up','Down')}. Default: orange for Up
#'   (\code{'#F5BB7C'}), blue for Down (\code{'#B4CBE0'}).
#' @param node_neutral Optional single colour overriding both node palettes.
#'   NULL (default) uses the light taxa and dark function palettes.
#' @param taxon_palette Light taxon colours; optionally named by colour category.
#' @param function_palette Dark function colours; optionally named by category.
#' @param taxon_groups Optional named character vector mapping each taxon name
#'   to a taxonomy category such as Phylum. Without it, taxa are coloured by name.
#' @param label_left_sources Logical. If \code{TRUE}, left-side node labels are
#'   right-aligned and wrapped via JavaScript for readability. Default \code{TRUE}.
#' @param margin Named list with pixel margins \code{top}, \code{right},
#'   \code{bottom}, \code{left}. Passed to \code{networkD3::sankeyNetwork()}.
#' @param height,width Numeric. Output size in pixels for the widget canvas.
#'
#'
#' @examples
#' \dontrun{
#' # Minimal toy example
#' sdaa_list <- list(
#'   K00001 = data.frame(feature = c('GenusA','GenusB'),
#'                       logFC = c( 1.1, -0.4),
#'                       p_adjust = c(0.02, 0.3)),
#'   K00002 = data.frame(feature = c('GenusB','GenusC'),
#'                       logFC = c(-1.6,  0.5),
#'                       p_adjust = c(0.01, 0.7))
#' )
#'
#' # Function to taxon (no grouping)
#' p1 <- fun2tax_sankey(sdaa_list, pathway = 'KO', group_functions = FALSE)
#' p1
#'
#' # Group functions to higher-level categories (requires reference files)
#' p2 <- fun2tax_sankey(sdaa_list, pathway = 'KO', group_functions = TRUE)
#' p2
#' }
#'
#'
#' @importFrom stats setNames
#' @importFrom jsonlite toJSON
#' @importFrom htmlwidgets JS onRender
#' @importFrom networkD3 sankeyNetwork
#' @param tree_mode Side tree: 'none' (networkD3), 'taxonomy', or 'phylogeny'.
#' @param physeq Optional phyloseq object with original ASV taxonomy/tree.
#' @param taxonomy Taxonomy data.frame/matrix; ASV row names or ASV/FeatureID column.
#' @param phylo_tree A phylo object or named list of independent phylo trees.
#' @param tree_file Newick file path(s). Each file must contain one tree.
#' @param picrust2_dir PICRUSt2 output folder. Prefers reduced domain trees, then out.tre.
#' @param taxon_level Rank inferred from recorded analysis settings; otherwise Genus. Use 'ASV' for ASV IDs.
#' @param taxonomy_ranks Ordered ranks used for the taxonomic side tree; NULL uses available ancestors.
#' @param tree_tip_level 'ASV' matches tree tips to taxonomy row IDs; 'taxon' requires already taxon-labelled tips.
#' @param show_branch_length Use available phylogenetic branch lengths; otherwise topology only.
#' @param phylogeny_order Taxon ordering in phylogeny mode: 'median' (default)
#'   uses the median displayed ASV tip position; 'first' uses the first tip position.
#'   Median ties use the first tip position. Tree topology determines tip order.
#'   Used only in phylogeny mode.
#' @param missing_taxa 'error' (default) or 'drop' for taxa without matched tree tips.
#' @param tree_options Named list: title, subtitle, font_size, font_family, tree_width, node_width,
#'   node_gap, show_tip_labels, tree_colour, link_opacity. Also accepts aggregate_links, collapse_clades, bundle_mapping and show_footer.
#' @export
fun2tax_sankey <- function(pathway_sdaa_result, pathway = c("KO", "EC", "MetaCyc"), group_functions = TRUE,
  taxa_list = NULL, node_color = c(Up = "#F5BB7C", Down = "#B4CBE0"), node_neutral = NULL, label_left_sources = TRUE,
  margin = list(top = 10, right = 200, bottom = 10, left = 180), height = 1200, width = 1000, taxon_palette = microfun_palette("taxa"),
  function_palette = microfun_palette("function"), taxon_groups = NULL, tree_mode = c("none", "taxonomy",
    "phylogeny"), physeq = NULL, taxonomy = NULL, phylo_tree = NULL, tree_file = NULL, picrust2_dir = NULL,
  taxon_level = NULL, taxonomy_ranks = NULL, tree_tip_level = c("ASV", "taxon"), show_branch_length = TRUE,
  missing_taxa = c("error", "drop"), tree_options = list(), phylogeny_order = c("median", "first")) {
  taxon_level <- .microfun_result_level(pathway_sdaa_result, taxon_level)
  .microfun_result_direction(pathway_sdaa_result, "function")
  tree_mode <- match.arg(tree_mode)
  stopifnot(is.list(pathway_sdaa_result))
  pathway <- match.arg(pathway)

  # Normalize identifiers and load reference tables.
  normalize_fun_ids <- function(ids, pathway = c("KO", "EC", "MetaCyc")) {
    pathway <- match.arg(pathway)
    ids <- trimws(as.character(ids))
    if (pathway == "KO") {
      ids <- toupper(gsub("^ko:", "", ids, ignore.case = TRUE))
    } else if (pathway == "EC") {
      ids <- gsub("^ec:", "", ids, ignore.case = TRUE)
      ids <- gsub("^EC[\\s_]*", "", ids, ignore.case = TRUE)
    } else {
      ids <- gsub("^metacyc:", "", ids, ignore.case = TRUE)
    }
    ids
  }

  load_reference_map <- function(pathway) .microfun_reference_map(pathway, "PathwayL1")

  get_qcol <- function(df) {
    cand <- c("p_adjust", "padj", "qval", "q_value", "qvalue", "p.adj", "FDR", "p_adj")
    hit <- intersect(cand, names(df))
    if (length(hit) == 0)
      return(NULL)
    hit[1]
  }

  # Normalize function identifiers.
  fun_ids_raw <- names(pathway_sdaa_result)
  if (is.null(fun_ids_raw))
    stop("pathway_sdaa_result must be a *named* list (names = function IDs).")

  fun_ids_norm <- normalize_fun_ids(fun_ids_raw, pathway)
  keep <- !is.na(fun_ids_norm) & nzchar(fun_ids_norm)
  if (!any(keep))
    stop("After normalization, no valid function IDs remain.")
  pathway_sdaa_result <- pathway_sdaa_result[keep]
  fun_ids_norm <- fun_ids_norm[keep]
  names(pathway_sdaa_result) <- fun_ids_norm

  # (Optional) load id -> group mapping
  id2group <- NULL
  if (isTRUE(group_functions)) {
    id2group <- load_reference_map(pathway)
  }

  # Build function-to-taxon associations.
  if (isTRUE(group_functions)) {
    # group by higher-level label first
    grp_lab <- unname(id2group[fun_ids_norm])
    valid_fun <- !is.na(grp_lab) & nzchar(grp_lab)
    fun_ids_norm <- fun_ids_norm[valid_fun]
    grp_lab <- grp_lab[valid_fun]
    if (!length(fun_ids_norm))
      stop("No function IDs matched the reference mapping after normalization.")
    groups <- split(fun_ids_norm, grp_lab)

    pieces <- lapply(names(groups), function(glab) {
      fids <- groups[[glab]]
      # collect rows from all member functions
      rows_list <- lapply(fids, function(fid) {
        df <- pathway_sdaa_result[[fid]]
        if (is.null(df) || !is.data.frame(df) || !nrow(df))
          return(NULL)
        qcol <- get_qcol(df)
        if (is.null(qcol) || !"feature" %in% names(df) || !"logFC" %in% names(df))
          return(NULL)
        keep <- !is.na(df$feature) & nzchar(df$feature) & !is.na(df[[qcol]])
        if (!any(keep))
          return(NULL)
        data.frame(function_id = fid, feature = as.character(df$feature[keep]), logFC = df$logFC[keep],
          p_adjust = df[[qcol]][keep], stringsAsFactors = FALSE)
      })
      rows_list <- Filter(Negate(is.null), rows_list)
      if (!length(rows_list))
        return(NULL)
      rows <- do.call(rbind, rows_list)
      if (!nrow(rows))
        return(NULL)

      # Count each function/taxon association within the group once.
      rows <- rows[order(rows$feature, rows$p_adjust, decreasing = FALSE), , drop = FALSE]
      rows <- rows[!duplicated(rows[c("function_id", "feature")]), , drop = FALSE]

      data.frame(source = rep(glab, nrow(rows)), target = rows$feature, logFC = rows$logFC, qval = rows$p_adjust,
        stringsAsFactors = FALSE)
    })

    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces))
      stop("No edges produced after grouping by functions.")
    sankey_df <- do.call(rbind, pieces)

  } else {
    # no grouping: function to taxon edges directly
    pieces <- lapply(names(pathway_sdaa_result), function(fid) {
      df <- pathway_sdaa_result[[fid]]
      if (is.null(df) || !is.data.frame(df) || !nrow(df))
        return(NULL)
      qcol <- get_qcol(df)
      if (is.null(qcol) || !"feature" %in% names(df) || !"logFC" %in% names(df))
        return(NULL)
      keep <- !is.na(df$feature) & nzchar(df$feature) & !is.na(df[[qcol]])
      if (!any(keep))
        return(NULL)
      data.frame(source = rep(fid, sum(keep)), target = as.character(df$feature[keep]), logFC = df$logFC[keep],
        qval = df[[qcol]][keep], stringsAsFactors = FALSE)
    })
    pieces <- Filter(Negate(is.null), pieces)
    if (!length(pieces))
      stop("No edges produced in function to taxon mode.")
    sankey_df <- do.call(rbind, pieces)
  }

  if (nrow(sankey_df) == 0)
    stop("No edges to plot after assembling sankey data.")

  # Filter source functions.
  if (!is.null(taxa_list)) {
    keep_sources <- unique(as.character(taxa_list))
    mask <- tolower(sankey_df$source) %in% tolower(keep_sources)
    if (!any(mask))
      stop("None of the requested sources in 'taxa_list' are present.")
    sankey_df <- sankey_df[mask, , drop = FALSE]
  }

  # Set link directions and count each retained association once.
  sankey_df$link_group <- ifelse(sankey_df$logFC > 0, "Up", "Down")
  sankey_df$value <- 1
  style <- .microfun_sankey_style(sankey_df, FALSE, taxon_palette, function_palette, taxon_groups,
    node_color, node_neutral)
  if (tree_mode != "none")
    return(.microfun_tree_sankey(style, FALSE, tree_mode, physeq, taxonomy, phylo_tree, tree_file,
      picrust2_dir, taxon_level, taxonomy_ranks, tree_tip_level, show_branch_length, missing_taxa,
      tree_options, height, width, node_color, phylogeny_order))
  nodes <- style$nodes
  sankey_df <- style$edges
  color_scale <- style$scale

  p <- networkD3::sankeyNetwork(Links = sankey_df, Nodes = nodes, Source = "source_id", Target = "target_id",
    Value = "value", NodeID = "name", LinkGroup = "link_group", NodeGroup = "node_group", colourScale = color_scale,
    fontSize = 18, nodeWidth = 30, sinksRight = TRUE, margin = margin, height = height, width = width)

  p$x$microfun_colours <- nodes[c("name", "side", "kind", "category", "colour")]

  if (label_left_sources) {
    p <- htmlwidgets::onRender(p, "
    function(el, x){
      var fontSize = 18,  // keep in sync with R
          maxWidth = 220; // wrap width in px for LEFT labels

      var svg = d3.select(el).select('svg');

      // set anchors & x offsets (as you had)
      svg.selectAll('.node text')
        .attr('text-anchor', function(d){
          var isLeft = (d.side === 'left') || (d.x0 === 0 || d.x === 0);
          return isLeft ? 'end' : 'start';
        })
        .attr('x', function(d){
          var w = (d.x1!==undefined && d.x0!==undefined) ? (d.x1 - d.x0)
                  : (d.dx!==undefined ? d.dx : 24);
          var isLeft = (d.side === 'left') || (d.x0 === 0 || d.x === 0);
          return isLeft ? -10 : (w + 10);
        });

      // wrap LEFT labels only
      svg.selectAll('.node')
        .filter(function(d){ return (d.side === 'left') || (d.x0 === 0 || d.x === 0); })
        .select('text')
        .each(function(d){
          var text = d3.select(this),
              words = (text.text() || '').split(/\\s+/).filter(Boolean),
              line = [], lineNumber = 0,
              lineHeight = 1.2, // ems
              y = (d.dy !== undefined) ? d.dy / 2 : 0,
              x = +text.attr('x'),
              anchor = text.attr('text-anchor');

          text.text(null);
          var tspan = text.append('tspan')
                          .attr('x', x)
                          .attr('y', y)
                          .attr('dy', 0 + 'em')
                          .attr('text-anchor', anchor);

          for (var i=0; i<words.length; ++i) {
            line.push(words[i]);
            tspan.text(line.join(' '));
            if (tspan.node().getComputedTextLength() > maxWidth) {
              line.pop();
              tspan.text(line.join(' '));
              line = [words[i]];
              ++lineNumber;
              tspan = text.append('tspan')
                          .attr('x', x)
                          .attr('y', y)
                          .attr('dy', lineNumber * lineHeight + 'em')
                          .attr('text-anchor', anchor)
                          .text(words[i]);
            }
          }
        });

      // tooltip with full label on hover
      svg.selectAll('.node').append('title').text(function(d){ return d.name; });
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
