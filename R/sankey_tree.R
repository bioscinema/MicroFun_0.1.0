# Shared tree preparation for both Sankey directions. Trees never carry flow weights.
.microfun_taxonomy <- function(physeq, taxonomy, taxon_level) {
  if (is.null(taxonomy) && !is.null(physeq))
    taxonomy <- as(phyloseq::tax_table(physeq), "matrix")
  if (is.null(taxonomy))
    return(NULL)
  tx <- as.data.frame(taxonomy, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(tx) || anyDuplicated(names(tx)))
    stop("taxonomy must have rows and unique rank columns.")
  idcol <- intersect(c("ASV", "FeatureID", "feature_id", "taxon_id"), names(tx))
  ids <- if (length(idcol))
    as.character(tx[[idcol[1]]]) else rownames(tx)
  if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("taxonomy requires unique ASV row names or an ASV/FeatureID column.")
  tx[] <- lapply(tx, as.character)
  rownames(tx) <- ids
  if (toupper(taxon_level) == "ASV")
    tx$.taxon <- ids else {
    k <- match(tolower(taxon_level), tolower(names(tx)))
    if (is.na(k))
      stop("taxonomy does not contain taxon_level: ", taxon_level)
    tx$.taxon <- tx[[k]]
  }
  tx
}

.microfun_taxonomy_tree <- function(taxa, tx, ranks, taxon_level) {
  if (is.null(tx))
    stop("tree_mode='taxonomy' requires taxonomy or physeq.")
  defaults <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  if (is.null(ranks)) {
    end <- match(tolower(taxon_level), tolower(defaults))
    ranks <- if (is.na(end))
      defaults else head(defaults, end - 1L)
    ranks <- ranks[tolower(ranks) %in% tolower(names(tx))]
  }
  idx <- match(tolower(ranks), tolower(names(tx)))
  if (anyNA(idx) || anyDuplicated(idx))
    stop("taxonomy_ranks must name distinct taxonomy columns.")
  ranks <- names(tx)[idx]
  if (any(tolower(ranks) == tolower(taxon_level)))
    ranks <- ranks[tolower(ranks) != tolower(taxon_level)]
  nodes <- data.frame(id = "root", parent = NA_character_, label = "Taxonomy", depth = 0, taxon = NA_character_,
    stringsAsFactors = FALSE)
  # A full path, rather than the displayed label, identifies an ancestor.
  for (taxon in sort(taxa)) {
    rows <- tx[!is.na(tx$.taxon) & tx$.taxon == taxon, ranks, drop = FALSE]
    rows[] <- lapply(rows, function(v) {
      v[is.na(v) | !nzchar(trimws(v))] <- "Unclassified"
      v
    })
    paths <- unique(rows)
    if (nrow(paths) != 1L)
      stop("Ambiguous or absent taxonomy for plotted taxon: ", taxon, ". Use unique taxon labels/full lineage identifiers.")
    parent <- "root"
    for (j in seq_along(ranks)) {
      label <- paths[[j]][1]
      id <- paste0(parent, "/", nchar(ranks[j]), ":", ranks[j], "=", nchar(label), ":", label)
      if (!id %in% nodes$id)
        nodes <- rbind(nodes, data.frame(id = id, parent = parent, label = paste0(ranks[j], ": ",
          label), depth = j, taxon = NA_character_))
      parent <- id
    }
    nodes <- rbind(nodes, data.frame(id = paste0("taxon:", taxon), parent = parent, label = taxon,
      depth = length(ranks) + 1L, taxon = taxon))
  }
  walk <- function(id) {
    kids <- nodes$id[!is.na(nodes$parent) & nodes$parent == id]
    if (!length(kids))
      return(id)
    unlist(lapply(kids, walk), use.names = FALSE)
  }
  leaves <- match(walk("root"), nodes$id)
  list(nodes = nodes, leaves = leaves, order = nodes$taxon[leaves], mapping = data.frame(tip = nodes$label[leaves],
    taxon = nodes$taxon[leaves], tree = "taxonomy"), diagnostics = list(missing_taxa = character(),
    missing_asvs = character()), branch_lengths = FALSE, trees = list(), mode = "taxonomy")
}

.microfun_read_trees <- function(phylo_tree, tree_file, picrust2_dir, physeq) {
  specified <- sum(c(!is.null(phylo_tree), !is.null(tree_file), !is.null(picrust2_dir)))
  if (specified > 1L)
    stop("Supply only one of phylo_tree, tree_file or picrust2_dir.")
  if (!is.null(picrust2_dir)) {
    if (!dir.exists(picrust2_dir))
      stop("picrust2_dir does not exist.")
    files <- list.files(picrust2_dir, pattern = "\\.(tre|tree|nwk|newick)$", full.names = TRUE)
    reduced <- files[basename(files) %in% c("bac_reduced.tre", "arc_reduced.tre")]
    old <- files[basename(files) == "out.tre"]
    if (length(reduced))
      tree_file <- reduced else if (length(old))
      tree_file <- old else {
      if (length(files) != 1L)
        stop("Tree selection is ambiguous. Supply tree_file explicitly; prefer reduced domain trees.")
      tree_file <- files
    }
  }
  if (!is.null(tree_file)) {
    if (!length(tree_file) || any(!file.exists(tree_file)))
      stop("tree_file not found.")
    phylo_tree <- lapply(tree_file, function(f) {
      tr <- ape::read.tree(f)
      if (!inherits(tr, "phylo"))
        stop("Each file must contain one Newick tree: ", f)
      tr
    })
    names(phylo_tree) <- basename(tree_file)
  }
  if (is.null(phylo_tree) && !is.null(physeq))
    phylo_tree <- phyloseq::phy_tree(physeq, errorIfNULL = FALSE)
  if (inherits(phylo_tree, "phylo"))
    phylo_tree <- list(phylogeny = phylo_tree)
  if (!is.list(phylo_tree) || !length(phylo_tree) || !all(vapply(phylo_tree, inherits, logical(1),
    "phylo")))
    stop("Provide a phylo object, list of phylo objects, tree_file, picrust2_dir, or a physeq tree.")
  if (is.null(names(phylo_tree)))
    names(phylo_tree) <- paste0("tree_", seq_along(phylo_tree))
  if (anyNA(names(phylo_tree)) || any(!nzchar(names(phylo_tree))) || anyDuplicated(names(phylo_tree)))
    stop("Tree names must be unique and non-empty.")
  phylo_tree
}

.microfun_phylogeny_order <- function(mapping, method = c("median", "first")) {
  method <- match.arg(method)
  # Order tips within each pruned tree, then concatenate trees in input order.
  mapping$tip_position <- seq_len(nrow(mapping))
  labels <- unique(mapping$taxon)
  positions <- lapply(labels, function(taxon) mapping$tip_position[mapping$taxon == taxon])
  summary <- data.frame(taxon = labels, median_tip_position = vapply(positions, stats::median, numeric(1)),
    first_tip_position = vapply(positions, min, numeric(1)), last_tip_position = vapply(positions,
      max, numeric(1)), n_tips = lengths(positions), stringsAsFactors = FALSE)
  summary$contiguous <- with(summary, last_tip_position - first_tip_position + 1L == n_tips)
  # Break median-position ties with the first tip position.
  key <- if (method == "median")
    summary$median_tip_position else summary$first_tip_position
  summary <- summary[order(key, summary$first_tip_position), , drop = FALSE]
  rownames(summary) <- NULL
  summary$taxon_position <- seq_len(nrow(summary))
  summary$method <- method
  mapping$taxon_position <- match(mapping$taxon, summary$taxon)
  list(order = summary$taxon, mapping = mapping, table = summary, method = method)
}

.microfun_phylogeny <- function(taxa, tx, trees, tip_level, show_branch_length, missing_taxa, phylogeny_order = c("median",
  "first")) {
  phylogeny_order <- match.arg(phylogeny_order)
  if (tip_level == "ASV" && is.null(tx))
    stop("ASV trees require original ASV taxonomy (or use tree_tip_level='taxon' for a taxon-labelled tree).")
  wanted <- if (tip_level == "taxon")
    stats::setNames(taxa, taxa) else stats::setNames(tx$.taxon[!is.na(tx$.taxon) & tx$.taxon %in% taxa], rownames(tx)[!is.na(tx$.taxon) &
    tx$.taxon %in% taxa])
  all_tips <- unlist(lapply(trees, `[[`, "tip.label"), use.names = FALSE)
  duplicate <- unique(all_tips[duplicated(all_tips) & all_tips %in% names(wanted)])
  if (length(duplicate))
    stop("Study tips occur repeatedly or in multiple trees: ", paste(head(duplicate, 8), collapse = ", "),
      ". Use the reduced PICRUSt2 domain trees or resolve domain assignments.")
  missing <- setdiff(taxa, unname(wanted[intersect(names(wanted), all_tips)]))
  if (length(missing) && missing_taxa == "error")
    stop("No matched tree tips for taxa: ", paste(missing, collapse = ", "))
  if (length(missing))
    warning("Dropping taxa without tree tips: ", paste(missing, collapse = ", "), call. = FALSE)
  absent_asvs <- setdiff(names(wanted), all_tips)
  if (length(absent_asvs))
    warning(length(absent_asvs), " requested ASV/taxon IDs absent from tree; see x$microfun_tree$diagnostics.",
      call. = FALSE)
  nodes <- NULL
  mapping <- NULL
  kept <- list()
  leaf_ids <- character()
  used_lengths <- logical()
  for (name in names(trees)) {
    tr <- trees[[name]]
    tips <- intersect(tr$tip.label, names(wanted))
    if (!length(tips))
      next
    if (anyDuplicated(tr$tip.label))
      stop("Tree tip labels must be unique: ", name)
    if (!is.null(tr$edge.length) && (any(!is.finite(tr$edge.length)) || any(tr$edge.length < 0)))
      stop("Tree branch lengths must be finite and non-negative: ", name)
    use_lengths <- show_branch_length && !is.null(tr$edge.length)
    if (show_branch_length && !use_lengths)
      warning("No branch lengths in ", name, "; displaying topology only.", call. = FALSE)
    # Retain the root edge when pruning removes the input root.
    tr <- ape::keep.tip(tr, tips)
    tr <- ape::reorder.phylo(tr, "cladewise")
    kept[[name]] <- tr
    n <- length(tr$tip.label)
    total <- n + tr$Nnode
    ids <- paste0(name, "::", seq_len(total))
    parent <- rep(NA_character_, total)
    parent[tr$edge[, 2]] <- ids[tr$edge[, 1]]
    depth <- numeric(total)
    if (use_lengths)
      depth <- ape::node.depth.edgelength(tr) else {
      for (i in seq_len(nrow(tr$edge))) depth[tr$edge[i, 2]] <- depth[tr$edge[i, 1]] + 1
    }
    label <- c(tr$tip.label, rep("Internal node", tr$Nnode))
    taxon <- c(unname(wanted[tr$tip.label]), rep(NA_character_, tr$Nnode))
    nodes <- rbind(nodes, data.frame(id = ids, parent = parent, label = label, depth = depth, taxon = taxon,
      tree = name))
    # reorder.phylo(cladewise) traverses edges in depth-first tip order.
    ord <- tr$edge[tr$edge[, 2] <= n, 2]
    leaf_ids <- c(leaf_ids, ids[ord])
    mapping <- rbind(mapping, data.frame(tip = tr$tip.label[ord], taxon = taxon[ord], tree = name))
    used_lengths <- c(used_lengths, use_lengths)
  }
  if (is.null(nodes))
    stop("No tree tips match the plotted data.")
  if (length(unique(used_lengths)) > 1L)
    stop("Some trees have no branch lengths. Set show_branch_length=FALSE for a comparable topology-only forest.")
  leaves <- match(leaf_ids, nodes$id)
  # One genus node per result label even when its ASVs occupy several branches.
  ordering <- .microfun_phylogeny_order(mapping, phylogeny_order)
  list(nodes = nodes, leaves = leaves, order = ordering$order, mapping = ordering$mapping, trees = kept,
    taxon_order = ordering$table, order_method = ordering$method, branch_lengths = all(used_lengths),
    mode = "phylogeny", diagnostics = list(missing_taxa = missing, missing_asvs = absent_asvs, removed_tips = setdiff(all_tips,
      names(wanted)), mapped_tips = nrow(mapping), multi_tip_taxa = names(which(table(mapping$taxon) >
      1)), noncontiguous_taxa = ordering$table$taxon[!ordering$table$contiguous], order_method = ordering$method))
}

# Retain maximal pure-taxon clades without inventing inter-clade ancestry.
.microfun_clade_display <- function(tree, collapse = TRUE) {
  n <- tree$nodes
  children <- split(seq_len(nrow(n)), n$parent)
  roots <- which(is.na(n$parent))
  queue <- roots
  cursor <- 1L
  while (cursor <= length(queue)) {
    queue <- c(queue, children[[n$id[queue[cursor]]]])
    cursor <- cursor + 1L
  }
  if (length(queue) != nrow(n) || anyDuplicated(queue))
    stop("Invalid tree structure.")
  pure <- rep(NA_character_, nrow(n))
  members <- vector("list", nrow(n))
  for (i in rev(queue)) {
    kids <- children[[n$id[i]]]
    if (!length(kids)) {
      pure[i] <- n$taxon[i]
      members[[i]] <- n$label[i]
    } else {
      members[[i]] <- unlist(members[kids], use.names = FALSE)
      v <- unique(pure[kids])
      if (length(v) == 1L && !is.na(v))
        pure[i] <- v
    }
  }
  parent <- match(n$parent, n$id)
  terminal <- if (collapse && tree$mode == "phylogeny")
    which(!is.na(pure) & (is.na(parent) | is.na(pure[parent]))) else tree$leaves
  position <- vapply(members[terminal], function(z) min(match(z, tree$mapping$tip)), numeric(1))
  terminal <- terminal[order(position)]
  visible <- logical(nrow(n))
  for (i in terminal) {
    while (!is.na(i) && !visible[i]) {
      visible[i] <- TRUE
      i <- parent[i]
    }
  }
  mapping <- data.frame(clade_id = n$id[terminal], taxon = pure[terminal], n_tips = lengths(members[terminal]),
    stringsAsFactors = FALSE)
  membership <- do.call(rbind, lapply(seq_along(terminal), function(j) data.frame(clade_id = mapping$clade_id[j],
    taxon = mapping$taxon[j], tip = members[[terminal[j]]], stringsAsFactors = FALSE)))
  list(visible = which(visible), terminals = terminal, mapping = mapping, membership = membership)
}


# ---- Merged from sankey_tree_plot.R ----
# SVG coordinates are computed once in R, shared by HTML and vector exports.
.microfun_svg_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

.microfun_tree_sankey <- function(style, taxa_on_left, tree_mode, physeq, taxonomy, phylo_tree, tree_file,
  picrust2_dir, taxon_level, taxonomy_ranks, tree_tip_level, show_branch_length, missing_taxa, tree_options,
  height = NULL, width = NULL, node_color = c(Up = "#F5BB7C", Down = "#B4CBE0"), phylogeny_order = c("median",
    "first"), prepared_tree = NULL) {
  phylogeny_order <- match.arg(phylogeny_order)
  tree_tip_level <- match.arg(tree_tip_level, c("ASV", "taxon"))
  missing_taxa <- match.arg(missing_taxa, c("error", "drop"))
  tx <- .microfun_taxonomy(physeq, taxonomy, taxon_level)
  taxa <- style$nodes$name[style$nodes$kind == "taxa"]
  tree <- if (!is.null(prepared_tree))
    prepared_tree else if (tree_mode == "taxonomy")
    .microfun_taxonomy_tree(taxa, tx, taxonomy_ranks, taxon_level) else .microfun_phylogeny(taxa, tx, .microfun_read_trees(phylo_tree, tree_file, picrust2_dir, physeq),
    tree_tip_level, show_branch_length, missing_taxa, phylogeny_order)
  edges <- style$edges
  edges$taxon <- if (taxa_on_left)
    edges$source else edges$target
  edges$function_name <- if (taxa_on_left)
    edges$target else edges$source
  edges <- edges[edges$taxon %in% tree$order, , drop = FALSE]
  if (!nrow(edges))
    stop("No plotted associations remain after tree matching.")
  if (any(!is.finite(edges$value)) || any(edges$value <= 0) || any(!is.finite(edges$logFC)))
    stop("Tree Sankey requires finite effects and finite positive link weights.")
  if (anyNA(edges$function_name) || any(!nzchar(edges$function_name)))
    stop("Function labels must be non-empty.")
  opts <- list(title = "", subtitle = "", font_size = 16, font_family = "Arial, Helvetica, sans-serif", tree_width = 450, node_width = 18, node_gap = 20,
    show_tip_labels = FALSE, tree_colour = "#6D7785", link_opacity = 0.48, aggregate_links = TRUE,
    collapse_clades = TRUE, bundle_mapping = TRUE, show_footer = FALSE)
  if (!is.list(tree_options) || (length(tree_options) && (is.null(names(tree_options)) || any(!nzchar(names(tree_options))))))
    stop("tree_options must be a named list.")
  unknown <- setdiff(names(tree_options), names(opts))
  if (length(unknown))
    stop("Unknown tree_options: ", paste(unknown, collapse = ", "))
  opts[names(tree_options)] <- tree_options
  for (k in c("font_size", "tree_width", "node_width", "node_gap")) if (length(opts[[k]]) != 1L ||
    !is.numeric(opts[[k]]) || !is.finite(opts[[k]]) || opts[[k]] <= 0)
    stop("Invalid tree_options$", k)
  if (!is.numeric(opts$link_opacity) || length(opts$link_opacity) != 1 || !is.finite(opts$link_opacity) ||
    opts$link_opacity < 0 || opts$link_opacity > 1)
    stop("link_opacity must be between 0 and 1.")
  if (length(opts$font_family) != 1L || is.na(opts$font_family) || !nzchar(opts$font_family))
    stop("font_family must be one non-empty font family.")
  edges <- edges[, setdiff(names(edges), c("taxon_y", "function_y", "pixel_width")), drop = FALSE]
  raw_edges <- edges
  edges$direction <- ifelse(edges$logFC > 0, "Up", "Down")
  if (opts$aggregate_links) {
    key <- paste(match(edges$taxon, unique(edges$taxon)), match(edges$function_name, unique(edges$function_name)),
      edges$direction, sep = ":")
    groups <- split(seq_len(nrow(edges)), factor(key, levels = unique(key)))
    raw_edges$draw_group <- match(key, unique(key))
    edges <- do.call(rbind, lapply(seq_along(groups), function(j) {
      ii <- groups[[j]]
      z <- edges[ii[1], c("source", "target", "taxon", "function_name", "direction"), drop = FALSE]
      z$value <- sum(edges$value[ii])
      z$n_records <- length(ii)
      z$draw_group <- j
      z
    }))
  } else {
    edges$n_records <- 1L
    edges$draw_group <- seq_len(nrow(edges))
    raw_edges$draw_group <- edges$draw_group
  }
  display <- .microfun_clade_display(tree, opts$collapse_clades)
  grDevices::col2rgb(opts$tree_colour)
  taxa <- tree$order
  functions <- sort(unique(edges$function_name))
  nt <- length(taxa)
  nf <- length(functions)
  nl <- length(display$terminals)
  fs <- opts$font_size
  min_width <- opts$tree_width + 1000
  label_lines <- function(labels, chars) max(vapply(labels, function(z) {
    w <- strwrap(z, width = chars)
    if (!length(w)) w <- z
    min(3L, sum(pmax(1L, ceiling(nchar(w)/chars))))
  }, numeric(1)))
  tax_lines <- label_lines(taxa, max(12, floor(255/(fs * 0.58))))
  fun_lines <- label_lines(functions, max(12, floor(265/(fs * 0.58))))
  gap <- max(opts$node_gap, fs * 1.12 * max(tax_lines, fun_lines) + 8)
  header <- 35 + if (nzchar(opts$title))
    38 else 0
  header <- header + if (nzchar(opts$subtitle))
    24 else 0
  top <- header + if (tree_mode == "taxonomy")
    35 else 25
  footer <- if (opts$show_footer)
    110 else 80
  min_height <- max(600, nt * (fs * 1.4 + gap) + top + footer, nf * (fs * 1.4 + gap) + top + footer,
    nl * (if (opts$show_tip_labels) fs + 4 else 9) + top + footer)
  for (v in list(height, width)) if (!is.null(v) && (length(v) != 1L || !is.numeric(v) || !is.finite(v) ||
    v <= 0))
    stop("height/width must be positive numeric pixels in tree modes.")
  W <- max(min_width, if (is.null(width)) 1500 else width)
  H <- max(min_height, if (is.null(height)) min_height else height)
  bottom <- H - footer
  usable <- bottom - top
  tw <- as.numeric(tapply(edges$value, factor(edges$taxon, levels = taxa), sum))
  fw <- as.numeric(tapply(edges$value, factor(edges$function_name, levels = functions), sum))
  if (!is.finite(sum(tw)) || !is.finite(sum(fw)))
    stop("Summed link weights overflow.")
  # Use one weight-to-pixel scale for both sides.
  scale <- min((usable - (nt - 1) * gap)/sum(tw), (usable - (nf - 1) * gap)/sum(fw))
  bands <- function(weights) {
    h <- weights * scale
    extra <- (usable - sum(h) - (length(h) - 1) * gap)/2
    y <- top + extra + c(0, head(cumsum(h + gap), -1))
    data.frame(y = y, h = h, cy = y + h/2)
  }
  raw_edges$pixel_width <- raw_edges$value * scale
  tb <- bands(tw)
  fb <- bands(fw)
  rownames(tb) <- taxa
  rownames(fb) <- functions
  tn <- tree$nodes
  tn$y <- NA_real_
  one_to_one <- !anyDuplicated(tree$mapping$taxon)
  if (one_to_one)
    tn$y[tree$leaves] <- tb[tn$taxon[tree$leaves], "cy"] else tn$y[tree$leaves] <- if (length(tree$leaves) == 1)
    mean(c(top, bottom)) else seq(top + 6, bottom - 6, length.out = length(tree$leaves))
  # Iterative postorder avoids recursion limits on deep ASV trees.
  children <- split(seq_len(nrow(tn)), tn$parent)
  roots <- which(is.na(tn$parent))
  queue <- integer(nrow(tn))
  queue[seq_along(roots)] <- roots
  nqueued <- length(roots)
  cursor <- 1L
  while (cursor <= nqueued) {
    kids <- children[[tn$id[queue[cursor]]]]
    if (length(kids)) {
      queue[nqueued + seq_along(kids)] <- kids
      nqueued <- nqueued + length(kids)
    }
    cursor <- cursor + 1L
  }
  if (nqueued != nrow(tn))
    stop("Tree has disconnected or cyclic nodes.")
  for (i in rev(queue)) {
    kids <- children[[tn$id[i]]]
    if (length(kids))
      tn$y[i] <- mean(range(tn$y[kids]))
  }
  max_depth <- max(tn$depth)
  if (max_depth == 0)
    max_depth <- 1
  tn$x <- 40 + opts$tree_width * tn$depth/max_depth
  tax_x <- opts$tree_width + 350
  fun_x <- W - 300
  nw <- opts$node_width
  flip <- function(x) if (taxa_on_left)
    x else W - x
  E <- .microfun_svg_escape
  fmt <- function(x) sprintf("%.2f", x)
  svg <- character()
  add <- function(s) svg <<- c(svg, s)
  text <- function(x, y, label, size = fs, anchor = "start", fill = "#273444", weight = "normal", tip = NULL) {
    if (!taxa_on_left)
      anchor <- switch(anchor, start = "end", end = "start", middle = "middle")
    add(sprintf("<text x=\"%s\" y=\"%s\" font-size=\"%s\" text-anchor=\"%s\" fill=\"%s\" font-weight=\"%s\">%s%s</text>",
      fmt(flip(x)), fmt(y), fmt(size), anchor, E(fill), weight, if (is.null(tip))
        "" else paste0("<title>", E(tip), "</title>"), E(label)))
  }
  line <- function(x1, y1, x2, y2, col = opts$tree_colour, dash = "", tooltip = NULL) {
    add(sprintf("<path d=\"M %s %s L %s %s\" fill=\"none\" stroke=\"%s\" stroke-width=\"1.3\" %s>%s</path>",
      fmt(flip(x1)), fmt(y1), fmt(flip(x2)), fmt(y2), E(col), if (nzchar(dash))
        paste0("stroke-dasharray=\"", dash, "\"") else "", if (is.null(tooltip))
        "" else paste0("<title>", E(tooltip), "</title>")))
  }
  # Document header is not mirrored.
  add(sprintf("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%s\" height=\"%s\" viewBox=\"0 0 %s %s\" role=\"img\" aria-label=\"%s\"><rect width=\"100%%\" height=\"100%%\" fill=\"white\"/><style>.mf-hover .mf-link{opacity:.08}.mf-hover .mf-link.mf-active{opacity:1}</style><g font-family=\"%s\">",
    fmt(W), fmt(H), fmt(W), fmt(H), E(opts$title), E(opts$font_family)))
  if (nzchar(opts$title))
    add(sprintf("<text x=\"40\" y=\"32\" font-size=\"27\" font-weight=\"bold\" fill=\"#233044\">%s</text>",
      E(opts$title)))
  if (nzchar(opts$subtitle))
    add(sprintf("<text x=\"40\" y=\"%s\" font-size=\"14\" fill=\"#677386\">%s</text>", fmt(header -
      22), E(opts$subtitle)))
  tree_heading <- if (tree_mode == "taxonomy")
    "Taxonomic tree" else "Phylogenetic tree"
  text(40, header, tree_heading, 18, weight = "bold")
  if (tree_mode == "taxonomy") {
    depths <- sort(unique(tn$depth[is.na(tn$taxon) & tn$depth > 0]))
    for (dep in depths) {
      ii <- which(tn$depth == dep)[1]
      rank_label <- sub(":.*", "", tn$label[ii])
      text(tn$x[ii], header + 20, rank_label, fs, "middle", fill = "#677386")
    }
  }
  text(tax_x, header, paste0("Taxon (", taxon_level, ")"), 18, "end", weight = "bold")
  text(fun_x, header, "Function", 18, weight = "bold")
  # Association ribbons first; labels and nodes render above them.
  a_used <- numeric(nt)
  b_used <- numeric(nf)
  ord <- order(match(edges$taxon, taxa), match(edges$function_name, functions), match(edges$direction,
    c("Up", "Down")))
  edges <- edges[ord, , drop = FALSE]
  rownames(edges) <- NULL
  edges$taxon_y <- edges$function_y <- edges$pixel_width <- NA_real_
  for (i in seq_len(nrow(edges))) {
    a <- match(edges$taxon[i], taxa)
    b <- match(edges$function_name[i], functions)
    y1 <- tb$y[a] + a_used[a]
    y2 <- fb$y[b] + b_used[b]
    h <- edges$value[i] * scale
    a_used[a] <- a_used[a] + h
    b_used[b] <- b_used[b] + h
    x1 <- flip(tax_x + nw)
    x2 <- flip(fun_x)
    c1 <- x1 + (x2 - x1) * 0.45
    c2 <- x1 + (x2 - x1) * 0.55
    d <- sprintf("M %.2f %.2f C %.2f %.2f %.2f %.2f %.2f %.2f L %.2f %.2f C %.2f %.2f %.2f %.2f %.2f %.2f Z",
      x1, y1, c1, y1, c2, y2, x2, y2, x2, y2 + h, c2, y2 + h, c1, y1 + h, x1, y1 + h)
    group <- edges$direction[i]
    tooltip <- paste0(edges$taxon[i], " | ", edges$function_name[i], "; ", group, "; records = ",
      edges$n_records[i], "; weight = ", signif(edges$value[i], 5))
    add(sprintf("<path class=\"mf-link\" data-taxon=\"%s\" data-function=\"%s\" d=\"%s\" fill=\"%s\" fill-opacity=\"%s\"><title>%s</title></path>",
      E(edges$taxon[i]), E(edges$function_name[i]), d, E(node_color[group]), fmt(opts$link_opacity),
      E(tooltip)))
    edges$taxon_y[i] <- y1
    edges$function_y[i] <- y2
    edges$pixel_width[i] <- h
  }
  # One path for visible tree branches; full nodes remain in the RDS object.
  segments <- character()
  for (i in display$visible) {
    if (is.na(tn$parent[i]))
      next
    parent_i <- match(tn$parent[i], tn$id)
    segments <- c(segments, sprintf("M %.2f %.2f V %.2f H %.2f", flip(tn$x[parent_i]), tn$y[parent_i],
      tn$y[i], flip(tn$x[i])))
  }
  add(sprintf("<path d=\"%s\" fill=\"none\" stroke=\"%s\" stroke-width=\"1.3\"/>", paste(segments,
    collapse = " "), E(opts$tree_colour)))
  if (tree_mode == "taxonomy") {
    for (i in which(is.na(tn$taxon) & !is.na(tn$parent))) {
      label <- sub("^[^:]+: ", "", tn$label[i])
      limit <- max(5, floor((opts$tree_width/max_depth - 10)/(fs * 0.58)))
      short <- if (nchar(label) > limit)
        paste0(substr(label, 1, max(2, limit - 3)), "...") else label
      text(tn$x[i] + 4, tn$y[i] - 7, short, fs, tip = tn$label[i])
    }
  }
  if (tree_mode == "phylogeny" && length(tree$trees) > 1L) {
    for (i in which(is.na(tn$parent))) text(tn$x[i] + 4, tn$y[i] - 8, tn$tree[i], 11, tip = tn$tree[i])
  }
  # Mapping junctions are annotation devices, never inferred tree ancestors.
  endpoints <- display$terminals
  map_taxa <- display$mapping$taxon
  for (t in taxa) {
    jj <- which(map_taxa == t)
    ii <- endpoints[jj]
    col <- style$nodes$colour[style$nodes$kind == "taxa" & style$nodes$name == t][1]
    hub_x <- tax_x - 290
    hub_y <- tb[t, "cy"]
    bundled <- opts$bundle_mapping && tree_mode == "phylogeny"
    for (j in seq_along(ii)) {
      i <- ii[j]
      count <- display$mapping$n_tips[jj[j]]
      line(tn$x[i], tn$y[i], if (bundled)
        hub_x else tax_x, hub_y, col = col, dash = "3 4")
      tip <- paste0(t, "; ", count, " matched tip", if (count != 1)
        "s" else "")
      if (count > 1L && opts$collapse_clades && tree_mode == "phylogeny") {
        add(sprintf("<path d=\"M %.2f %.2f l %.2f -4 v 8 Z\" fill=\"%s\"><title>%s</title></path>",
          flip(tn$x[i]), tn$y[i], if (taxa_on_left)
          7 else -7, E(col), E(tip)))
      } else add(sprintf("<circle cx=\"%s\" cy=\"%s\" r=\"3\" fill=\"%s\"><title>%s</title></circle>",
        fmt(flip(tn$x[i])), fmt(tn$y[i]), E(col), E(tn$label[i])))
      if (opts$show_tip_labels && tree_mode == "phylogeny") {
        label <- if (count > 1)
          paste0(t, " (", count, ")") else tn$label[i]
        short <- if (nchar(label) > 16)
          paste0(substr(label, 1, 13), "...") else label
        text(tn$x[i] + 8, tn$y[i] - 5, short, fs, tip = label)
      }
    }
    if (bundled) {
      line(hub_x, hub_y, tax_x, hub_y, col = col, dash = "3 4")
      add(sprintf("<rect x=\"%.2f\" y=\"%.2f\" width=\"4\" height=\"4\" fill=\"%s\"><title>%s</title></rect>",
        flip(hub_x) - 2, hub_y - 2, E(col), E(paste(t, "annotation junction"))))
    }
  }
  wrap <- function(s, limit) {
    lines <- strwrap(s, width = limit)
    if (!length(lines))
      lines <- s
    # Also wrap long labels without spaces, rather than overflow into another panel.
    lines <- unlist(lapply(lines, function(z) if (nchar(z) <= limit)
      z else substring(z, seq(1, nchar(z), limit), pmin(nchar(z), seq(1, nchar(z), limit) + limit - 1))))
    if (length(lines) > 3)
      lines <- c(lines[1:2], paste0(substr(lines[3], 1, limit - 3), "..."))
    lines
  }
  plot_nodes <- function(labels, band, x, kind, anchor, label_x) {
    for (i in seq_along(labels)) {
      label <- labels[i]
      col <- style$nodes$colour[style$nodes$kind == kind & style$nodes$name == label][1]
      xleft <- min(flip(x), flip(x + nw))
      add(sprintf("<g class=\"mf-node\" data-%s=\"%s\"><rect x=\"%s\" y=\"%s\" width=\"%s\" height=\"%s\" fill=\"%s\" stroke=\"#FFFFFF\" stroke-width=\"1\"><title>%s</title></rect>",
        if (kind == "taxa")
          "taxon" else "function", E(label), fmt(xleft), fmt(band$y[i]), fmt(nw), fmt(band$h[i]), E(col),
        E(label)))
      lines <- wrap(label, if (kind == "taxa")
        max(12, floor(255/(fs * 0.58))) else max(12, floor(265/(fs * 0.58))))
      if (kind == "taxa") {
        mask_w <- min(258, max(nchar(lines)) * fs * 0.64 + 10)
        mask_left <- min(flip(label_x - mask_w), flip(label_x + 3))
        add(sprintf("<rect x=\"%s\" y=\"%s\" width=\"%s\" height=\"%s\" fill=\"white\" fill-opacity=\"0.94\"/>",
          fmt(mask_left), fmt(band$cy[i] - length(lines) * fs * 0.58 - 3), fmt(mask_w + 3), fmt(length(lines) *
          fs * 1.16 + 6)))
      }
      for (j in seq_along(lines)) {
        yy <- band$cy[i] + (j - (length(lines) + 1)/2) * fs * 1.12 + fs * 0.34
        # White backing keeps mapping lines/ribbons out of label glyphs.
        text(label_x, yy, lines[j], fs, anchor, tip = label)
      }
      add("</g>")
    }
  }
  plot_nodes(taxa, tb, tax_x, "taxa", "end", tax_x - 12)
  plot_nodes(functions, fb, fun_x, "function", "start", fun_x + nw + 12)
  # Footer and scale remain left-to-right regardless of graph direction.
  note <- "Ribbon width: retained association count"
  if (tree_mode == "phylogeny")
    note <- paste0(note, "  |  ", if (tree$branch_lengths)
      "Tree: branch lengths retained" else "Tree: topology only")
  if (!one_to_one)
    note <- paste0(note, "  |  Dashed lines: ASV-to-taxon mapping")
  if (opts$show_footer)
    add(sprintf("<text x=\"40\" y=\"%s\" font-size=\"13\" fill=\"#677386\">%s</text>", fmt(H - 63),
      E(note)))
  for (i in 1:2) {
    xx <- 40 + (i - 1) * 210
    add(sprintf("<rect x=\"%s\" y=\"%s\" width=\"18\" height=\"12\" fill=\"%s\"/><text x=\"%s\" y=\"%s\" font-size=\"13\" fill=\"#455468\">%s</text>",
      xx, fmt(H - 40), E(node_color[c("Up", "Down")[i]]), xx + 26, fmt(H - 29), c("Positive effect",
        "Non-positive effect")[i]))
  }
  if (tree$branch_lengths) {
    bar <- signif(max_depth/4, 1)
    px <- opts$tree_width * bar/max_depth
    line(40, H - 58, 40 + px, H - 58)
    text(40 + px/2, H - 65, paste0(bar, " branch-length units"), 11, "middle")
  }
  add("</g></svg>")
  tree$display_nodes <- tn[display$visible, , drop = FALSE]
  tree$display_mapping <- display$mapping
  tree$clade_membership <- display$membership
  tree$diagnostics$display_tips <- length(display$terminals)
  tree$diagnostics$raw_ribbons <- nrow(raw_edges)
  tree$diagnostics$display_ribbons <- nrow(edges)
  tree$nodes <- tn
  tree$taxon_bands <- cbind(taxon = taxa, tb, weight = tw)
  tree$function_bands <- cbind(function_name = functions, fb, weight = fw)
  tree$diagnostics$weight_sum <- sum(edges$value)
  x <- list(svg = paste(svg, collapse = "\n"), microfun_tree = tree, microfun_edges = raw_edges, microfun_ribbons = edges,
    microfun_colours = style$nodes[c("name", "side", "kind", "category", "colour")], dimensions = list(width = W,
      height = H), settings = list(tree_mode = tree_mode, taxa_on_left = taxa_on_left, node_color = node_color,
      taxon_level = taxon_level, tree_options = opts, phylogeny_order = if (tree_mode == "phylogeny") phylogeny_order else NULL),
    title = opts$title)
  htmlwidgets::createWidget("microfunTreeSankey", x, width = W, height = H, package = "MicroFun", sizingPolicy = htmlwidgets::sizingPolicy(padding = 0,
    browser.fill = FALSE, viewer.fill = FALSE), preRenderHook = function(widget) {
    # HTML receives display data only; RDS retains trees and statistical records.
    widget$x <- widget$x[c("svg", "dimensions", "title")]
    widget
  })
}

#' Export a Sankey (standard or tree-integrated)
#'
#' Standard Sankey supports HTML and RDS; tree Sankey supports HTML, SVG, PDF, PNG and RDS.
#' Tree HTML and SVG use the same fixed layout as PDF/PNG. Hovering over a taxon,
#' function or ribbon highlights its connections in HTML. No browser is required.
#' @param x A widget returned by either Sankey function (any tree_mode).
#' @param file Output path ending in .html, .svg, .pdf, .png or .rds.
#' @param width,height Optional dimensions: pixels for PNG, points for PDF.
#' @param pdf_size Maximum PDF width and height in inches.
#' @param pdf_mode Fit one page or export vertical continuation tiles.
#' @param selfcontained For standard HTML widgets, create one self-contained
#'   HTML file (requires Pandoc); FALSE creates an HTML file with a dependency directory.
#' @return Invisibly, the normalized output file path.
#' @export
save_sankey <- function(x, file, width = NULL, height = NULL, pdf_size = c(11.7, 16.5), pdf_mode = c("fit",
  "tiles"), selfcontained = TRUE) {
  pdf_mode <- match.arg(pdf_mode)
  if (!inherits(x, "htmlwidget") || is.null(x$x))
    stop("Expected a Sankey htmlwidget returned by tax2fun_sankey() or fun2tax_sankey().")
  tree_enabled <- !is.null(x$x$svg) && !is.null(x$x$microfun_tree)
  ext <- tolower(tools::file_ext(file))
  if (!ext %in% c("html", "svg", "pdf", "png", "rds"))
    stop("Supported extensions: html, svg, pdf, png, rds.")
  if (!tree_enabled && !ext %in% c("html", "rds"))
    stop("Standard Sankey (tree_mode='none') supports HTML or RDS only. ",
         "Use tree_mode='taxonomy' or 'phylogeny' for SVG/PDF/PNG.")
  if (!dir.exists(dirname(file)))
    dir.create(dirname(file), recursive = TRUE)
  if (!tree_enabled) {
    if (ext == "rds")
      saveRDS(x, file)
    else
      htmlwidgets::saveWidget(x, file = file, selfcontained = selfcontained)
    return(invisible(normalizePath(file, mustWork = TRUE)))
  }
  if (ext == "svg")
    writeLines(enc2utf8(x$x$svg), file, useBytes = TRUE)
  if (ext == "html") {
    script <- .microfun_sankey_interaction()
    html <- paste0("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\"><meta name=\"viewport\" content=\"width=device-width,initial-scale=1\"><title>",
      .microfun_svg_escape(x$x$title), '</title><style>body{margin:0;background:#edf0f4;font-family:Arial,sans-serif}.canvas{background:white;margin:20px auto;padding:16px;overflow:auto;max-width:1600px}svg{display:block;max-width:100%;height:auto}.mf-node,.mf-link{cursor:crosshair}text{user-select:text}</style></head><body><main class="canvas">', x$x$svg, "</main><script>(",
      script, ")(document.querySelector(\"svg\"));</script></body></html>")
    writeLines(enc2utf8(html), file, useBytes = TRUE)
  }
  if (ext %in% c("pdf", "png")) {
    if (!requireNamespace("rsvg", quietly = TRUE))
      stop("Install 'rsvg' for PDF/PNG export.")
    if (ext == "png")
      rsvg::rsvg_png(charToRaw(enc2utf8(x$x$svg)), file = file, width = width, height = height) else {
      if (!is.numeric(pdf_size) || length(pdf_size) != 2L || any(!is.finite(pdf_size)) || any(pdf_size <=
        0))
        stop("pdf_size must contain positive width and height in inches.")
      target <- pdf_size * 72
      if (!is.null(width))
        target[1] <- width
      if (!is.null(height))
        target[2] <- height
      if (any(!is.finite(target)) || any(target <= 0))
        stop("Invalid PDF dimensions.")
      dims <- unlist(x$x$dimensions[c("width", "height")])
      if (pdf_mode == "fit") {
        factor <- min(target/dims)
        rsvg::rsvg_pdf(charToRaw(enc2utf8(x$x$svg)), file = file, width = dims[1] * factor, height = dims[2] *
          factor)
      } else {
        if (!requireNamespace("qpdf", quietly = TRUE))
          stop("Install 'qpdf' for tiled PDF export.")
        # Consecutive vertical viewports preserve all geometry and weights.
        viewport_height <- target[2] * dims[1]/target[1]
        starts <- seq(0, dims[2] - 1, by = viewport_height)
        temporary <- tempfile("microfun_pdf_")
        dir.create(temporary)
        on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
        pages <- file.path(temporary, paste0(seq_along(starts), ".pdf"))
        for (j in seq_along(starts)) {
          body <- sub("^<svg[^>]*>", "", x$x$svg)
          svg <- paste0("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"", dims[1], "\" height=\"",
          viewport_height, "\" viewBox=\"0 ", starts[j], " ", dims[1], " ", viewport_height,
          "\">", body)
          rsvg::rsvg_pdf(charToRaw(enc2utf8(svg)), file = pages[j], width = target[1], height = target[2])
        }
        qpdf::pdf_combine(pages, output = file)
      }
    }
  }
  if (ext == "rds")
    saveRDS(x, file)
  invisible(normalizePath(file, mustWork = TRUE))
}

.microfun_sankey_interaction <- function() {
  "function(svg){\n  if(!svg)return;\n  const byTaxon=new Map(),byFunction=new Map();\n  const add=(m,k,p)=>{if(!m.has(k))m.set(k,[]);m.get(k).push(p);};\n  svg.querySelectorAll('.mf-link').forEach(p=>{add(byTaxon,p.dataset.taxon,p);add(byFunction,p.dataset.function,p);});\n  let active=[],current=null;\n  const clear=()=>{active.forEach(p=>p.classList.remove('mf-active'));active=[];svg.classList.remove('mf-hover');current=null;};\n  svg.addEventListener('mouseover',ev=>{\n    const el=ev.target.closest('.mf-node,.mf-link');\n    if(!el||el===current)return;\n    clear();current=el;\n    const t=el.dataset.taxon,f=el.dataset.function;\n    active=(t?byTaxon.get(t):byFunction.get(f))||[];\n    if(t&&f)active=active.filter(p=>p.dataset.function===f);\n    svg.classList.add('mf-hover');active.forEach(p=>p.classList.add('mf-active'));\n  });\n  svg.addEventListener('mouseout',ev=>{if(current&&!current.contains(ev.relatedTarget))clear();});\n  svg.addEventListener('mouseleave',clear);\n}"
}

#' Redraw a saved tree Sankey with the current visualization engine
#' @param x A tree-enabled MicroFun widget or an RDS filename.
#' @param tree_options Named appearance options. Titles/footnotes are off by default.
#' @param width,height Optional logical canvas size.
#' @return A newly rendered htmlwidget retaining the original statistical records.
#' @export
redraw_sankey <- function(x, tree_options = list(), width = NULL, height = NULL) {
  if (is.character(x) && length(x) == 1L)
    x <- readRDS(x)
  if (is.null(x$x$microfun_tree) || is.null(x$x$microfun_edges))
    stop("A saved tree Sankey is required.")
  old <- x$x
  settings <- old$settings
  opts <- utils::modifyList(settings$tree_options, list(title = "", subtitle = "", show_footer = FALSE,
    aggregate_links = TRUE, collapse_clades = TRUE, bundle_mapping = TRUE))
  opts <- utils::modifyList(opts, tree_options)
  colours <- settings$node_color
  if (is.null(colours))
    colours <- c(Up = "#F5BB7C", Down = "#B4CBE0")
  .microfun_tree_sankey(style = list(nodes = old$microfun_colours, edges = old$microfun_edges), taxa_on_left = settings$taxa_on_left,
    tree_mode = settings$tree_mode, physeq = NULL, taxonomy = NULL, phylo_tree = NULL, tree_file = NULL,
    picrust2_dir = NULL, taxon_level = settings$taxon_level, taxonomy_ranks = NULL, tree_tip_level = "ASV",
    show_branch_length = old$microfun_tree$branch_lengths, missing_taxa = "error", tree_options = opts,
    height = height, width = width, node_color = colours, phylogeny_order = if (is.null(settings$phylogeny_order))
      "median" else settings$phylogeny_order, prepared_tree = old$microfun_tree)
}
