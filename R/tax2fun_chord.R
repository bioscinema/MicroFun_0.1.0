#' Chord diagram of taxon-resolved functional contributions
#'
#' Connects stratified PICRUSt2 contributions to a named pathway_sdaa result.
#' Width represents contribution, colour represents effect direction, and the
#' outer bars represent absolute model effect (not necessarily log fold change).
#' The returned object is redrawn by print() or plot(); it is not a ggplot.
#' @param pathway_sdaa_result Named list returned by pathway_sdaa().
#' @param stratified_path PICRUSt2 contribution TSV or TSV.GZ file.
#' @param physeq phyloseq object containing taxonomy and sample metadata.
#' @param group Metadata column used in the differential analysis.
#' @param taxon_level Taxonomy rank; inferred from recorded settings, otherwise Genus.
#' @param direction Result orientation: taxa or function. Inferred if recorded.
#' @param function_level auto, KEGG, KO, EC or MetaCyc. KEGG uses PICRUSt2
#'   pathway contributions directly or sums raw contributions of member KOs.
#' @param contrast Two group labels, c(reference, target). Positive effects
#'   indicate target relative to reference; must match the result contrast.
#' @param taxa_list Optional taxon names, matched case-insensitively, exactly.
#' @param functions Optional functional IDs to retain.
#' @param padj_cutoff Retain finite p_adjust <= cutoff. Default 0.05.
#' @param min_abs_effect Minimum absolute model effect.
#' @param max_taxa_per_function Maximum links per function, ranked by width.
#' @param max_functions_per_taxon Maximum links per taxon, ranked by width.
#' @param contribution balanced_mean (equal weight to two group means) or mean
#'   (mean over all selected samples). Uses original contribution units.
#' @param function_table Optional unique function_id/function_group table, with an optional function_name column. Otherwise
#'   missing and inconsistent annotations are reconciled; genuine conflicts are
#'   labelled Ambiguous annotation and retained in observed_groups.
#' @param output_pdf,output_png,output_html Optional export paths. No file is written by default.
#' @param draw Draw on the current device when no export paths are supplied.
#' @param ... Graphical arguments to the internal renderer; see Details in help.
#' @return A microfun_chord object with links (after filtering), all_links (before
#'   top-N), tax_table, function_table, settings and diagnostics.
#' @export
tax2fun_chord <- function(pathway_sdaa_result, stratified_path, physeq, group, taxon_level = NULL, direction = NULL,
  function_level = c("auto", "KEGG", "KO", "EC", "MetaCyc"), contrast = NULL, taxa_list = NULL, functions = NULL,
  padj_cutoff = 0.05, min_abs_effect = 0, max_taxa_per_function = Inf, max_functions_per_taxon = Inf,
  contribution = c("balanced_mean", "mean"), function_table = NULL, output_pdf = NULL, output_png = NULL,
  draw = TRUE, output_html = NULL, ...) {
  contribution <- match.arg(contribution)
  function_level <- match.arg(function_level)
  info <- attr(pathway_sdaa_result, "microfun")
  taxon_level <- .microfun_result_level(pathway_sdaa_result, taxon_level)
  if (is.null(direction))
    direction <- if (!is.null(info$direction))
      info$direction else "taxa"
  direction <- match.arg(direction, c("taxa", "function"))
  for (nm in c("taxon_level", "direction", "group")) {
    if (!is.null(info[[nm]]) && !identical(as.character(info[[nm]]), as.character(get(nm))))
      stop(nm, " differs from the settings recorded in pathway_sdaa_result.")
  }
  if (length(padj_cutoff) != 1L || !is.finite(padj_cutoff) || padj_cutoff < 0 || padj_cutoff > 1)
    stop("padj_cutoff must be between 0 and 1.")
  if (length(min_abs_effect) != 1L || !is.finite(min_abs_effect) || min_abs_effect < 0)
    stop("min_abs_effect must be finite and non-negative.")
  for (n in c(max_taxa_per_function, max_functions_per_taxon)) if (is.na(n) || n < 1 || (is.finite(n) &&
    n != floor(n)))
    stop("Top-N limits must be positive integers or Inf.")
  if (!is.list(pathway_sdaa_result) || !length(pathway_sdaa_result) || is.null(names(pathway_sdaa_result)) ||
    anyDuplicated(names(pathway_sdaa_result)))
    stop("pathway_sdaa_result must be a non-empty, uniquely named list.")
  meta <- data.frame(phyloseq::sample_data(physeq), check.names = FALSE)
  tx <- as.data.frame(as(phyloseq::tax_table(physeq), "matrix"), stringsAsFactors = FALSE)
  if (!group %in% names(meta) || !taxon_level %in% names(tx))
    stop("group or taxon_level is absent from physeq.")
  if (!is.null(info$samples) && !setequal(info$samples, rownames(meta)))
    stop("physeq samples differ from the recorded analysis samples. Use the same subset.")
  groups <- as.character(meta[[group]])
  if (is.null(contrast)) {
    contrast <- if (is.factor(meta[[group]]))
      levels(droplevels(meta[[group]])) else sort(unique(groups[!is.na(groups)]))
    if (length(contrast) != 2L)
      stop("Supply contrast = c(reference, target) for a multi-group dataset.")
  }
  if (length(contrast) != 2L || anyNA(contrast) || anyDuplicated(contrast) || !all(contrast %in% groups))
    stop("contrast must contain two different observed group labels.")
  contrast <- as.character(contrast)
  selected <- !is.na(groups) & groups %in% contrast
  meta <- meta[selected, , drop = FALSE]
  groups <- groups[selected]
  sample_ids <- rownames(meta)
  n_group <- as.numeric(table(factor(groups, levels = contrast)))

  pieces <- lapply(names(pathway_sdaa_result), function(nm) {
    z <- as.data.frame(pathway_sdaa_result[[nm]])
    if (!nrow(z))
      return(NULL)
    needed <- c("feature", "p_adjust", "group1", "group2")
    if (!all(needed %in% names(z)))
      stop("Result '", nm, "' lacks: ", paste(setdiff(needed, names(z)), collapse = ", "))
    z <- z[!is.na(z$group1) & !is.na(z$group2) & z$group1 == contrast[1] & z$group2 == contrast[2],
      , drop = FALSE]
    if (!nrow(z))
      return(NULL)
    effect_col <- intersect(c("effect", "logFC", "log2FoldChange"), names(z))
    if (!length(effect_col))
      stop("Result '", nm, "' has no effect estimate. This backend cannot colour a chord diagram.")
    eff <- suppressWarnings(as.numeric(z[[effect_col[1]]]))
    q <- suppressWarnings(as.numeric(z$p_adjust))
    typ <- if ("effect_type" %in% names(z))
      as.character(z$effect_type) else if ("method" %in% names(z) && all(z$method == "Maaslin2"))
      rep("Maaslin2 AST coefficient", nrow(z)) else rep(effect_col[1], nrow(z))
    data.frame(taxon = if (direction == "taxa")
      nm else as.character(z$feature), function_id = if (direction == "taxa")
      as.character(z$feature) else nm, effect = eff, effect_type = typ, p_adjust = q, function_group = if ("pathway_class" %in%
      names(z))
      sub(";.*", "", z$pathway_class) else NA_character_, stringsAsFactors = FALSE)
  })
  links <- do.call(rbind, pieces)
  if (is.null(links) || !nrow(links))
    stop("No results match contrast. Check group1/group2 and reference order.")
  if (function_level == "auto") {
    function_level <- if (!is.null(info$function_level))
      info$function_level else .microfun_function_level(links$function_id)
  }
  if (!is.null(info$function_level) && !identical(function_level, info$function_level))
    stop("function_level differs from the differential analysis.")
  links$function_id <- .chord_ids(links$function_id, function_level)
  if (anyDuplicated(links[c("taxon", "function_id")]))
    stop("Duplicate taxon/function edges. Select one method and contrast before plotting.")
  count_input <- nrow(links)
  links <- links[is.finite(links$effect) & is.finite(links$p_adjust) & links$p_adjust >= 0 & links$p_adjust <=
    padj_cutoff & abs(links$effect) >= min_abs_effect, , drop = FALSE]
  if (!is.null(taxa_list))
    links <- links[tolower(links$taxon) %in% tolower(taxa_list), , drop = FALSE]
  if (!is.null(functions))
    links <- links[links$function_id %in% .chord_ids(functions, function_level), , drop = FALSE]
  if (!nrow(links))
    stop("No edges remain after effect, significance and selection filters.")

  # Keep one lineage per taxon, with only ranks above the selected level.
  ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  rank_i <- match(taxon_level, ranks)
  if (is.na(rank_i))
    stop("Use a standard taxon_level: ", paste(ranks, collapse = ", "))
  ancestors <- intersect(ranks[seq_len(rank_i - 1L)], names(tx))
  ancestors <- setdiff(ancestors, "Kingdom")
  tax_anno <- data.frame(taxon_id = unique(links$taxon), stringsAsFactors = FALSE)
  for (rank in ancestors) {
    tax_anno[[rank]] <- vapply(tax_anno$taxon_id, function(tn) {
      x <- unique(as.character(tx[!is.na(tx[[taxon_level]]) & tx[[taxon_level]] == tn, rank]))
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x) > 1L)
        stop("Ambiguous lineage for ", tn, " at ", rank, ". Resolve duplicated taxon names first.")
      if (length(x))
        x else "Unclassified"
    }, character(1))
  }
  if (!all(links$taxon %in% tx[[taxon_level]]))
    stop("Some result taxa do not match physeq at taxon_level.")

  # Read compressed or uncompressed contribution tables.
  raw <- if (grepl("\\.gz$", stratified_path, ignore.case = TRUE) && !requireNamespace("R.utils", quietly = TRUE)) {
    con <- gzfile(stratified_path, "rt")
    on.exit(close(con), add = TRUE)
    data.table::as.data.table(utils::read.delim(con, check.names = FALSE))
  } else data.table::fread(stratified_path)
  map_ko <- function_level == "KEGG" && .microfun_kegg_uses_ko(raw[["function"]], info$pathway)
  value_col <- .microfun_contribution_column(names(raw), info$value_col, require_abundance = function_level == "KEGG")
  req <- c("sample", "taxon", "function", value_col)
  if (!all(req %in% names(raw)))
    stop("Stratified table must contain: ", paste(req, collapse = ", "))
  if (!all(sample_ids %in% raw$sample))
    stop("Some selected samples are completely absent from the stratified file; cannot assume these are zero.")
  raw <- as.data.frame(raw[, req, with = FALSE])
  raw <- raw[raw$sample %in% sample_ids, , drop = FALSE]
  raw$taxon <- as.character(tx[match(raw$taxon, rownames(tx)), taxon_level])
  raw <- raw[!is.na(raw$taxon) & raw$taxon %in% links$taxon, , drop = FALSE]
  value <- suppressWarnings(as.numeric(raw[[value_col]]))
  if (any(!is.finite(value) | value < 0))
    stop("Selected contributions must be finite and non-negative.")
  raw$value <- value
  if (map_ko) {
    raw$ko_id <- .microfun_ko_ids(raw[["function"]])
    map <- .microfun_kegg_links(.microfun_kegg_members(info$kegg_map), links$function_id)
    raw <- merge(raw, map, by = "ko_id", sort = FALSE)
    raw$function_id <- raw$pathway_id
  } else {
    raw$function_id <- .chord_ids(raw[["function"]], function_level)
  }
  raw <- merge(raw, unique(links[c("taxon", "function_id")]), by = c("taxon", "function_id"), sort = FALSE)
  if (!nrow(raw))
    stop("No stratified contributions match the selected taxon/function edges.")
  raw$group_id <- match(groups[match(raw$sample, sample_ids)], contrast)
  sums <- stats::aggregate(value ~ taxon + function_id + group_id, raw, sum)
  out <- links
  for (g in 1:2) {
    part <- sums[sums$group_id == g, c("taxon", "function_id", "value"), drop = FALSE]
    names(part)[3] <- paste0("mean_group", g)
    part[[3]] <- part[[3]]/n_group[g]
    out <- merge(out, part, by = c("taxon", "function_id"), all.x = TRUE, sort = FALSE)
    out[[paste0("mean_group", g)]][is.na(out[[paste0("mean_group", g)]])] <- 0
  }
  out$contribution <- if (contribution == "balanced_mean")
    (out$mean_group1 + out$mean_group2)/2 else (out$mean_group1 * n_group[1] + out$mean_group2 * n_group[2])/sum(n_group)
  out$edge_direction <- ifelse(out$effect > 0, "Up", ifelse(out$effect < 0, "Down", "Neutral"))
  zero_width <- sum(out$contribution <= 0)
  out <- out[out$contribution > 0, , drop = FALSE]
  if (!nrow(out))
    stop("All retained edges have zero contribution.")
  all_links <- out
  # Deterministic sequential top-N (function first, then taxon).
  out <- out[order(-out$contribution, out$taxon, out$function_id), , drop = FALSE]
  if (is.finite(max_taxa_per_function))
    out <- out[ave(seq_len(nrow(out)), out$function_id, FUN = seq_along) <= max_taxa_per_function,
      , drop = FALSE]
  if (is.finite(max_functions_per_taxon))
    out <- out[ave(seq_len(nrow(out)), out$taxon, FUN = seq_along) <= max_functions_per_taxon, ,
      drop = FALSE]
  tax_anno <- tax_anno[tax_anno$taxon_id %in% out$taxon, , drop = FALSE]
  if (is.null(function_table)) {
    # Resolve from all significant positive-width edges, before Top-N, so overview filtering cannot
    # change a function's displayed category.
    function_table <- .chord_function_groups(all_links)
    missing <- function_table$annotation_status == "missing"
    if (any(missing) && function_level == "KO") {
      env <- new.env(parent = emptyenv())
      ref <- .microfun_chord_reference("KO_reference")
      if (all(c("KO", "PathwayL1") %in% names(ref))) {
        ref_ids <- .chord_ids(ref$KO, "KO")
        needed <- ref_ids %in% function_table$function_id[missing]
        reference_groups <- .chord_function_groups(data.frame(function_id = ref_ids[needed],
          function_group = as.character(ref$PathwayL1[needed])))
        idx <- match(function_table$function_id, reference_groups$function_id)
        use <- which(missing & !is.na(idx))
        function_table[use, c("function_group", "observed_groups", "annotation_status")] <- reference_groups[idx[use],
          c("function_group", "observed_groups", "annotation_status")]
      }
    }
    function_table <- function_table[function_table$function_id %in% out$function_id, , drop = FALSE]
    conflicts <- function_table$annotation_status == "conflicting"
    if (any(conflicts))
      warning(sum(conflicts), " function(s) have conflicting annotations; displayed as Ambiguous annotation. ",
        "Inspect $function_table$observed_groups or supply function_table.", call. = FALSE)
  } else {
    if (!all(c("function_id", "function_group") %in% names(function_table)))
      stop("function_table needs function_id and function_group.")
    function_table$function_id <- .chord_ids(function_table$function_id, function_level)
    if (anyDuplicated(function_table$function_id) || !all(out$function_id %in% function_table$function_id))
      stop("function_table must match every function exactly once.")
    if (anyNA(function_table$function_group) || any(!nzchar(function_table$function_group)))
      stop("function_group must not be missing.")
    function_table <- function_table[function_table$function_id %in% out$function_id, , drop = FALSE]
  }
  function_table$function_group <- sub(";.*", "", function_table$function_group)
  if (!"function_name" %in% names(function_table)) function_table$function_name <- NA_character_
  inferred_names <- .microfun_function_names(pathway_sdaa_result, direction, function_level, function_table$function_id)
  blank <- is.na(function_table$function_name) | !nzchar(trimws(function_table$function_name))
  function_table$function_name[blank] <- unname(inferred_names[function_table$function_id[blank]])
  settings <- list(taxon_level = taxon_level, direction = direction, function_level = function_level,
    contrast = contrast, contribution = contribution, sample_ids = sample_ids, group_sizes = stats::setNames(n_group,
      contrast), padj_cutoff = padj_cutoff, effect_type = unique(out$effect_type), graphical = list(...))
  if (function_level == "KEGG")
    settings$width_label <- paste("Width:", contribution, if (map_ko)
      "member KO contribution" else "PICRUSt2 pathway contribution")
  if (length(settings$effect_type) != 1L)
    stop("Cannot compare different effect scales in one chord diagram.")
  ans <- structure(list(links = out, all_links = all_links, tax_table = tax_anno, function_table = function_table,
    settings = settings, diagnostics = list(input_edges = count_input, zero_width_edges = zero_width,
      plotted_edges = nrow(out), retained_width_fraction = sum(out$contribution)/sum(all_links$contribution))),
    class = "microfun_chord")
  if (!is.null(output_html))
    save_chord(ans, output_html)
  if (!is.null(output_pdf))
    save_chord(ans, output_pdf)
  if (!is.null(output_png))
    save_chord(ans, output_png)
  if (isTRUE(draw) && is.null(output_pdf) && is.null(output_png) && is.null(output_html))
    print(ans)
  invisible(ans)
}

.chord_ids <- function(x, level) {
  x <- trimws(as.character(x))
  switch(level, KO = toupper(sub("^ko:", "", x, ignore.case = TRUE)), KEGG = .microfun_kegg_ids(x),
    EC = sub("^(EC|ec)[:_ -]*", "", x), MetaCyc = sub("^metacyc:", "", x, ignore.case = TRUE), x)
}

.microfun_function_level <- function(ids) {
  ids <- trimws(ids[!is.na(ids) & nzchar(ids)])
  if (length(ids) && all(grepl("^K[0-9]{5}$", .microfun_ko_ids(ids))))
    "KO" else if (length(ids) && all(grepl("^ko[0-9]{5}$", .microfun_kegg_ids(ids))))
    "KEGG" else if (length(ids) && all(grepl("^[0-9-]+\\.[0-9-]+\\.[0-9-]+\\.[0-9-]+$", .chord_ids(ids, "EC"))))
    "EC" else "MetaCyc"
}

.microfun_kegg_uses_ko <- function(ids, pathway) {
  if (!is.null(pathway)) return(pathway == "KO")
  any(grepl("^K[0-9]{5}$", .microfun_ko_ids(ids)))
}

#' @rdname tax2fun_chord
#' @param x A microfun_chord object.
#' @export
print.microfun_chord <- function(x, ...) {
  opts <- utils::modifyList(x$settings$graphical, list(...))
  rendered <- do.call(.draw_tax2fun_chord, c(list(x = x), opts))
  if (length(rendered$legend_pages))
    .microfun_chord_legend_pages(rendered$legend_pages)
  invisible(rendered)
}

#' @rdname tax2fun_chord
#' @export
plot.microfun_chord <- function(x, ...) print(x, ...)


# Collapse hierarchical annotations to one auditable display category per ID.
.chord_function_groups <- function(x) {
  ids <- trimws(as.character(x$function_id))
  groups <- trimws(sub(";.*$", "", as.character(x$function_group)))
  canonical <- c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
    "Cellular Processes", "Organismal Systems", "Human Diseases")
  hit <- match(tolower(groups), tolower(canonical))
  groups[!is.na(hit)] <- canonical[hit[!is.na(hit)]]
  groups[is.na(groups) | tolower(groups) %in% c("", "na", "n/a", "null", "unclassified function")] <- NA_character_
  valid <- !is.na(ids) & nzchar(ids)
  records <- split(groups[valid], ids[valid])
  result <- lapply(names(records), function(id) {
    labels <- sort(unique(records[[id]][!is.na(records[[id]])]))
    n <- length(labels)
    data.frame(function_id = id, function_group = if (n == 0L)
      "Unclassified function" else if (n == 1L)
      labels else "Ambiguous annotation", observed_groups = paste(labels, collapse = " | "), annotation_status = if (n ==
      0L)
      "missing" else if (n == 1L)
      "resolved" else "conflicting", stringsAsFactors = FALSE)
  })
  if (!length(result))
    return(data.frame(function_id = character(), function_group = character(), observed_groups = character(),
      annotation_status = character()))
  do.call(rbind, result)
}

# Load a reference from a source directory or the installed package.
.microfun_chord_reference <- function(name) {
  source_dir <- get0(".microfun_source_dir", envir = environment(), inherits = TRUE)
  if (!is.null(source_dir)) {
    path <- if (name == "ko_to_kegg_reference")
      file.path(source_dir, "R", "sysdata.rda") else file.path(source_dir, "inst", "extdata", paste0(name, ".RData"))
    e <- new.env(parent = emptyenv())
    load(path, envir = e)
    if (!exists(name, e, inherits = FALSE))
      stop("Missing reference object: ", name)
    return(get(name, e, inherits = FALSE))
  }
  if (name == "ko_to_kegg_reference")
    return(get(name, envir = asNamespace("MicroFun")))
  e <- new.env(parent = emptyenv())
  load(system.file("extdata", paste0(name, ".RData"), package = "MicroFun"), envir = e)
  get(name, e, inherits = FALSE)
}


# ---- Merged from tax2fun_chord_plot.R (implementation preserved) ----
.draw_tax2fun_chord <- function(x, start_degree = 180, side_gap = 10, sector_gap = 0.9, link_alpha = 0.32,
  bar_alpha = 0.95, effect_scale = 1, effect_cap = Inf, outer_bar_track_height = 0.135, label_track_height = 0.086,
  taxon_label_cex = 0.9, function_label_cex = 0.85, show_legend = TRUE, show_title = FALSE, title_text = NULL,
  direction_cols = c(Down = "#B4CBE0", Up = "#F5BB7C", Neutral = "#999999"), phylum_base_cols = microfun_palette("taxa"),
  function_group_palette = microfun_palette("function"), taxonomy_tint_step = 0, show_category_legend = show_legend,
  legend_cex = 0.72, show_outer_border = TRUE, outer_border_col = "#333333", outer_border_lwd = 0.8,
  narrow_labels = c("auto", "letters", "full"), label_codes = NULL, label_fit = 0.85, label_mapping_file = NULL,
  legend_min_cex = 0.65, legend_position = c("separate", "right"), taxonomy_label_cex = 0.7, narrow_sector_degrees = 0.5,
  prefer_full_taxa = TRUE, html_mode = FALSE) {
  if (length(narrow_sector_degrees) != 1L || !is.finite(narrow_sector_degrees) || narrow_sector_degrees <
    0)
    stop("narrow_sector_degrees must be non-negative.")
  legend_position <- match.arg(legend_position)
  narrow_labels <- match.arg(narrow_labels)
  if (length(label_fit) != 1L || !is.finite(label_fit) || label_fit <= 0 || label_fit > 1)
    stop("label_fit must be in (0,1].")
  if (!is.null(label_codes) && (is.null(names(label_codes)) || anyDuplicated(names(label_codes)) ||
    anyDuplicated(unname(label_codes)) || anyNA(label_codes) || any(!grepl("^[A-Z]+$", label_codes))))
    stop("label_codes must be a uniquely named vector of unique uppercase letter codes.")
  if (!all(c("Up", "Down", "Neutral") %in% names(direction_cols)))
    stop("direction_cols needs Up, Down and Neutral.")
  if (!is.finite(effect_scale) || effect_scale <= 0 || is.na(effect_cap) || effect_cap <= 0)
    stop("effect_scale and effect_cap must be positive.")
  if (!is.finite(side_gap) || side_gap < 0 || !is.finite(sector_gap) || sector_gap < 0)
    stop("Gaps must be finite and non-negative.")
  tx <- x$tax_table
  ft <- x$function_table
  links <- x$links
  ranks <- setdiff(names(tx), "taxon_id")
  tx <- tx[do.call(order, c(unname(tx[c(ranks, "taxon_id")]), list(na.last = TRUE))), , drop = FALSE]
  ft <- ft[order(ft$function_group, ft$function_id), , drop = FALSE]
  tax_ids <- paste0("tax_", seq_len(nrow(tx)))
  fun_ids <- paste0("fun_", seq_len(nrow(ft)))
  ids <- c(tax_ids, fun_ids)
  fnames <- if ("function_name" %in% names(ft)) as.character(ft$function_name) else ft$function_id
  missing_names <- is.na(fnames) | !nzchar(trimws(fnames))
  fnames[missing_names] <- ft$function_id[missing_names]
  flabels <- ifelse(fnames == ft$function_id, ft$function_id, paste0(fnames, "[", ft$function_id, "]"))
  labels <- stats::setNames(c(tx$taxon_id, flabels), ids)
  # Full lineage keys prevent homonymous classifications from sharing a code.
  lineage_key <- function(k, i) paste(c("Taxa", ranks[k], unlist(tx[i, ranks[seq_len(k)], drop = FALSE])),
    collapse = "|")
  terminal_key <- stats::setNames(c(paste0("Taxa|terminal|", tx$taxon_id), paste0("Function|terminal|",
    ft$function_id)), ids)
  keys <- c(unname(terminal_key), paste0("Function|group|", unique(ft$function_group)))
  if (length(ranks))
    for (k in seq_along(ranks)) keys <- c(keys, vapply(seq_len(nrow(tx)), function(i) lineage_key(k,
      i), character(1)))
  keys <- sort(unique(keys))
  codes <- label_codes
  if (is.null(codes))
    codes <- stats::setNames(character(), character())
  next_code <- 1L
  for (key in keys) if (!key %in% names(codes)) {
    repeat {
      candidate <- .microfun_letter_code(next_code)
      next_code <- next_code + 1L
      if (!candidate %in% unname(codes))
        break
    }
    codes[key] <- candidate
  }
  taxon_rank <- if (is.null(x$settings$taxon_level)) "Taxon" else x$settings$taxon_level
  function_rank <- if (is.null(x$settings$function_level)) "Function" else x$settings$function_level
  label_rows <- list()
  label_geometry <- list()
  choose_label <- function(key, label, side, rank, colour, info, cex, prefer_full = FALSE) {
    span <- info$start - info$end
    available <- abs(span) * pi/180 * mean(info$radii) * label_fit
    width <- graphics::strwidth(label, cex = cex, units = "user")
    short <- narrow_labels == "letters" ||
      (narrow_labels == "auto" && (width > available ||
        (prefer_full && abs(span) < narrow_sector_degrees)))
    code <- if (short) unname(codes[key]) else ""
    shown <- if (short) code else label
    visible <- if (short) graphics::strwidth(code, cex = cex, units = "user") <= available else
      width <= available
    if (!visible) shown <- ""
    # Keep the legacy leader column for callers consuming the label mapping.
    label_rows[[length(label_rows) + 1L]] <<- data.frame(
      key = key, side = side, rank = rank, full_name = label, code = code,
      display_label = shown, colour = unname(colour), abbreviated = short,
      leader = FALSE, visible = visible, stringsAsFactors = FALSE)
    if (html_mode) label_geometry[[length(label_geometry) + 1L]] <<- info
    shown
  }
  span_info <- function(sectors, track, last_track = track) {
    list(
      start = circlize::get.cell.meta.data("cell.start.degree", sector.index = sectors[1], track.index = track),
      end = circlize::get.cell.meta.data("cell.end.degree", sector.index = tail(sectors, 1), track.index = track),
      radii = c(
        circlize::get.cell.meta.data("cell.bottom.radius", sector.index = sectors[1], track.index = last_track),
        circlize::get.cell.meta.data("cell.top.radius", sector.index = sectors[1], track.index = track))
    )
  }
  links$t_id <- tax_ids[match(links$taxon, tx$taxon_id)]
  links$f_id <- fun_ids[match(links$function_id, ft$function_id)]
  links <- links[order(links$t_id, -links$contribution, links$f_id), , drop = FALSE]
  links$t_end <- ave(links$contribution, links$t_id, FUN = cumsum)
  links$t_start <- links$t_end - links$contribution
  links <- links[order(links$f_id, -links$contribution, links$t_id), , drop = FALSE]
  links$f_end <- ave(links$contribution, links$f_id, FUN = cumsum)
  links$f_start <- links$f_end - links$contribution
  widths <- c(stats::setNames(vapply(tax_ids, function(id) sum(links$contribution[links$t_id == id]),
    numeric(1)), tax_ids), stats::setNames(vapply(fun_ids, function(id) sum(links$contribution[links$f_id ==
    id]), numeric(1)), fun_ids))
  # Global rescaling preserves all width ratios and avoids near-zero sector ranges.
  width_scale <- max(widths)
  widths <- widths/width_scale
  for (col in c("t_start", "t_end", "f_start", "f_end")) links[[col]] <- links[[col]]/width_scale
  gap <- rep(sector_gap, length(ids))
  gap[c(length(tax_ids), length(ids))] <- side_gap
  if (sum(gap) >= 300) {
    if (2 * side_gap >= 300)
      stop("side_gap is too large.")
    gap[] <- min(sector_gap, (300 - 2 * side_gap)/max(1, length(ids) - 2))
    gap[c(length(tax_ids), length(ids))] <- side_gap
  }
  links$bar <- pmin(abs(links$effect) * effect_scale, effect_cap)
  ymax <- max(links$bar)
  if (ymax <= 0)
    ymax <- 1
  root <- if (length(ranks))
    tx[[ranks[1]]] else tx$taxon_id
  root_cols <- .microfun_category_colours(root, phylum_base_cols, "taxa")
  f_cols <- .microfun_category_colours(ft$function_group, function_group_palette, "function")
  if (length(taxonomy_tint_step) != 1L || !is.finite(taxonomy_tint_step) || taxonomy_tint_step < 0 ||
    taxonomy_tint_step > 1)
    stop("taxonomy_tint_step must be between 0 and 1.")
  mix_white <- function(col, amount) {
    rgb <- grDevices::col2rgb(col)/255
    rgb <- rgb * (1 - amount) + amount
    grDevices::rgb(rgb[1], rgb[2], rgb[3])
  }
  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  old_par <- graphics::par(c("mar", "xpd"))
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = rep(0.6, 4), xpd = FALSE)
  has_legend <- legend_position == "right" && (isTRUE(show_legend) || isTRUE(show_category_legend) ||
    narrow_labels != "full")
  # The browser gets fewer polygon vertices. Static exports retain native detail.
  if (html_mode) circlize::circos.par(unit.circle.segments = 240)
  circlize::circos.par(start.degree = start_degree, clock.wise = TRUE, gap.after = gap, track.margin = c(0.003,
    0.003), cell.padding = c(0, 0, 0, 0), canvas.xlim = if (has_legend)
    c(-1.4, 2.65) else c(-1.4, 1.4), canvas.ylim = c(-1.14, 1.14))
  circlize::circos.initialize(factors = factor(ids, levels = ids), xlim = cbind(rep(0, length(ids)),
    as.numeric(widths[ids])))
  circlize::circos.trackPlotRegion(ylim = c(0, ymax), track.height = outer_bar_track_height, bg.border = NA,
    panel.fun = function(xx, yy) {
      id <- circlize::get.cell.meta.data("sector.index")
      is_tax <- id %in% tax_ids
      z <- links[if (is_tax)
        links$t_id == id else links$f_id == id, , drop = FALSE]
      circlize::circos.lines(circlize::get.cell.meta.data("xlim"), c(0, 0), col = "grey40", lwd = 0.4)
      for (j in seq_len(nrow(z))) circlize::circos.rect(if (is_tax)
        z$t_start[j] else z$f_start[j], 0, if (is_tax)
        z$t_end[j] else z$f_end[j], z$bar[j], col = grDevices::adjustcolor(direction_cols[z$edge_direction[j]],
        alpha.f = bar_alpha), border = "white", lwd = 0.2)
    })
  # Group taxonomy tracks by lineage and share radial bands with function groups.
  n_tracks <- max(1L, length(ranks))
  tracks <- integer(n_tracks)
  for (k in seq_len(n_tracks)) {
    circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.04, bg.border = NA)
    tracks[k] <- circlize::get.current.track.index()
    if (length(ranks)) {
      lineage <- do.call(paste, c(tx[ranks[seq_len(k)]], sep = " / "))
      for (key in unique(lineage)) {
        idx <- which(lineage == key)
        colour <- mix_white(root_cols[root[idx[1]]], min(1, (k - 1) * taxonomy_tint_step))
        info <- span_info(tax_ids[idx], tracks[k])
        shown <- choose_label(lineage_key(k, idx[1]), tx[[ranks[k]]][idx[1]], "Taxa", ranks[k],
          colour, info, taxonomy_label_cex)
        circlize::highlight.sector(tax_ids[idx], track.index = tracks[k], col = colour, border = "white",
          text = if (nzchar(shown)) shown else NULL, cex = taxonomy_label_cex, niceFacing = TRUE)
      }
    } else {
      # Colour terminal taxa when ancestor ranks are absent.
      for (j in seq_along(tax_ids)) circlize::highlight.sector(tax_ids[j], track.index = tracks[k],
        col = root_cols[tx$taxon_id[j]], border = "white")
    }
  }
  for (fg in unique(ft$function_group)) {
    sectors <- fun_ids[ft$function_group == fg]
    info <- span_info(sectors, tracks[1], tail(tracks, 1))
    shown <- choose_label(paste0("Function|group|", fg), fg, "Function", "Group", f_cols[fg], info,
      taxonomy_label_cex)
    circlize::highlight.sector(sectors, track.index = tracks, col = f_cols[fg], border = "white",
      text = if (nzchar(shown)) shown else NULL, cex = taxonomy_label_cex, niceFacing = TRUE)
  }
  circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = label_track_height, bg.border = "#E6E6E6",
    bg.col = "white", panel.fun = function(xx, yy) {
      id <- circlize::get.cell.meta.data("sector.index")
      is_tax <- id %in% tax_ids
      info <- span_info(id, circlize::get.current.track.index())
      cex <- if (is_tax)
        taxon_label_cex else function_label_cex
      colour <- if (is_tax)
        root_cols[root[match(id, tax_ids)]] else f_cols[ft$function_group[match(id, fun_ids)]]
      shown <- choose_label(terminal_key[id], labels[id], if (is_tax)
        "Taxa" else "Function", if (is_tax) taxon_rank else function_rank, colour, info, cex,
        prefer_full = is_tax && prefer_full_taxa)
      if (nzchar(shown)) {
        circlize::circos.text(mean(circlize::get.cell.meta.data("xlim")), 0.5, shown,
          facing = "bending.inside", niceFacing = TRUE, cex = cex)
      }
    })
  links <- links[order(abs(links$effect), links$contribution), , drop = FALSE]
  for (j in seq_len(nrow(links))) circlize::circos.link(links$t_id[j], c(links$t_start[j], links$t_end[j]),
    links$f_id[j], c(links$f_start[j], links$f_end[j]), col = grDevices::adjustcolor(direction_cols[links$edge_direction[j]],
      alpha.f = link_alpha), border = NA)
  if (isTRUE(show_outer_border)) {
    radius <- circlize::get.cell.meta.data("cell.top.radius", sector.index = ids[1], track.index = 1)
    angle <- seq(0, 2 * pi, length.out = if (html_mode) 241L else 1441L)
    graphics::lines(radius * cos(angle), radius * sin(angle), col = outer_border_col, lwd = outer_border_lwd)
  }
  if (isTRUE(show_title) && !is.null(title_text))
    graphics::text(0, 1.1, title_text, cex = 0.9, font = 2)
  # Match legend swatches to taxonomy and function colours.
  blocks <- list()
  if (isTRUE(show_category_legend)) {
    blocks <- list(list(title = "Taxa", legend = names(root_cols), fill = unname(root_cols)), list(title = "Function",
      legend = names(f_cols), fill = unname(f_cols)))
  }
  label_map <- do.call(rbind, label_rows)
  if (html_mode) {
    x$label_geometry <- label_geometry
    # R user coordinates are cartesian. SVG coordinates point downwards.
    x$svg_axes <- list(
      x = graphics::grconvertX(c(0, 1), from = "user", to = "ndc"),
      y = 1 - graphics::grconvertY(c(0, 1), from = "user", to = "ndc"))
  }
  # Static legends explain visible abbreviations. The HTML panel lists all labels.
  aliases <- label_map[label_map$visible & label_map$abbreviated &
    !duplicated(label_map$key) & !html_mode, , drop = FALSE]
  for (side in c("Taxa", "Function")) {
    z <- aliases[aliases$side == side, , drop = FALSE]
    if (!nrow(z))
      next
    entry <- paste0(z$code, " - ", z$full_name)
    k <- which(vapply(blocks, function(b) identical(b$title, side), logical(1)))
    if (length(k)) {
      blocks[[k]]$legend <- c(blocks[[k]]$legend, entry)
      blocks[[k]]$fill <- c(blocks[[k]]$fill, z$colour)
    } else blocks <- c(blocks, list(list(title = side, legend = entry, fill = z$colour)))
  }
  if (length(blocks))
    blocks <- blocks[order(match(vapply(blocks, `[[`, character(1), "title"), c("Taxa", "Function")))]
  if (isTRUE(show_legend)) {
    ct <- x$settings$contrast
    width_text <- if (!is.null(x$settings$width_label))
      x$settings$width_label else paste("Width:", x$settings$contribution, "contribution")
    blocks <- c(blocks, list(list(title = "Effect direction", legend = c(paste(ct[2], ">", ct[1]),
      paste(ct[2], "<", ct[1]), "Zero effect"), fill = unname(direction_cols[c("Up", "Down", "Neutral")])),
      list(title = "Encoding", legend = c(width_text, paste("Outer bars: |", x$settings$effect_type,
        "|", sep = ""), paste0("Scale = ", effect_scale, "; cap = ", effect_cap), paste0("Bar maximum = ",
        signif(ymax, 3))))))
  }
  x$legend_blocks <- blocks
  x$legend_pages <- list()
  if (legend_position == "separate") {
    x$legend_pages <- blocks
  } else if (length(blocks)) {
    # Fit legends in the right panel, with Taxa above Function.
    for (k in seq_along(blocks)) blocks[[k]]$legend <- vapply(blocks[[k]]$legend, function(s) paste(strwrap(s,
      width = 35), collapse = "\n"), character(1))
    draw_block <- function(b, y, cex, plot) do.call(graphics::legend, c(list(x = 1.12, y = y, xjust = 0,
      yjust = 1, bty = "n", border = NA, cex = cex, title.adj = 0, x.intersp = 0.7, y.intersp = 1.05,
      plot = plot), b))
    sizes <- lapply(blocks, draw_block, y = 1.04, cex = legend_cex, plot = FALSE)
    total_h <- sum(vapply(sizes, function(z) z$rect$h, numeric(1))) + 0.045 * (length(blocks) - 1L)
    max_w <- max(vapply(sizes, function(z) z$rect$w, numeric(1)))
    fitted_cex <- legend_cex * min(1, 2.04/total_h, 1.08/max_w)
    if (fitted_cex < legend_min_cex) {
      # Keep readable typography: export the complete legend on following pages.
      x$legend_pages <- blocks
      graphics::text(1.12, 1, "Taxa / Function", adj = c(0, 1), cex = legend_cex, font = 2)
    } else {
      x$legend_pages <- list()
      y <- 1.04
      for (b in blocks) {
        z <- draw_block(b, y, fitted_cex, TRUE)
        y <- y - z$rect$h - 0.045
      }
    }
  }
  x$label_mapping <- label_map
  x$label_codes <- codes
  if (!is.null(label_mapping_file))
    utils::write.csv(label_map, label_mapping_file, row.names = FALSE)
  invisible(x)
}

.microfun_letter_code <- function(i) {
  out <- ""
  while (i > 0L) {
    i <- i - 1L
    out <- paste0(LETTERS[i%%26L + 1L], out)
    i <- i%/%26L
  }
  out
}

# Paginated readable legends, in Taxa then Function order; no silent truncation.
.microfun_chord_legend_pages <- function(blocks) {
  for (b in blocks) {
    labels <- vapply(b$legend, function(z) paste(strwrap(gsub("\n", " ", z), width = 65), collapse = "\n"),
      character(1))
    costs <- vapply(strsplit(labels, "\n", fixed = TRUE), length, integer(1)) + 1L
    groups <- list()
    current <- integer()
    used <- 0
    for (i in seq_along(labels)) {
      if (used + costs[i] > 64L && length(current)) {
        groups[[length(groups) + 1L]] <- current
        current <- integer()
        used <- 0
      }
      current <- c(current, i)
      used <- used + costs[i]
    }
    if (length(current))
      groups[[length(groups) + 1L]] <- current
    for (ii in groups) {
      graphics::par(mar = c(1, 1, 1, 1))
      graphics::plot.new()
      graphics::plot.window(c(0, 1), c(0, 1))
      args <- list(x = 0.03, y = 0.97, legend = labels[ii], title = b$title, bty = "n", cex = 1.1,
        xjust = 0, yjust = 1, x.intersp = 0.7, y.intersp = 1.1, title.adj = 0)
      if (!is.null(b$fill))
        args$fill <- b$fill[ii]
      do.call(graphics::legend, args)
    }
  }
}


# ---- Merged from chord_html.R (implementation preserved) ----
# Read installed assets, or the source tree when running source-based checks.
.microfun_chord_asset <- function(name) {
  source_dir <- get0(".microfun_source_dir", inherits = TRUE)
  path <- if (is.null(source_dir)) system.file("chord", name, package = "MicroFun") else
    file.path(source_dir, "inst", "chord", name)
  if (!nzchar(path) || !file.exists(path)) stop("Missing Chord viewer asset: ", name)
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

# Convert a circlize annular sector to an SVG path in the exported page space.
.microfun_chord_sector_path <- function(info, origin, scale) {
  point <- function(angle, radius) {
    radians <- angle * pi/180
    origin + scale * radius * c(cos(radians), sin(radians))
  }
  inner <- info$radii[1]
  outer <- info$radii[2]
  outer_start <- point(info$start, outer)
  outer_end <- point(info$end, outer)
  inner_end <- point(info$end, inner)
  inner_start <- point(info$start, inner)
  large <- as.integer(abs(info$start - info$end) > 180)
  radius <- abs(scale)
  list(
    path = sprintf(paste("M %.4f %.4f A %.4f %.4f 0 %d 1 %.4f %.4f",
                         "L %.4f %.4f A %.4f %.4f 0 %d 0 %.4f %.4f Z"),
                   outer_start[1], outer_start[2], outer * radius[1], outer * radius[2],
                   large, outer_end[1], outer_end[2],
                   inner_end[1], inner_end[2], inner * radius[1], inner * radius[2],
                   large, inner_start[1], inner_start[2]),
    center = unname(point(mean(c(info$start, info$end)), mean(info$radii)))
  )
}

.microfun_chord_html <- function(svg, rendered, width, height) {
  escape <- function(x) {
    x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    gsub('"', "&quot;", x, fixed = TRUE)
  }
  # The Cairo SVG uses points. Convert the captured R device fractions to
  # the actual viewBox, including its origin if it is nonzero.
  viewbox <- sub('(?s)^<svg[^>]*viewBox="([^"]+)"[^>]*>.*$', "\\1", svg, perl = TRUE)
  values <- suppressWarnings(as.numeric(strsplit(trimws(viewbox), "[,[:space:]]+")[[1]]))
  if (length(values) != 4L || any(!is.finite(values))) values <- c(0, 0, width * 72, height * 72)
  origin <- values[1:2] + values[3:4] * c(rendered$svg_axes$x[1], rendered$svg_axes$y[1])
  scale <- values[3:4] * c(diff(rendered$svg_axes$x), diff(rendered$svg_axes$y))
  keep <- !duplicated(rendered$label_mapping$key)
  labels <- rendered$label_mapping[keep, , drop = FALSE]
  sectors <- lapply(rendered$label_geometry[keep], .microfun_chord_sector_path,
                    origin = origin, scale = scale)
  blocks <- vapply(rendered$legend_blocks, function(b) {
    fill <- if (is.null(b$fill)) rep("#999999", length(b$legend)) else b$fill
    rows <- paste0('<li><span class="swatch" style="background:', escape(fill), '"></span>',
                   escape(b$legend), "</li>")
    paste0("<section><h2>", escape(b$title), "</h2><ul>", paste0(rows, collapse = ""), "</ul></section>")
  }, character(1))
  lists <- vapply(c("Taxa", "Function"), function(side) {
    i <- which(labels$side == side)
    names <- escape(labels$full_name[i])
    code <- ifelse(nzchar(labels$code[i]), paste0('<span class="code">', escape(labels$code[i]), "</span>"), "")
    rank <- escape(labels$rank[i])
    rows <- paste0(
      '<li class="sector-row"><button type="button" class="sector-entry" data-sector="', i,
      '" data-name="', names, '" data-rank="', rank, '" aria-pressed="false">', code,
      '<span><span class="sector-name">', names, '</span><small>', rank, "</small></span></button></li>"
    )
    paste0("<details open><summary>", side, " sectors (", length(i), ")</summary><ul>",
           paste0(rows, collapse = ""), "</ul></details>")
  }, character(1))
  paste0(
    '<!doctype html><html lang="en"><head><meta charset="UTF-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1"><title>MicroFun Chord</title><style>',
    .microfun_chord_asset("viewer.css"),
    '</style></head><body><main><div id="panel"><nav>',
    '<button id="plus" type="button" aria-label="Zoom in">+</button>',
    '<button id="minus" type="button" aria-label="Zoom out">&minus;</button>',
    '<button id="fit" type="button">Fit</button><span id="scale"></span></nav>',
    '<div id="chart">', svg, '</div></div><aside>',
    '<input id="search" type="search" aria-label="Search sectors" placeholder="Search names, ranks or codes">',
    '<p class="help">Select a sector below to highlight it on the chart.</p>',
    '<div id="selection" aria-live="polite"></div>',
    '<button id="clear" type="button" hidden>Clear selection</button>',
    paste0(lists, collapse = ""), paste0(blocks, collapse = ""),
    '</aside></main><script type="application/json" id="sector-geometry">',
    jsonlite::toJSON(sectors, auto_unbox = TRUE, digits = 8),
    '</script><script>', .microfun_chord_asset("viewer.js"), '</script></body></html>'
  )
}


# ---- Merged from chord_export.R (implementation preserved) ----
#' Export a chord diagram
#'
#' With no file specified, an HTML diagram is saved as microfun_chord.html.
#' Additional formats may be requested together with formats; all formats share
#' the same output basename. Rendering logic and the companion legend/mapping
#' outputs are unchanged.
#' @param x A microfun_chord object.
#' @param file Output path. Defaults to microfun_chord.html. When formats is
#'   supplied, the file basename is used for every requested format.
#' @param width,height Figure dimensions in inches.
#' @param formats Optional additional extensions: html, pdf, png, svg, rds.
#'   For example, formats = c("pdf", "png") saves HTML plus PDF and PNG
#'   when file is omitted. If file is explicitly set, its extension is
#'   always included.
#' @param ... Graphical arguments forwarded to the original chord renderer.
#' @return Invisibly, one normalized path for a single file or a named
#'   character vector of output paths for multiple formats.
#' @seealso tax2fun_chord
#' @export
save_chord <- function(x, file = "microfun_chord.html", width = 18, height = 16,
                       formats = NULL, ...) {
  if (!inherits(x, "microfun_chord"))
    stop("Expected a microfun_chord object.")
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file))
    stop("file must be a non-empty output filename.")
  supported <- c("html", "pdf", "png", "svg", "rds")
  primary <- tolower(tools::file_ext(file))
  if (!primary %in% supported)
    stop("Supported formats: html, pdf, png, svg, rds.")
  if (is.null(formats))
    return(.microfun_save_chord_one(x, file, width = width, height = height, ...))
  if (!is.character(formats) || anyNA(formats) ||
      any(!tolower(sub("^\\.", "", formats)) %in% supported))
    stop("formats must contain only: html, pdf, png, svg, rds.")
  formats <- unique(c(primary, tolower(sub("^\\.", "", formats))))
  base <- sub("\\.[^.]+$", "", file)
  paths <- vapply(formats, function(ext) {
    .microfun_save_chord_one(x, paste0(base, ".", ext), width = width,
                             height = height, ...)
  }, character(1))
  invisible(stats::setNames(paths, formats))
}

.microfun_save_chord_one <- function(x, file, width = 18, height = 16, ...) {
  if (!inherits(x, "microfun_chord"))
    stop("Expected a microfun_chord object.")
  if (any(!is.finite(c(width, height))) || length(width) != 1L || length(height) != 1L || min(width,
    height) <= 0)
    stop("width and height must be positive inches.")
  ext <- tolower(tools::file_ext(file))
  if (!ext %in% c("pdf", "png", "svg", "html", "rds"))
    stop("Supported formats: pdf, png, svg, html, rds.")
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (ext == "rds") {
    saveRDS(x, file)
    return(invisible(normalizePath(file)))
  }
  opts <- utils::modifyList(x$settings$graphical, list(...))
  opts <- utils::modifyList(list(legend_position = "separate", show_title = FALSE), opts)
  # HTML always places the legend in its own scrollable panel.
  if (ext == "html")
    opts$legend_position <- "separate"
  opts$html_mode <- ext == "html"
  temporary <- if (ext == "html")
    tempfile(fileext = ".svg") else file
  if (ext == "html")
    on.exit(unlink(temporary), add = TRUE)
  if (ext == "pdf") {
    .microfun_pdf_device(file, width = width, height = height)
  } else if (ext == "png") {
    grDevices::png(file, width = width * 180, height = height * 180, res = 180)
  } else {
    .microfun_svg_device(temporary, width = width, height = height)
  }
  dev <- grDevices::dev.cur()
  on.exit(if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE)
  rendered <- do.call(.draw_tax2fun_chord, c(list(x = x), opts))
  if (ext == "pdf" && length(rendered$legend_pages))
    .microfun_chord_legend_pages(rendered$legend_pages)
  grDevices::dev.off(dev)
  if (ext %in% c("png", "svg") && length(rendered$legend_pages)) {
    .microfun_pdf_device(sub("\\.[^.]+$", "_legend.pdf", file), width = 11, height = 11)
    tryCatch(.microfun_chord_legend_pages(rendered$legend_pages), finally = grDevices::dev.off())
  }
  utils::write.csv(rendered$label_mapping, sub("\\.[^.]+$", "_label_mapping.csv", file), row.names = FALSE)
  if (ext == "html") {
    if (!file.exists(temporary) || is.na(file.info(temporary)$size) ||
        file.info(temporary)$size == 0) {
      stop("Chord HTML export could not create a temporary SVG. ",
           "Install 'svglite' or check the graphics device on this system.")
    }
    svg <- paste(readLines(temporary, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    svg <- sub("(?s)^.*?(<svg)", "\\1", svg, perl = TRUE)
    svg <- .microfun_compact_chord_svg(svg)
    html <- .microfun_chord_html(svg, rendered, width, height)
    writeLines(enc2utf8(html), file, useBytes = TRUE)
  }
  invisible(normalizePath(file, mustWork = TRUE))
}

# Round geometry to 0.01 point; keep XML IDs and text untouched.
.microfun_compact_chord_svg <- function(svg) {
  attributes <- gregexpr("\\b(?:d|points|transform)=\"[^\"]*\"", svg, perl = TRUE)
  chunks <- regmatches(svg, attributes)[[1]]
  if (!length(chunks))
    return(svg)
  number <- "[-+]?(?:[0-9]*\\.[0-9]+|[0-9]+\\.?)(?:[eE][-+]?[0-9]+)?"
  chunks <- vapply(chunks, function(z) {
    indices <- gregexpr(number, z, perl = TRUE)
    tokens <- regmatches(z, indices)[[1]]
    if (length(tokens)) {
      values <- formatC(as.numeric(tokens), format = "f", digits = 2)
      values <- sub("0+$", "", values)
      values <- sub("[.]$", "", values)
      values[values == "-0"] <- "0"
      regmatches(z, indices) <- list(values)
    }
    z
  }, character(1))
  regmatches(svg, attributes) <- list(chunks)
  svg
}

# The native PDF device avoids Cairo device-closing failures on macOS.
.microfun_pdf_device <- function(file, width, height) {
  if (identical(Sys.info()[["sysname"]], "Darwin") ||
      !isTRUE(capabilities("cairo"))) {
    grDevices::pdf(file = file, width = width, height = height,
                   useDingbats = FALSE)
  } else {
    grDevices::cairo_pdf(filename = file, width = width, height = height,
                         family = "sans")
  }
}

# Prefer an SVG device that does not require Cairo; preserve the original
# device as fallback to keep behavior unchanged where it works.
.microfun_svg_device <- function(file, width, height) {
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(filename = file, width = width, height = height)
  } else if (isTRUE(capabilities("cairo"))) {
    grDevices::svg(filename = file, width = width, height = height,
                   onefile = FALSE)
  } else {
    stop("SVG/HTML chord export needs 'svglite' or a Cairo SVG device. ",
         "Install 'svglite' to enable this format.")
  }
}


# ---- Merged from fun2tax_complete.R ----
#' Chord diagram from function-oriented differential results
#'
#' Uses the same contribution weighting, filtering, taxonomy bands, labels and
#' exports as tax2fun_chord. Requires results from direction = 'function'.
#' @param pathway_sdaa_result Named list keyed by function ID.
#' @param stratified_path PICRUSt2 contribution TSV or TSV.GZ.
#' @param physeq A phyloseq object used for the differential analysis.
#' @param group Metadata grouping column.
#' @param ... Arguments passed to tax2fun_chord, except direction.
#' @return A microfun_chord object.
#' @export
fun2tax_chord <- function(pathway_sdaa_result, stratified_path, physeq, group, ...) {
  .microfun_result_direction(pathway_sdaa_result, "function")
  tax2fun_chord(pathway_sdaa_result, stratified_path, physeq, group,
    direction = "function", ...)
}
