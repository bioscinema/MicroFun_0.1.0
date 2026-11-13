#' Boxplots of taxon contributions within functions or function groups
#'
#' Aggregate stratified contributions for the function IDs in a function-oriented result list and compare their taxonomic contributors across sample groups.
#'
#' @param stratified_path Character string. Path to a stratified PICRUSt2 output file
#'   (typically `*_stratified.tsv`), containing columns for function, taxon, sample, and contribution values.
#'
#' @param physeq A \code{\link[phyloseq]{phyloseq}} object containing sample metadata
#'   and taxonomy table. Sample IDs must match those in \code{stratified_path}.
#'
#' @param group Character scalar. Column name in sample metadata to use for grouping (e.g. treatment).
#'
#' @param taxon_level Taxonomy rank inferred from the result metadata; defaults to
#' Genus for legacy lists without metadata. Conflicting explicit ranks are rejected.
#'
#' @param pathway Character scalar. Functional annotation scheme used: one of
#'   \code{"KO"}, \code{"EC"}, or \code{"MetaCyc"}. Defaults to "KO".
#'
#' @param value_col Character scalar. Name of the numeric contribution column
#'   in the stratified table. Defaults to the column recorded in the SDAA
#'   result or the PICRUSt2 abundance column when available.
#'
#' @param wrap_width Integer. Width (in characters) for line-wrapping taxon labels
#'   on the x-axis. Defaults to 36.
#'
#' @importFrom data.table %chin% :=
#'
#' @import ggplot2
#'
#' @importFrom stringr str_wrap
#'
#' @importFrom phyloseq sample_data tax_table
#'
#' @export
#'
#' @param pathway_sdaa_result Named result list from `pathway_sdaa(direction = "function")`, keyed by function ID. List names select functions; rows within each result table do not restrict the plotted taxa.
#'
#' @param group_functions Logical. Group KO IDs by the bundled `PathwayL2` annotation, or EC/MetaCyc IDs by description. Unmatched functions are pooled into `Unmapped`. If `FALSE`, return a plot for each function ID.
#'
#' @param top_n Positive integer or `Inf` (default, retain all). Within each plot, rank features by the mean of their group mean contributions, giving equal weight to each group. Ties are ordered by feature ID. Missing group assignments are excluded when selecting a finite number of features.
#'
#' @return A named list of ggplot objects, one per function or reference group with matching contributions. Empty plots are omitted.
#'
#' @details ASVs are aggregated to `taxon_level`. Missing taxon/sample combinations are filled with zero. The y-axis is `log10(pmax(contribution, 0) + 1e-6)`. Input contribution units are retained before this transformation. This function neither reruns differential analysis nor filters taxa by adjusted p-value.
#'
#' @seealso [pathway_sdaa()], [tax2fun_box()]
#'
#' @md
fun2tax_box <- function(
    stratified_path,
    pathway_sdaa_result,
    physeq,
    group,
    taxon_level = NULL,
    pathway     = c("KO","EC","MetaCyc"),
    value_col   = NULL,
    wrap_width  = 36,
    group_functions = TRUE,
    top_n = Inf
){
  .microfun_check_box_top_n(top_n)
  taxon_level <- .microfun_result_level(pathway_sdaa_result, taxon_level)
  .microfun_result_direction(pathway_sdaa_result, "function")
  # Match the functional namespace.
  pathway <- match.arg(pathway)
  eps <- 1e-6

  # Normalize identifiers.
  normalize_fun_ids <- function(x, pathway){
    x <- as.character(x)
    if (pathway == "KO") {
      toupper(gsub("^ko:", "", x, ignore.case = TRUE))
    } else if (pathway == "EC") {
      x <- gsub("^ec:", "", x, ignore.case = TRUE)
      toupper(gsub("^EC[\\s_]*", "", x, ignore.case = TRUE))
    } else { # MetaCyc
      toupper(gsub("^metacyc:", "", x, ignore.case = TRUE))
    }
  }
  is_bad_lvl <- function(z){
    z <- as.character(z)
    is.na(z) | !nzchar(z) | z %in% c("Unassigned","Unclassified","Unknown")
  }

  # Read stratified contributions.
  df <- if (grepl("\\.gz$", stratified_path, ignore.case = TRUE) && !requireNamespace("R.utils", quietly = TRUE)) {
    con <- gzfile(stratified_path, "rt")
    on.exit(close(con), add = TRUE)
    data.table::as.data.table(utils::read.delim(con, check.names = FALSE))
  } else data.table::fread(stratified_path)
  data.table::setDT(df)
  data.table::setnames(df, tolower(names(df)))
  if (is.null(value_col)) {
    info <- attr(pathway_sdaa_result, "microfun")
    value_col <- .microfun_contribution_column(names(df), info$value_col)
  }
  # ensure required columns exist (case-insensitive for value_col)
  req <- c("function","taxon","sample")
  if (!all(req %in% names(df))) {
    stop("Stratified file must include columns: function, taxon, sample, and ", value_col)
  }
  if (!value_col %in% names(df)) {
    if (tolower(value_col) %in% names(df)) value_col <- tolower(value_col)
  }
  if (!value_col %in% names(df)) {
    stop("Column '", value_col, "' not found in stratified file.")
  }

  # Extract sample metadata and taxonomy.
  meta_df <- data.frame(phyloseq::sample_data(physeq), check.names = FALSE)
  meta_df$sample <- rownames(meta_df)
  meta <- data.table::as.data.table(meta_df)
  data.table::setDT(meta)
  if (!group %in% names(meta)) {
    stop("Grouping variable '", group, "' not found in sample_data(physeq).")
  }
  sample_ids <- as.character(meta$sample)

  tax_tab <- as.data.frame(as(phyloseq::tax_table(physeq), "matrix"), stringsAsFactors = FALSE)
  tax_tab$ASV <- rownames(tax_tab)
  taxtab <- data.table::as.data.table(tax_tab)

  # Map functions to reference categories.
  load_reference_map <- function(pathway) .microfun_reference_map(pathway, "PathwayL2")

  id2group <- NULL
  if (isTRUE(group_functions)) {
    id2group <- load_reference_map(pathway)
  }

  # Assign function group labels.
  fid_vec <- names(pathway_sdaa_result)
  if (is.null(fid_vec) || !length(fid_vec)) {
    stop("pathway_sdaa_result must be a named list keyed by function IDs.")
  }

  grp_vec <- if (is.null(id2group)) {
    fid_vec
  } else {
    unname(id2group[ normalize_fun_ids(fid_vec, pathway) ])
  }
  grp_vec[is.na(grp_vec) | !nzchar(grp_vec)] <- "Unmapped"

  # Normalize function IDs.
  df[, fid_norm := normalize_fun_ids(`function`, pathway)]

  # Join taxonomy at the requested rank.
  taxtab[, level := .SD[[1]], .SDcols = taxon_level]
  data.table::setkey(taxtab, ASV)
  # ensure types
  df[, sample := as.character(sample)]
  # restrict to samples in physeq early
  df <- df[sample %chin% sample_ids]

  # helpful keys
  data.table::setkey(df, taxon)           # df$taxon matches taxtab$ASV
  # bring level into df, drop bad levels
  df_annot <- taxtab[df, on = .(ASV = taxon), nomatch = 0L][
    !is_bad_lvl(level)
  ]

  # Select the contribution column.
  df_annot[, value := .SD[[1]], .SDcols = value_col]

  # Sum contributions by function, taxon and sample.
  df2 <- df_annot[
    , .(value = sum(value, na.rm = TRUE)),
    by = .(fid_norm, level, sample)
  ]
  # indices for fast lookups
  data.table::setindexv(df2, c("fid_norm","level","sample"))
  data.table::setindexv(meta, "sample")

  # Group normalized function identifiers.
  fid_vec_norm <- normalize_fun_ids(fid_vec, pathway)
  names(fid_vec_norm) <- fid_vec
  groups_norm <- split(fid_vec_norm, grp_vec)

  # Sum contributions across functions.
  build_long_for_fids_fast <- function(fids_norm) {
    if (!length(fids_norm)) return(NULL)
    dt <- df2[fid_norm %chin% fids_norm,
              .(abundance = sum(value, na.rm = TRUE)),
              by = .(level, sample)]
    if (!nrow(dt)) return(NULL)
    dt[]
  }

  # Plot contributions by group.
  plots <- vector("list", length(groups_norm))
  names(plots) <- names(groups_norm)

  for (glab in names(groups_norm)) {
    fids_n <- groups_norm[[glab]]
    long_sum <- build_long_for_fids_fast(fids_n)
    if (is.null(long_sum) || !nrow(long_sum)) next

    sample_grid <- data.table::CJ(level = unique(long_sum$level), sample = sample_ids, unique = TRUE)
    long_sum <- merge(sample_grid, long_sum, by = c("level", "sample"), all.x = TRUE)
    long_sum[is.na(abundance), abundance := 0]
    # attach metadata and keep rows that have the grouping variable
    long_sum <- meta[long_sum, on = "sample", nomatch = 0L]
    if (!group %in% names(long_sum)) next
    if (!nrow(long_sum)) next
    if (is.finite(top_n)) {
      top_levels <- .microfun_box_top_levels(long_sum, "level", group, top_n)
      if (!length(top_levels)) next
      keep_rows <- long_sum$level %chin% top_levels & !is.na(long_sum[[group]])
      long_sum <- long_sum[keep_rows, ]
      long_sum[, level := factor(level, levels = top_levels)]
    }

    # cosmetics
    long_sum[, abundance_log := log10(pmax(abundance, 0) + eps)]
    long_sum[[group]] <- as.factor(long_sum[[group]])

    p <- ggplot(
      long_sum,
      aes(x = level, y = abundance_log, fill = .data[[group]])
    ) +
      geom_boxplot(outlier.shape = 16, outlier.size = 0.8, alpha = 1) +
      scale_fill_manual(values = .microfun_group_colours(meta[[group]], pathway_sdaa_result)) +
      scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = wrap_width)) +
      labs(
        title = glab,
        x = taxon_level,
        y = "Log10 contribution",
        fill = group
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title.position = "plot"
      )

    plots[[glab]] <- p
  }

  Filter(Negate(is.null), plots)
}
