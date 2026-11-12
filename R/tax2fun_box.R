#' Plot Functional Abundance Boxplots for Driving Taxa
#'
#' This function generates boxplots of stratified functional pathway abundances
#' (from PICRUSt2 output) at a specified taxonomic level for each driving taxon,
#' grouped by a metadata variable.
#'
#' @param stratified_path Path to the stratified output file from PICRUSt2 (typically `*_strat.tsv`).
#' @param pathway_sdaa_result A named list of data.frames containing differential results
#'   (with `feature` and optionally `pathway_name`) for each taxon.
#' @param physeq A `phyloseq` object containing the microbial abundance table, taxonomy, and sample metadata.
#' @param taxon_level A string indicating the taxonomic level to aggregate functions by
#'   (e.g., `"Genus"`, `"Family"`). Inferred from recorded analysis settings; otherwise `"Genus"`.
#' @param group A string specifying the sample metadata column used for grouping in boxplots.
#' @param kegg_map Optional PICRUSt2 mapping TSV or named member list. Defaults
#'   to the mapping recorded in the SDAA result or the package reference.
#'
#' @return A named list of `ggplot2` boxplot objects, one for each taxon with matching functional results.
#'
#' @details This function:
#' \itemize{
#'   \item Parses PICRUSt2 stratified output and maps ASVs to their taxonomic levels.
#'   \item Aggregates function abundances at the taxon level.
#'   \item Filters and matches functions found in the user-provided differential result per taxon.
#'   \item Generates and returns a list of ggplot boxplots of function abundance by group.
#' }
#'
#' @import ggplot2 phyloseq
#' @export
#'
#' @examples
#' \dontrun{
#' plot_list <- tax2fun_box(
#'   stratified_path = "ko_metagenome_out/pred_metagenome_contrib.tsv",
#'   pathway_sdaa_result = diff_results,
#'   physeq = my_phyloseq,
#'   taxon_level = NULL,
#'   group = "disease_status"
#' )
#' plot_list$Bacteroides  # View boxplot for Bacteroides genus
#' }
tax2fun_box <- function(stratified_path, pathway_sdaa_result, physeq, taxon_level = NULL, group, kegg_map = NULL) {

  taxon_level <- .microfun_result_level(pathway_sdaa_result, taxon_level)
  .microfun_result_direction(pathway_sdaa_result, "taxa")
  # Load stratified Picrust2 output
  df <- if (grepl("\\.gz$", stratified_path, ignore.case = TRUE) && !requireNamespace("R.utils", quietly = TRUE)) {
    con <- gzfile(stratified_path, "rt")
    on.exit(close(con), add = TRUE)
    data.table::as.data.table(utils::read.delim(con, check.names = FALSE))
  } else fread(stratified_path, sep = "\t")
  info <- attr(pathway_sdaa_result, "microfun")
  feature_ids <- unlist(lapply(pathway_sdaa_result, function(x) as.character(x$feature)), use.names = FALSE)
  function_level <- if (!is.null(info$function_level)) info$function_level else
    .microfun_function_level(feature_ids)
  map_ko <- function_level == "KEGG" && .microfun_kegg_uses_ko(df[["function"]], info$pathway)
  value_col <- .microfun_contribution_column(names(df), info$value_col, require_abundance = function_level == "KEGG")
  if (map_ko) {
    if (is.null(kegg_map)) kegg_map <- info$kegg_map
    map <- .microfun_kegg_links(.microfun_kegg_members(kegg_map), .microfun_kegg_ids(feature_ids))
  }

  # Merge taxonomy to desired level
  tax_table_df <- as.data.frame(tax_table(physeq))
  tax_table_df$ASV <- rownames(tax_table_df)
  tax_table_df <- tax_table_df[, c("ASV", taxon_level)]

  df_split <- split(df, df$taxon)
  result_list <- lapply(df_split, function(dt) {
    data.table::dcast(dt, `function` ~ sample, value.var = value_col, fun.aggregate = sum, fill = 0)
  })

  sample_names <- sample_names(physeq)
  result_list_subset <- lapply(result_list, function(dt) {
    cols_to_keep <- c("function", intersect(names(dt), sample_names))
    dt[, ..cols_to_keep]
  })

  result_list_merged <- list()
  for (i in seq_along(result_list_subset)) {
    asv   <- names(result_list_subset)[i]
    level <- tax_table_df[tax_table_df$ASV == asv, taxon_level]
    if (length(level) == 0 || is.na(level) || level == "" ) {
      next
    }
    if (tolower(level) == "unknown" || tolower(level) == "uncultured") next

    dt <- data.table::as.data.table(result_list_subset[[i]])
    if (!"function" %in% names(dt) || ncol(dt) <= 1L) next

    # long format for this ASV/taxon
    dt_long <- reshape2::melt(
      dt,
      id.vars       = "function",
      variable.name = "sample",
      value.name    = "value"
    )

    if (level %in% names(result_list_merged)) {
      # existing table for this level → melt too
      existing      <- data.table::as.data.table(result_list_merged[[level]])
      existing_long <- reshape2::melt(
        existing,
        id.vars       = "function",
        variable.name = "sample",
        value.name    = "value"
      )

      # combine and sum duplicates (function × sample)
      combined <- data.table::rbindlist(list(existing_long, dt_long), use.names = TRUE, fill = TRUE)
      combined <- combined[, .(value = sum(value, na.rm = TRUE)), by = .(`function`, sample)]

      # back to wide: rows = function, cols = sample
      result_list_merged[[level]] <- data.table::dcast(
        combined,
        `function` ~ sample,
        value.var = "value",
        fill = 0
      )
    } else {
      result_list_merged[[level]] <- dt  # keep as data.table
    }
  }

  meta <- data.frame(sample_data(physeq))
  meta$sample <- rownames(meta)
  meta <- meta[meta$sample %in% df$sample, , drop = FALSE]

  plot_list <- list()
  genera <- names(pathway_sdaa_result)

  for (genus in genera) {
    genus_df <- result_list_merged[[genus]]
    genus_diff <- pathway_sdaa_result[[genus]]

    if (is.null(genus_df) || !nrow(genus_diff)) next

    # Match full feature IDs, mapping member KOs before summing pathway contributions.
    genus_df_long <- data.table::melt(genus_df, id.vars = "function", variable.name = "sample", value.name = "abundance")
    if (map_ko) {
      genus_df_long$ko_id <- .microfun_ko_ids(genus_df_long[["function"]])
      genus_df_long <- merge(genus_df_long, map, by = "ko_id", allow.cartesian = TRUE)
      genus_df_long$feature <- genus_df_long$pathway_id
    } else {
      genus_df_long$feature <- .chord_ids(genus_df_long[["function"]], function_level)
    }
    genus_diff <- as.data.frame(genus_diff)
    genus_diff$feature <- .chord_ids(genus_diff$feature, function_level)
    genus_df_long <- genus_df_long[feature %in% genus_diff$feature,
      .(abundance = sum(abundance)), by = .(feature, sample)]
    if (!nrow(genus_df_long)) next

    labels <- if ("pathway_name" %in% names(genus_diff)) genus_diff$pathway_name else genus_diff$description
    if (is.null(labels)) labels <- genus_diff$feature
    labels <- as.character(labels)
    missing_labels <- is.na(labels) | !nzchar(labels)
    labels[missing_labels] <- genus_diff$feature[missing_labels]
    labels <- setNames(labels, genus_diff$feature)
    sample_grid <- expand.grid(feature = unique(genus_df_long$feature), sample = meta$sample, stringsAsFactors = FALSE)
    genus_df_long <- merge(sample_grid, genus_df_long, by = c("feature", "sample"), all.x = TRUE)
    genus_df_long$abundance[is.na(genus_df_long$abundance)] <- 0
    genus_df_long <- merge(genus_df_long, meta, by = "sample")
    genus_df_long$pathway_name <- unname(labels[genus_df_long$feature])

    p <- ggplot(genus_df_long, aes(x = feature, y = abundance, fill = .data[[group]])) +
      geom_boxplot(alpha = 1) +
      scale_fill_manual(values = .microfun_group_colours(meta[[group]], pathway_sdaa_result)) +
      scale_x_discrete(labels = labels) +
      labs(title = genus, x = "Function", y = "Contribution",fill=group) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    plot_list[[genus]] <- p
  }

  return(plot_list)
}

