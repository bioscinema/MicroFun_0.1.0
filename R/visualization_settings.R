# Shared defaults: category colours are separate from signed effect colours.
.microfun_direction_colours <- function() c(Up = "#F5BB7C", Down = "#B4CBE0", Neutral = "#999999")

.microfun_result_level <- function(result, taxon_level = NULL, fallback = "Genus") {
  recorded <- attr(result, "microfun")$taxon_level
  if (is.null(taxon_level)) taxon_level <- if (is.null(recorded)) fallback else recorded
  if (length(taxon_level) != 1L || is.na(taxon_level) || !nzchar(taxon_level))
    stop("taxon_level must be a non-empty rank name.")
  if (!is.null(recorded) && tolower(taxon_level) != tolower(recorded))
    stop("taxon_level differs from the settings recorded in pathway_sdaa_result.")
  standard <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
  hit <- match(tolower(taxon_level), tolower(standard))
  if (!is.na(hit)) taxon_level <- standard[hit]
  taxon_level
}

.microfun_result_direction <- function(result, expected) {
  recorded <- attr(result, "microfun")$direction
  if (!is.null(recorded) && !identical(recorded, expected))
    stop("Result direction is '", recorded, "'; this plot requires direction = '", expected, "'.")
  invisible(expected)
}

# Groups are not effect signs: reference/target order determines the two colours.
.microfun_group_colours <- function(x, result = NULL) {
  lev <- if (is.factor(x)) levels(droplevels(x)) else sort(unique(as.character(x[!is.na(x)])))
  info <- attr(result, "microfun")
  reference <- info$reference
  if (is.null(reference) && length(result)) {
    first <- result[[1L]]
    if (is.data.frame(first) && "group1" %in% names(first)) reference <- unique(as.character(first$group1))
  }
  if (length(reference) == 1L && reference %in% lev) lev <- c(reference, setdiff(lev, reference))
  cols <- c(unname(.microfun_direction_colours()[c("Down", "Up")]),
    setdiff(microfun_palette("taxa"), c("#B4CBE0", "#F5BB7C")))
  stats::setNames(rep(cols, length.out = length(lev)), lev)
}

# Prefer explicit names, then result annotations, then bundled references.
.microfun_function_names <- function(result, direction, level, ids) {
  records <- data.frame(id = character(), label = character())
  if (direction == "taxa") {
    for (z in result) {
      col <- intersect(c("pathway_name", "description", "function_name"), names(z))
      if (length(col) && "feature" %in% names(z))
        records <- rbind(records, data.frame(id = .chord_ids(z$feature, level), label = as.character(z[[col[1]]])))
    }
  }
  reference <- if (level == "KEGG") "kegg_reference" else paste0(level, "_reference")
  ref <- tryCatch(.microfun_chord_reference(reference), error = function(e) NULL)
  if (level == "KEGG") {
    ko <- tryCatch(.microfun_chord_reference("KO_reference"), error = function(e) NULL)
    if (is.data.frame(ko) && "Pathway" %in% names(ko)) {
      valid <- grepl("\\[PATH:ko[0-9]{5}\\]", ko$Pathway)
      records <- rbind(records, data.frame(
        id = sub(".*\\[PATH:(ko[0-9]{5})\\].*", "\\1", ko$Pathway[valid]),
        label = trimws(sub("\\s*\\[PATH:ko[0-9]{5}\\].*", "", ko$Pathway[valid]))))
    }
  }
  if (is.data.frame(ref)) {
    idcol <- intersect(c(level, "KO", "KEGG", "EC", "MetaCyc", "ID", "pathway", "pathway_id", "id"), names(ref))
    labelcol <- intersect(c("Description", "description", "pathway_name", "Name", "name", "PathwayName", "KoDescription"), names(ref))
    if (length(idcol) && length(labelcol))
      records <- rbind(records, data.frame(id = .chord_ids(ref[[idcol[1]]], level), label = as.character(ref[[labelcol[1]]])))
  }
  records <- records[!is.na(records$label) & nzchar(trimws(records$label)) & records$label != records$id, , drop = FALSE]
  records <- records[!duplicated(records$id), , drop = FALSE]
  out <- records$label[match(ids, records$id)]
  out[is.na(out)] <- ids[is.na(out)]
  stats::setNames(out, ids)
}

# The bundled EC and MetaCyc references use 'id', not their namespace name.
.microfun_reference_map <- function(pathway, ko_group = "PathwayL1") {
  ref <- .microfun_chord_reference(paste0(pathway, "_reference"))
  idcol <- intersect(c(pathway, "id", "ID"), names(ref))
  groupcol <- if (pathway == "KO") ko_group else "description"
  if (!length(idcol) || !groupcol %in% names(ref)) stop("Reference lacks ID/group columns for ", pathway)
  ids <- .chord_ids(ref[[idcol[1]]], pathway)
  labels <- as.character(ref[[groupcol]])
  valid <- !is.na(ids) & nzchar(ids) & !is.na(labels) & nzchar(labels)
  mapping <- stats::setNames(labels[valid], ids[valid])
  # Keep the original first packaged membership for multi-category KOs.
  mapping[!duplicated(names(mapping))]
}
