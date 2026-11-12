#' MicroFun categorical palettes
#'
#' @param side Either taxa (light muted nine colours) or function (dark muted
#'   nine colours). Direction colours belong to the plotting functions.
#' @return A character vector of nine hexadecimal colours.
#' @export
microfun_palette <- function(side = c("taxa", "function"), n = 9L) {
  side <- match.arg(side)
  base <- if (side == "taxa") {
    c("#B5C3D7", "#B6C9C0", "#F5D8B7", "#F1D0C6", "#DDAEAB", "#D0E0EF", "#E9CDDF", "#C8BFD9", "#E4E3BF")
  } else {
    c("#6E8FB2", "#7DA494", "#EAB67A", "#E5A79A", "#C16E71", "#ABC8E5", "#D8A0C1", "#9F8DB8", "#D0D08A")
  }
  .microfun_extend_colours(base, n, side)
}

.microfun_extend_colours <- function(palette, n, side = "taxa") {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 0 || n != floor(n))
    stop("n must be a non-negative integer.")
  rgb <- grDevices::col2rgb(palette)
  hex <- grDevices::rgb(rgb[1L, ], rgb[2L, ], rgb[3L, ], maxColorValue = 255)
  keep <- !duplicated(hex)
  palette <- unname(palette[keep])
  if (n <= length(palette))
    return(palette[seq_len(n)])

  to_luv <- function(colours) {
    grDevices::convertColor(t(grDevices::col2rgb(colours)) / 255, from = "sRGB", to = "Luv")
  }
  base_luv <- to_luv(hex[keep])
  hue <- (atan2(base_luv[, 3L], base_luv[, 2L]) * 180 / pi) %% 360
  lightness <- if (side == "taxa") seq(54, 89, by = 5) else seq(30, 75, by = 5)
  grid <- expand.grid(hue = hue, chroma = c(18, 30, 44, 58, 72), lightness = lightness)
  candidates <- unique(grDevices::hcl(grid$hue, grid$chroma, grid$lightness, fixup = FALSE))
  candidates <- setdiff(candidates[!is.na(candidates)], hex[keep])
  select_colours <- function(candidates, existing_luv, count) {
    if (count == 0L)
      return(character())
    candidate_luv <- to_luv(candidates)
    distance <- function(colour) rowSums(sweep(candidate_luv, 2L, colour, "-")^2)
    nearest <- rep(Inf, length(candidates))
    for (i in seq_len(nrow(existing_luv)))
      nearest <- pmin(nearest, distance(existing_luv[i, ]))
    selected <- character(count)
    for (i in seq_along(selected)) {
      best <- which.max(nearest)
      selected[i] <- candidates[best]
      nearest <- pmin(nearest, distance(candidate_luv[best, ]))
      nearest[best] <- -Inf
    }
    selected
  }
  count <- n - length(palette)
  added <- select_colours(candidates, base_luv, min(count, length(candidates)))
  if (length(added) == count)
    return(c(palette, added))

  grid <- expand.grid(hue = seq(0, 354, by = 6), chroma = seq(12, 84, by = 6),
    lightness = seq(min(lightness), max(lightness), by = 2))
  extra <- grDevices::hcl(grid$hue, grid$chroma, grid$lightness, fixup = FALSE)
  extra <- setdiff(unique(extra[!is.na(extra)]), c(hex[keep], added))
  remaining <- count - length(added)
  if (length(extra) < remaining)
    stop("Too many categories for the available colour range.")
  c(palette, added, select_colours(extra, to_luv(c(hex[keep], added)), remaining))
}

# Named palettes pin category colours across independently subsetted figures.
.microfun_category_colours <- function(labels, palette, side = "taxa") {
  labels <- sort(unique(as.character(labels)))
  if (!length(labels))
    return(stats::setNames(character(), character()))
  if (anyNA(labels) || any(!nzchar(labels)))
    stop("Colour categories must be non-empty.")
  if (!is.character(palette) || !length(palette) || anyNA(palette))
    stop("A palette must contain valid colour strings.")
  grDevices::col2rgb(palette)
  if (!is.null(names(palette))) {
    if (anyNA(names(palette)) || any(!nzchar(names(palette))) || anyDuplicated(names(palette)))
      stop("A named palette must have unique, non-empty category names.")
    missing <- setdiff(labels, names(palette))
    if (length(missing))
      stop("Named palette lacks categories: ", paste(missing, collapse = ", "))
    return(palette[labels])
  }
  # Reserve fixed slots for the six common KEGG level-1 categories.
  canonical <- c("Metabolism", "Genetic Information Processing", "Environmental Information Processing",
    "Cellular Processes", "Organismal Systems", "Human Diseases")
  slots <- match(tolower(labels), tolower(canonical))
  if (side == "function" && any(!is.na(slots))) {
    extra <- which(is.na(slots))
    slots[extra] <- length(canonical) + seq_along(extra)
    n_needed <- max(slots)
  } else {
    slots <- seq_along(labels)
    n_needed <- length(labels)
  }
  stats::setNames(.microfun_extend_colours(palette, n_needed, side)[slots], labels)
}

# Colour identities are separate from displayed names and edge directions.
.microfun_sankey_style <- function(edges, taxa_on_left, taxon_palette, function_palette, taxon_groups,
  node_color, node_neutral = NULL) {
  if (!all(c("Up", "Down") %in% names(node_color)))
    stop("node_color needs Up and Down.")
  grDevices::col2rgb(node_color[c("Up", "Down")])
  src <- unique(as.character(edges$source))
  dst <- unique(as.character(edges$target))
  nodes <- data.frame(name = c(src, dst), side = c(rep("left", length(src)), rep("right", length(dst))),
    stringsAsFactors = FALSE)
  nodes$kind <- ifelse(nodes$side == if (taxa_on_left)
    "left" else "right", "taxa", "function")
  taxa <- nodes$name[nodes$kind == "taxa"]
  if (is.null(taxon_groups))
    taxon_groups <- stats::setNames(taxa, taxa)
  if (is.null(names(taxon_groups)) || anyDuplicated(names(taxon_groups)) || any(!taxa %in% names(taxon_groups)))
    stop("taxon_groups must uniquely name every plotted taxon.")
  groups <- as.character(taxon_groups[taxa])
  # Use the complete supplied mapping as the colour universe, even after filtering.
  tcols <- .microfun_category_colours(unname(taxon_groups), taxon_palette, "taxa")
  fcols <- .microfun_category_colours(nodes$name[nodes$kind == "function"], function_palette, "function")
  nodes$category <- nodes$name
  nodes$category[nodes$kind == "taxa"] <- groups
  nodes$colour <- NA_character_
  nodes$colour[nodes$kind == "taxa"] <- unname(tcols[groups])
  nodes$colour[nodes$kind == "function"] <- unname(fcols[nodes$name[nodes$kind == "function"]])
  if (!is.null(node_neutral)) {
    if (length(node_neutral) != 1L || is.na(node_neutral))
      stop("node_neutral must be one colour or NULL.")
    grDevices::col2rgb(node_neutral)
    nodes$colour[] <- node_neutral
  }
  nodes$node_group <- paste0("node_", seq_len(nrow(nodes)))
  edges$source_id <- match(edges$source, src) - 1L
  edges$target_id <- length(src) + match(edges$target, dst) - 1L
  domain <- c(nodes$node_group, "Up", "Down")
  colours <- c(nodes$colour, unname(node_color[c("Up", "Down")]))
  scale <- htmlwidgets::JS(sprintf("d3.scaleOrdinal().domain(%s).range(%s)", jsonlite::toJSON(domain,
    auto_unbox = FALSE), jsonlite::toJSON(colours, auto_unbox = FALSE)))
  list(nodes = nodes, edges = edges, scale = scale)
}
