#' Convert KO abundance in picrust2 export files to KEGG pathway abundance
#'
#' This function takes a file containing KO (KEGG Orthology) abundance data in picrust2 export format and converts it to KEGG pathway abundance data.
#' The input file should be in .tsv, .txt, or .csv format.
#'
#' @param file A character string representing the file path of the input file containing KO abundance data in picrust2 export format. The input file should have KO identifiers in the first column and sample identifiers in the first row. The remaining cells should contain the abundance values for each KO-sample pair.
#' @param data An optional data.frame containing KO abundance data in the same format as the input file. If provided, the function will use this data instead of reading from the file. By default, this parameter is set to NULL.
#' @param kegg_map Optional PICRUSt2 mapping TSV (one pathway and its member KOs
#'   per line, without a header), or a named list of member KOs by pathway.
#'   Defaults to the package reference.
#'
#' @return
#' A data frame with KEGG pathway abundance values. Rows represent KEGG pathways, identified by their KEGG pathway IDs. Columns represent samples, identified by their sample IDs from the input file. Each cell contains the abundance of a specific KEGG pathway in a given sample, calculated by summing the abundances of the corresponding KOs in the input file.
#' @examples
#' \dontrun{
#' library(readr)
#'
#' # Example 1: Demonstration with a hypothetical input file
#'
#' # Prepare an input file path
#' input_file <- "path/to/your/picrust2/results/pred_metagenome_unstrat.tsv"
#'
#' # Run ko2kegg_abundance function
#' kegg_abundance <- ko2kegg_abundance(file = input_file)
#'
#' # Alternatively, read the data from a file and use the data argument
#' file_path <- "path/to/your/picrust2/results/pred_metagenome_unstrat.tsv"
#' ko_abundance <- read_delim(file_path, delim = "\t")
#' kegg_abundance <- ko2kegg_abundance(data = ko_abundance)
#' }
#' @export
ko2kegg_abundance <- function (file = NULL, data = NULL, kegg_map = NULL) {
  # Validate basic arguments.
  if (is.null(file) & is.null(data)) {
    stop("Error: Please provide either a file or a data.frame.")
  }

  if (!is.null(file) && !is.null(data)) {
    warning("Both file and data provided. Using data and ignoring file.")
  }

  # Validate the file format.
  validate_file_format <- function(file_path) {
    valid_extensions <- c(".txt", ".tsv", ".csv")
    ext <- tolower(tools::file_ext(file_path))
    if (!paste0(".", ext) %in% valid_extensions) {
      stop(sprintf("Error: Input file should be in %s format.",
                  paste(valid_extensions, collapse = ", ")))
    }
  }

  # Validate the data frame.
  validate_abundance_data <- function(df) {
    if (ncol(df) < 2) {
      stop("Error: Data must contain at least one KO column and one sample column.")
    }
    if (!is.character(df[[1]]) && !is.factor(df[[1]])) {
      stop("Error: First column must contain KO identifiers.")
    }
    if (!all(sapply(df[-1], is.numeric))) {
      stop("Error: All sample columns must contain numeric values.")
    }
  }

  # Define input validation.
  validate_input <- function(data) {
    # Validate the basic structure.
    if (!is.data.frame(data)) {
      stop("The provided data must be a data.frame")
    }

    if (ncol(data) < 2) {
      stop("Data must contain at least one KO column and one sample column")
    }

    # Check the first column (KO IDs).
    ko_col <- .microfun_ko_ids(data[[1]])

    # Check the KO ID format.
    ko_pattern <- "^K\\d{5}$"
    invalid_kos <- ko_col[!grepl(ko_pattern, ko_col)]
    if (length(invalid_kos) > 0) {
      warning(sprintf(
        "Found %d KO IDs that don't match the expected format (K#####).\nFirst few invalid IDs: %s",
        length(invalid_kos),
        paste(head(invalid_kos, 3), collapse = ", ")
      ))
    }

    # Check numeric columns.
    numeric_cols <- data[,-1, drop = FALSE]

    # Require numeric values in all sample columns.
    non_numeric_cols <- names(numeric_cols)[!sapply(numeric_cols, is.numeric)]
    if (length(non_numeric_cols) > 0) {
      stop(sprintf(
        "The following columns contain non-numeric values: %s",
        paste(non_numeric_cols, collapse = ", ")
      ))
    }

    # Check for negative values.
    neg_values <- sapply(numeric_cols, function(x) any(x < 0, na.rm = TRUE))
    if (any(neg_values)) {
      stop(sprintf(
        "Negative abundance values found in columns: %s",
        paste(names(numeric_cols)[neg_values], collapse = ", ")
      ))
    }

    # Check for missing values.
    na_counts <- sapply(numeric_cols, function(x) sum(is.na(x)))
    cols_with_na <- names(na_counts[na_counts > 0])
    if (length(cols_with_na) > 0) {
      warning(sprintf(
        "Missing values found in columns: %s\nTotal NA count per column: %s",
        paste(cols_with_na, collapse = ", "),
        paste(paste(cols_with_na, na_counts[cols_with_na], sep = ": "), collapse = ", ")
      ))
    }

    # Check for all-zero rows.
    zero_rows <- rowSums(numeric_cols == 0) == ncol(numeric_cols)
    if (any(zero_rows)) {
      warning(sprintf(
        "%d KOs have zero abundance across all samples",
        sum(zero_rows)
      ))
    }

    # Check that column names are unique.
    if (any(duplicated(colnames(data)))) {
      stop("Duplicate column names found in the input data")
    }

    return(TRUE)
  }

  # Validate the file format.
  if (!is.null(file)) {
    validate_file_format(file)
    message("Loading data from file...")
    file_format <- substr(file, nchar(file) - 3, nchar(file))
    abundance <- switch(
      file_format,
      ".txt" = readr::read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
      ".tsv" = readr::read_delim(file, delim = "\t", escape_double = FALSE, trim_ws = TRUE),
      ".csv" = readr::read_delim(file, delim = ",", escape_double = FALSE, trim_ws = TRUE)
    )
  } else {
    message("Processing provided data frame...")
    abundance <- data
  }

  # Run input validation.
  tryCatch({
    validate_input(abundance)
  }, warning = function(w) {
    warning("Warning during data validation: ", w$message, call. = FALSE)
  }, error = function(e) {
    stop("Data validation failed: ", e$message, call. = FALSE)
  })

  message("Loading KEGG reference data. This might take a while...")
  members <- .microfun_kegg_members(kegg_map)
  kegg_names <- names(members)
  ko_ids <- .microfun_ko_ids(abundance[[1]])

  message(sprintf("Processing %d KEGG pathways...", length(unique(kegg_names))))

  # Create the result data frame.
  kegg_abundance <- data.frame(
    row.names = unique(kegg_names)  # Use pathway IDs as row names.
  )

  # Add a column for each sample.
  sample_names <- colnames(abundance)[-1]
  for (sample in sample_names) {
    kegg_abundance[[sample]] <- 0
  }

  # Handle errors with tryCatch.
  tryCatch({
    pb <- txtProgressBar(min = 0, max = nrow(kegg_abundance), style = 3)
    start_time <- Sys.time()

    total_matches <- 0

    for (i in seq_len(nrow(kegg_abundance))) {
      setTxtProgressBar(pb, i)
      current_kegg <- rownames(kegg_abundance)[i]

      # Get all KOs associated with the current KEGG pathway.
      relevant_kos <- members[[current_kegg]]

      # Find matching KOs.
      matching_rows <- ko_ids %in% relevant_kos
      if (any(matching_rows)) {
        total_matches <- total_matches + sum(matching_rows)
        kegg_abundance[i, ] <- colSums(abundance[matching_rows, -1, drop = FALSE])
      }
    }

    end_time <- Sys.time()
    close(pb)

    # Add diagnostic information.
    message(sprintf("Total KO matches found: %d", total_matches))
    message(sprintf("Number of non-zero pathways before filtering: %d",
                   sum(rowSums(kegg_abundance) > 0)))

    # Remove pathways with zero abundance.
    message("Removing KEGG pathways with zero abundance across all samples...")
    zero_rows <- rowSums(kegg_abundance) == 0
    kegg_abundance <- kegg_abundance[!zero_rows, , drop = FALSE]

    message(sprintf("Final number of KEGG pathways: %d", nrow(kegg_abundance)))

    if (nrow(kegg_abundance) == 0) {
      stop("No non-zero KEGG pathways found after processing")
    }

    return(kegg_abundance)

  }, error = function(e) {
    if (exists("pb")) close(pb)
    message("Error occurred during processing: ", e$message)
    stop(e)
  })
}

.microfun_ko_ids <- function(x) {
  sub("^ko:", "", trimws(as.character(x)))
}

.microfun_kegg_ids <- function(x) {
  sub("^map", "ko", sub("^path:", "", trimws(as.character(x))))
}

.microfun_kegg_members <- function(kegg_map = NULL) {
  if (is.null(kegg_map)) {
    ref <- as.data.frame(.microfun_chord_reference("ko_to_kegg_reference"))
    kegg_map <- setNames(lapply(seq_len(nrow(ref)), function(i) {
      unlist(ref[i, -1, drop = FALSE], use.names = FALSE)
    }), ref[[1]])
  } else if (is.character(kegg_map) && length(kegg_map) == 1L) {
    lines <- readLines(kegg_map, warn = FALSE)
    fields <- strsplit(lines[nzchar(trimws(lines))], "\t", fixed = TRUE)
    if (any(lengths(fields) < 2L))
      stop("Each KEGG mapping line must have a pathway ID and at least one KO.")
    kegg_map <- setNames(lapply(fields, function(x) x[-1L]), vapply(fields, `[[`, "", 1L))
  }
  if (!is.list(kegg_map) || !length(kegg_map) || is.null(names(kegg_map)))
    stop("kegg_map must be a PICRUSt2 mapping file or a named list of KOs.")
  ids <- .microfun_kegg_ids(names(kegg_map))
  if (any(!grepl("^ko[0-9]{5}$", ids)))
    stop("KEGG mapping requires reference pathway IDs (koNNNNN or mapNNNNN).")
  rows <- split(kegg_map, factor(ids, levels = unique(ids)))
  lapply(rows, function(x) {
    kos <- .microfun_ko_ids(unlist(x, use.names = FALSE))
    unique(kos[grepl("^K[0-9]{5}$", kos)])
  })
}

.microfun_kegg_links <- function(members, pathways) {
  members <- members[intersect(names(members), pathways)]
  data.frame(ko_id = as.character(unlist(members, use.names = FALSE)),
    pathway_id = rep(names(members), lengths(members)), stringsAsFactors = FALSE)
}

.microfun_contribution_column <- function(columns, value_col = NULL, require_abundance = FALSE) {
  if (is.null(value_col)) {
    value_col <- if ("taxon_function_abun" %in% columns || require_abundance)
      "taxon_function_abun" else "norm_taxon_function_contrib"
  }
  if (!value_col %in% columns)
    stop("Missing PICRUSt2 contribution column: ", value_col, ".")
  if (require_abundance && value_col != "taxon_function_abun")
    stop("KEGG aggregation requires PICRUSt2 taxon_function_abun.")
  value_col
}
