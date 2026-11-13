#' Filter significant function-oriented differential results
#'
#' Internal helper that retains rows with non-missing adjusted p-values below a threshold and sorts them by adjusted p-value. It does not add functional annotations.
#'
#' @param daa_result A data frame containing the results of a differential abundance analysis. 
#'   Must include at least the columns \code{feature} and \code{p_adjust}.
#'
#' @param threshold A numeric value specifying the adjusted p-value cutoff for significance. 
#'   Defaults to \code{0.05}.
#'
#' @keywords internal
#'
#' @return A data frame with the input columns, restricted to `p_adjust < threshold` and sorted in ascending adjusted p-value order.
#'
#' @md
result_annotation <- function(daa_result, threshold = 0.05) {
  # Check required columns
  required_cols <- c("feature", "p_adjust")
  if (!all(required_cols %in% colnames(daa_result))) {
    stop("Input must contain at least 'feature' and 'p_adjust' columns.")
  }
  
  # Subset to significant features
  significant <- daa_result[!is.na(daa_result$p_adjust) & daa_result$p_adjust < threshold, ]
  
  # Optional: sort by adjusted p-value
  significant <- significant[order(significant$p_adjust), ]
  
  return(significant)
}

