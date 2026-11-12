#' @importFrom grDevices col2rgb
#' @importFrom stats as.formula ave model.matrix
#' @importFrom utils flush.console head tail
NULL

# Column names used through data.table, dplyr and ggplot2 non-standard evaluation.
utils::globalVariables(c(
  ".", "..cols", "..cols_to_keep", ".SD", "ASV", "abundance",
  "abundance_log", "angle", "fid_norm", "genus", "hjust", "is_leaf",
  "label", "label_raw", "level", "level_wrapped", "logFC", "node",
  "size_val", "taxon", "x", "y"
))
